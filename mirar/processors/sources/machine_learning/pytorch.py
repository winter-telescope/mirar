"""
Module with classes to use apply an ML score from pytorch
"""

import logging
from pathlib import Path
from typing import Callable

import pandas as pd
import requests
import torch
from torch import nn

from mirar.data import SourceBatch
from mirar.paths import ml_models_dir
from mirar.processors.base_processor import BaseSourceProcessor

logger = logging.getLogger(__name__)


class Pytorch(BaseSourceProcessor):
    """
    Class to apply a pytorch model to a source table
    """

    base_key = "pytorch"

    def __init__(
        self,
        model: nn.Module,
        model_weights_url: str,
        apply_to_table: Callable[[nn.Module, pd.DataFrame], pd.DataFrame],
    ):
        super().__init__()
        self._model = model
        self.model_weights_url = model_weights_url
        self.model_name = Path(self.model_weights_url).name
        self.apply_to_table = apply_to_table

        self.model = None

    def description(self) -> str:
        return f"Processor to use Pytorch model {self.model_name} to score sources"

    def get_ml_path(self) -> Path:
        """
        Get the path to the ML model

        :return: Path to the ML model
        """
        return ml_models_dir.joinpath(self.model_name)

    def download_model(self):
        """
        Download the ML model
        """

        url = self.model_weights_url
        local_path = self.get_ml_path()

        logger.info(
            f"Downloading model {self.model_name} " f"from {url} to {local_path}"
        )

        with requests.get(url, stream=True, timeout=120.0) as r:
            r.raise_for_status()
            with open(local_path, "wb") as f:
                for chunk in r.iter_content(chunk_size=8192):
                    # If you have chunk encoded response uncomment if
                    # and set chunk_size parameter to None.
                    # if chunk:
                    f.write(chunk)

        if not local_path.exists():
            err = f"Model {self.model_name} not downloaded"
            logger.error(err)
            raise FileNotFoundError(err)

    @staticmethod
    def load_model(path):
        """
        Function to load a pytorch model dict from a path

        :param path: Path to the model
        :return: Pytorch model dict
        """
        if not path.exists():
            err = f"Model {path} not found"
            logger.error(err)
            raise FileNotFoundError(err)

        if path.suffix in [".pth", ".pt"]:
            return torch.load(path)

        raise ValueError(f"Unknown model type {path.suffix}")

    def get_model(self):
        """
        Load the ML model weights. Download it if it doesn't exist.

        :return: ML model
        """

        if self.model is None:

            model = self._model

            local_path = self.get_ml_path()

            if not local_path.exists():
                self.download_model()

            model.load_state_dict(torch.load(local_path))
            model.eval()

            self.model = model

        return self.model

    def _apply_to_sources(
        self,
        batch: SourceBatch,
    ) -> SourceBatch:

        model = self.get_model()

        for source_table in batch:
            sources = source_table.get_data()
            new = self.apply_to_table(model, sources)
            source_table.set_data(new)

        return batch
