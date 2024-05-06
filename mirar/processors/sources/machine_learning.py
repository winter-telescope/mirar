"""
Module with classes to use apply an ML score
"""

from typing import Callable
import logging
import pickle
from pathlib import Path

import pandas as pd
import requests
from mirar.data import SourceBatch
from mirar.paths import ml_models_dir
from mirar.processors.base_processor import BaseSourceProcessor
import tensorflow as tf

logger = logging.getLogger(__name__)


class MLScore(BaseSourceProcessor):
    """
    Class to apply an ML model to a source table
    """

    base_key = "MLScore"

    def __init__(
        self,
        model_url: str,
        apply_to_row: Callable[[object, pd.Series], pd.Series] = None,
    ):
        super().__init__()
        self.model_url = model_url
        self.model_name = Path(self.model_url).name
        self.apply_to_row = apply_to_row

        self._model = None

    def __str__(self) -> str:
        return (
            f"Processor to use ML model {self.model_name} to score sources"
        )

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

        url = self.model_url
        local_path = self.get_ml_path()

        logger.info(
            f"Downloading model {self.model_name} "
            f"from {url} to {local_path}"
        )

        with requests.get(url, stream=True, timeout=120.) as r:
            r.raise_for_status()
            with open(local_path, 'wb') as f:
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
        if not path.exists():
            err = f"Model {path} not found"
            logger.error(err)
            raise FileNotFoundError(err)

        if path.suffix in [".h5", ".hdf5", ".keras"]:
            return tf.keras.models.load_model(path)
        if path.suffix == ".pkl":
            with open(path, "rb") as model_file:
                return pickle.load(model_file)

        raise ValueError(f"Unknown model type {path.suffix}")

    def get_model(self):
        """
        Load the ML model. Download it if it doesn't exist.

        :return: ML model
        """

        if self._model is None:

            local_path = self.get_ml_path()

            if not local_path.exists():
                self.download_model()

            self._model = self.load_model(local_path)

        return self._model

    def _apply_to_sources(
        self,
        batch: SourceBatch,
    ) -> SourceBatch:

        print("Applying ML model to sources")

        model = self.get_model()

        print(model)

        for source_table in batch:

            sources = source_table.get_data()

            new = []

            for _, source in sources.iterrow():
                row = self.apply_to_row(model, source)
                new.append(row)

            new = pd.DataFrame(new)
            source_table.set_data(new)

        return batch


ml = MLScore(model_url='https://github.com//winter-telescope/winter_rb_models/raw/v1.0.0/models/winterdrb_VGG6_20240410_051006.keras') #FIXME: Update URL
ml.apply(SourceBatch())
