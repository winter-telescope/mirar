import logging
import os

import astropy.io.fits
import numpy as np
import copy
from winterdrp.paths import saturate_key
from winterdrp.errors import ErrorStack
from winterdrp.processors.base_processor import BaseProcessor
from winterdrp.processors.utils.error_annotator import ErrorStackAnnotator
from winterdrp.data import DataBatch, Dataset, Image
from winterdrp.processors.base_chain import BaseChain

logger = logging.getLogger(__name__)


class Pipeline:
    pipelines = {}

    default_cal_requirements = None
    
    @property
    def name(self):
        raise NotImplementedError()

    @property
    def all_pipeline_configurations(self):
        raise NotImplementedError()

    @property
    def gain(self):
        raise NotImplementedError()

    @property
    def non_linear_level(self):
        raise NotImplementedError()

    def __init__(
            self,
            selected_configurations: str | list[str] = "default",
            night: int | str = "",
    ):

        self.night_sub_dir = os.path.join(self.name, night)
        if not isinstance(selected_configurations, list):
            selected_configurations = [selected_configurations]
        self.selected_configurations = selected_configurations

    @classmethod
    def __init_subclass__(cls, **kwargs):
        super().__init_subclass__(**kwargs)
        if cls.name in cls.pipelines.keys():
            err = f"Pipeline name '{cls.name}' is already found in the pipeline registered keys. " \
                  f"The Pipeline class variable 'name' must be unique!"
            logger.error(err)
            raise ValueError(err)
        cls.pipelines[cls.name] = cls

    def load_pipeline_configuration(
            self,
            configuration: str = "default",
    ) -> BaseChain:
        return copy.copy(self.all_pipeline_configurations[configuration])

    @staticmethod
    def _load_raw_image(path: str) -> tuple[np.ndarray, astropy.io.fits.header]:
        raise NotImplementedError

    def load_raw_image(self, path: str) -> Image:
        data, header = self._load_raw_image(path)
        return Image(data=data, header=header)

    def unpack_raw_image(self, path: str) -> tuple[np.ndarray, astropy.io.fits.Header]:
        return self._load_raw_image(path)

    @staticmethod
    def configure_processors(
            chain: BaseChain,
            sub_dir: str = ""
    ) -> BaseChain:
        for processor in chain.get_child_processors():
            processor.set_night(night_sub_dir=sub_dir)
        return chain

    def add_configuration(
            self,
            configuration_name: str,
            configuration: BaseChain
    ):
        self.all_pipeline_configurations[configuration_name] = configuration

    def set_configuration(
            self,
            new_configuration: str = "default"
    ) -> BaseChain:
        logger.debug(f"Setting pipeline configuration to {new_configuration}.")

        chain = self.load_pipeline_configuration(new_configuration)
        chain = self.configure_processors(chain, sub_dir=self.night_sub_dir)
        processors = chain.get_child_processors()
        for i, (processor) in enumerate(chain.get_child_processors()):
            logger.debug(f"Initialising processor {processor.__class__}")
            processor.set_preceding_steps(previous_steps=processors[:i])
            processor.check_prerequisites()
        logger.debug("Pipeline initialisation complete.")
        return chain

    @staticmethod
    def download_raw_images_for_night(
            night: str | int
    ):
        raise NotImplemented

    def reduce_images(
            self,
            dataset: Dataset,
            output_error_path: str = None,
            catch_all_errors: bool = True,
            selected_configurations: str | list[str] = None
    ) -> tuple[Dataset, ErrorStack]:

        err_stack = ErrorStack()

        if selected_configurations is None:
            selected_configurations = self.selected_configurations

        if not isinstance(selected_configurations, list):
            selected_configurations = [selected_configurations]

        for j, configuration in enumerate(selected_configurations):

            logger.info(f"Using pipeline configuration {configuration} "
                        f"({j+1}/{len(selected_configurations)})")

            chain = self.set_configuration(configuration)

            dataset, new_err_stack = chain.base_apply(
                dataset
            )
            err_stack += new_err_stack

            if np.logical_and(not catch_all_errors, len(err_stack.reports) > 0):
                raise err_stack.reports[0].error

        err_stack.summarise_error_stack(output_path=output_error_path)
        return dataset, err_stack

    def postprocess_configuration(
            self,
            errorstack: ErrorStack,
            selected_configurations: str | list[str],
            processed_images: list[str] = None
    ) -> list[BaseProcessor]:

        cleanup_config = [
            ErrorStackAnnotator(
                errorstack=errorstack,
                processed_images=processed_images),
        ]

        if isinstance(selected_configurations, str):
            cleanup_config += self.all_pipeline_configurations[selected_configurations]
        else:
            for config in selected_configurations:
                cleanup_config += self.all_pipeline_configurations[config]

        return cleanup_config


