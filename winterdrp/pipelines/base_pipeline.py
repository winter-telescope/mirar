import logging
import os

import astropy.io.fits
import numpy as np
import copy
from winterdrp.paths import saturate_key
from winterdrp.errors import ErrorStack
from winterdrp.processors.base_processor import BaseProcessor

logger = logging.getLogger(__name__)

core_fields = ["OBSCLASS", "TARGET", "UTCTIME"]


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
    ):
        return copy.copy(self.all_pipeline_configurations[configuration])

    @staticmethod
    def load_raw_image(path: str) -> tuple[np.ndarray, astropy.io.fits.header]:
        raise NotImplementedError

    @staticmethod
    def configure_processors(
            processors: list[BaseProcessor],
            sub_dir: str = ""
    ) -> list[BaseProcessor]:

        for processor in processors:
            processor.set_night(night_sub_dir=sub_dir)
        return processors

    def add_configuration(
            self,
            configuration_name: str,
            configuration: str | list[BaseProcessor]
    ):
        self.all_pipeline_configurations[configuration_name] = configuration

    def set_configuration(
            self,
            new_configuration: str = "default"
    ) -> list[BaseProcessor]:
        logger.debug(f"Setting pipeline configuration to {new_configuration}.")

        processors = self.load_pipeline_configuration(new_configuration)
        processors = self.configure_processors(processors, sub_dir=self.night_sub_dir)
        for i, (processor) in enumerate(processors):
            logger.debug(f"Initialising processor {processor.__class__}")
            processor.set_preceding_steps(previous_steps=processors[:i])
            processor.check_prerequisites()
        logger.debug("Pipeline initialisation complete.")
        return processors

    @staticmethod
    def download_raw_images_for_night(
            night: str | int
    ):
        raise NotImplemented

    def reduce_images(
            self,
            batches: list[list[list[np.ndarray], list[astropy.io.fits.header]]],
            output_error_path: str = None,
            catch_all_errors: bool = True,
            selected_configurations: str | list[str] = None
    ):
        err_stack = ErrorStack()

        if selected_configurations is None:
            selected_configurations = self.selected_configurations

        if not isinstance(selected_configurations, list):
            selected_configurations = [selected_configurations]

        for j, configuration in enumerate(selected_configurations):

            logger.info(f"Using pipeline configuration {configuration} "
                        f"({j+1}/{len(selected_configurations)})")

            processors = self.set_configuration(configuration)

            for i, processor in enumerate(processors):
                logger.debug(f"Applying '{processor.__class__}' processor to {len(batches)} batches. "
                             f"(Step {i+1}/{len(processors)})")

                batches, new_err_stack = processor.base_apply(
                    batches
                )
                err_stack += new_err_stack

                if np.logical_and(not catch_all_errors, len(err_stack.reports) > 0):
                    raise err_stack.reports[0].error

        err_stack.summarise_error_stack(output_path=output_error_path)
        return batches, err_stack

    def set_saturation(
            self,
            header: astropy.io.fits.Header
    ) -> astropy.io.fits.Header:
        # update the SATURATE keyword in the header for subsequent sextractor runs
        co_add_head = header['COADDS']
        num_co_adds = int(co_add_head)
        saturation_level = self.non_linear_level * num_co_adds
        if "SKMEDSUB" in header.keys():
            saturation_level -= header['SKMEDSUB']
        header.append((saturate_key, saturation_level, 'Saturation level'), end=True)
        return header