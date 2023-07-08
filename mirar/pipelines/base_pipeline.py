"""
Module containing the base of the :class:`~mirar.pipelines.base_pipeline.Pipeline`
 class.

Each :class:`~mirar.pipelines.base_pipeline.Pipeline` will have several
`configurations`. A configuration corresponds to a list of
:class:`~mirar.processors.BaseProcessor` objects.

The pipeline will process data using a chosen list of these individual
:class:`~mirar.processors.BaseProcessor` objects.
"""
import copy
import logging
import os
from pathlib import Path
from typing import Optional

import astropy.io.fits
import numpy as np

from mirar.data import Dataset, Image, ImageBatch
from mirar.errors import ErrorStack
from mirar.paths import get_output_path
from mirar.processors.base_processor import BaseProcessor
from mirar.processors.utils.error_annotator import ErrorStackAnnotator

logger = logging.getLogger(__name__)


class Pipeline:
    """
    Base class for pipelines.

    Each pipeline must have the following class variables:
     * a name (the name of the instrument
     * pipeline configurations
     * gain
     * a :func:`~mirar.pipelines.base_pipeline.Pipeline._load_raw_image` function
     to load raw images and modify the headers etc as required
    """

    pipelines = {}

    default_cal_requirements = None

    @property
    def name(self):
        """
        Unique name of pipeline , used to call it from the command via
        :func:`~mirar.pipelines.get_pipeline`.
        Should be the name of the instrument.
        """
        raise NotImplementedError()

    @property
    def all_pipeline_configurations(self):
        """Dictionary containing all pipeline configurations"""
        raise NotImplementedError()

    @property
    def gain(self):
        """Gain of instrument"""
        raise NotImplementedError()

    @property
    def non_linear_level(self):
        """Non-linear level of instrument"""
        raise NotImplementedError()

    def __init__(
        self,
        selected_configurations: str | list[str] = "default",
        night: int | str = "",
    ):
        self.night_sub_dir = os.path.join(self.name, night)
        self.night = night
        if not isinstance(selected_configurations, list):
            selected_configurations = [selected_configurations]
        self.selected_configurations = selected_configurations

    @classmethod
    def __init_subclass__(cls, **kwargs):
        super().__init_subclass__(**kwargs)
        if cls.name in cls.pipelines:
            err = (
                f"Pipeline name '{cls.name}' is already found "
                f"in the pipeline registered keys. "
                f"The Pipeline class variable 'name' must be unique!"
            )
            logger.error(err)
            raise ValueError(err)
        cls.pipelines[cls.name] = cls

    def load_pipeline_configuration(
        self,
        configuration: str = "default",
    ) -> list[BaseProcessor]:
        """
        Load a particular named configuration from self.all_pipeline_configurations

        :param configuration: configuration to be used
        :return: list of processors
        """
        return copy.copy(self.all_pipeline_configurations[configuration])

    @staticmethod
    def _load_raw_image(path: str | Path) -> tuple[np.ndarray, astropy.io.fits.header]:
        """
        Function to load in a raw image and ensure it
        has the correct format for the code.

        :param path: path of raw image
        :return: The image data and image header
        """
        raise NotImplementedError

    def load_raw_image(self, path: str) -> Image:
        """
        Function to load in a raw image and create an
        :class:`~mirar.data.image_data.Image` object which
        can then be processed further

        :param path: path of raw image
        :return: Image object
        """
        data, header = self._load_raw_image(path)
        return Image(data=data, header=header)

    def unpack_raw_image(self, path: str) -> tuple[np.ndarray, astropy.io.fits.Header]:
        """
        Function to load in a raw image and ensure it has
         the correct format for the code.
        This function, unlike
        :func:`~mirar.pipelines.base_pipeline.Pipeline._load_raw_image`
        is not protected.

        :param path: path of raw image
        :return: The image data and image header
        """
        return self._load_raw_image(path)

    @staticmethod
    def configure_processors(
        processors: list[BaseProcessor], sub_dir: str = ""
    ) -> list[BaseProcessor]:
        """
        Propagates the correct nightly setting to a list of processors.

        :param processors: Processors to configure
        :param sub_dir: night sub directory to use
        :return: Updated processors
        """
        for processor in processors:
            processor.set_night(night_sub_dir=sub_dir)
        return processors

    def add_configuration(
        self, configuration_name: str, configuration: list[BaseProcessor]
    ):
        """
        Add a new configuration to the pipeline.

        :param configuration_name: Name of new configuration
        :param configuration: the list of processors
        :return: None
        """
        self.all_pipeline_configurations[configuration_name] = configuration

    def set_configuration(
        self, new_configuration: str = "default"
    ) -> list[BaseProcessor]:
        """
        Loads a new configuration for the pipeline

        :param new_configuration: name of configuration
        :return: list of corresponding processors
        """
        logger.debug(f"Setting pipeline configuration to {new_configuration}.")

        processors = self.load_pipeline_configuration(new_configuration)
        processors = self.configure_processors(processors, sub_dir=self.night_sub_dir)
        for i, processor in enumerate(processors):
            logger.debug(f"Initialising processor {processor.__class__}")
            processor.set_preceding_steps(previous_steps=processors[:i])
            processor.check_prerequisites()
        logger.debug("Pipeline initialisation complete.")
        return processors

    @staticmethod
    def download_raw_images_for_night(night: str | int):
        """
        Function to download images from a remote server

        :param night: Night of data to download
        :return: None
        """
        raise NotImplementedError()

    def get_error_output_path(self) -> Path:
        """
        Generates a unique path for the error summary,
        in the output data directory.
        Makes the parent directory structure if needed.

        :return: path for error summary
        """
        error_output_path = Path(
            get_output_path(
                base_name=f"{self.night}_error_stack.txt",
                dir_root=self.night_sub_dir,
            )
        )

        if not error_output_path.parent.exists():
            error_output_path.parent.mkdir(parents=True)

        return error_output_path

    def reduce_images(
        self,
        dataset: Optional[Dataset] = None,
        output_error_path: Optional[str] = None,
        catch_all_errors: bool = True,
        selected_configurations: Optional[str | list[str]] = None,
    ) -> tuple[Dataset, ErrorStack]:
        """
        Function to process a given dataset.

        :param dataset: dataset to process (can  be empty)
        :param output_error_path: optional path to write error summary
        :param catch_all_errors: Either catch errors, or just immediately raise them
        :param selected_configurations: Configuration to use
        :return: Post-processing dataset and summary of errors caught
        """

        if dataset is None:
            dataset = Dataset([ImageBatch()])

        if output_error_path is None:
            output_error_path = self.get_error_output_path()

        err_stack = ErrorStack()

        if selected_configurations is None:
            selected_configurations = self.selected_configurations

        if not isinstance(selected_configurations, list):
            selected_configurations = [selected_configurations]

        for j, configuration in enumerate(selected_configurations):
            logger.info(
                f"Using pipeline configuration {configuration} "
                f"({j+1}/{len(selected_configurations)})"
            )

            processors = self.set_configuration(configuration)

            for i, processor in enumerate(processors):
                logger.info(
                    f"Applying '{processor.__class__} to {len(dataset)} batches "
                    f"(Step {i + 1}/{len(processors)})"
                )
                logger.info(f"[{str(processor)}]")

                dataset, new_err_stack = processor.base_apply(dataset)
                err_stack += new_err_stack

                if np.logical_and(not catch_all_errors, len(err_stack.reports) > 0):
                    raise err_stack.reports[0].error

                if len(dataset) == 0:
                    logger.error(
                        f"No images left in dataset. "
                        f"Terminating early, after step {i + 1}/{len(processors)} "
                        f"({processor.__class__.__name__})."
                    )
                    break

        err_stack.summarise_error_stack(output_path=output_error_path)
        return dataset, err_stack

    def postprocess_configuration(
        self,
        errorstack: ErrorStack,
        selected_configurations: str | list[str],
        processed_images: Optional[list[str]] = None,
    ) -> list[BaseProcessor]:
        """
        Generate a postprocessing/cleanup processor sequence,
        Used by  :class:`~mirar.monitor.base_monitor.Monitor` class
        for realtime processing.
        The first step is to update the header of images
        with any saved errors in errorstack.

        :param errorstack: Caught errors
        :param selected_configurations: Configurations to use.
        :param processed_images: list of processed images
        :return: list of postprocess processors
        """

        cleanup_config = [
            ErrorStackAnnotator(
                errorstack=errorstack, processed_images=processed_images
            ),
        ]

        if isinstance(selected_configurations, str):
            cleanup_config += self.all_pipeline_configurations[selected_configurations]
        else:
            for config in selected_configurations:
                cleanup_config += self.all_pipeline_configurations[config]

        return cleanup_config
