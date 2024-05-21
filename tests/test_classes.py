"""
Module to check that all combinations of pipelines and configurations are valid
"""

import logging
from abc import ABC
from datetime import datetime

from mirar.data.base_data import Dataset
from mirar.errors.exceptions import BaseProcessorError
from mirar.pipelines import Pipeline, get_pipeline
from mirar.processors.base_processor import (
    BaseImageProcessor,
    BaseProcessor,
    BaseSourceGenerator,
    BaseSourceProcessor,
)
from mirar.testing import BaseTestCase

logger = logging.getLogger(__name__)


def all_subclasses(cls):
    """
    Return all subclasses of a class

    :param cls: Parent class
    :return: Set of subclasses
    """
    return set(cls.__subclasses__()).union(
        [s for c in cls.__subclasses__() for s in all_subclasses(c)]
    )


def is_abstract(cls) -> bool:
    """
    Get if ABC is in the object's __bases__ attribute.

    :param cls: Class to check
    :return: Boolean
    """
    try:
        return ABC in cls.__bases__
    except AttributeError:
        return False


class InvalidClassError(BaseProcessorError):
    """
    Class is missing a required method
    """


def test_processor_class(processor: type[BaseProcessor]):
    """
    Test a processor

    :param processor: Processor to test
    :return: None
    """

    # Skip if the processor is abstract
    if is_abstract(processor):
        return

    # Check if the processor has a description
    if processor.description == BaseProcessor.description:
        raise InvalidClassError(f"{processor} is missing a description")

    if processor._apply == BaseProcessor._apply:
        raise InvalidClassError(f"{processor} is missing an _apply method")

    if issubclass(processor, BaseSourceProcessor):
        if processor._apply_to_sources == BaseSourceProcessor._apply_to_sources:
            raise InvalidClassError(
                f"{processor} is missing an _apply_to_sources method"
            )

    elif issubclass(processor, BaseImageProcessor):

        if processor._apply_to_images == BaseImageProcessor._apply_to_images:
            raise InvalidClassError(
                f"{processor} is missing an _apply_to_images method"
            )

    elif issubclass(processor, BaseSourceGenerator):

        if processor._apply_to_images == BaseSourceGenerator._apply_to_images:
            raise InvalidClassError(
                f"{processor} is missing an _apply_to_images method"
            )

    else:
        raise InvalidClassError(
            f"{processor} does not inherit from a valid processor class"
        )


def test_processor_instance(processor: BaseProcessor):
    """
    Test a processor instance

    :param processor: Processor instance
    :return: None
    """
    processor.description()
    processor.base_apply(Dataset())


def test_pipeline(pipeline: Pipeline):
    """
    Test a pipeline

    :param pipeline: Pipeline to test
    :return: None
    """
    try:
        str(pipeline.non_linear_level)
    except NotImplementedError as exc:
        raise InvalidClassError(
            f"{pipeline.__class__.__name__} is missing a non_linear_level method"
        ) from exc

    try:
        all_configs = pipeline.all_pipeline_configurations
    except NotImplementedError as exc:
        raise InvalidClassError(
            f"{pipeline.__class__.__name__} is missing "
            f"an all_pipeline_configurations method"
        ) from exc

    if "default" not in all_configs:
        raise InvalidClassError(
            f"{pipeline.__class__.__name__} is missing a default configuration"
        )


class TestClasses(BaseTestCase):
    """
    Class to test all pipeline configurations
    """

    def setUp(self):
        self.logger = logging.getLogger(__name__)
        self.logger.setLevel(logging.INFO)

    def test_processor_has_functions(self):
        """
        Test the processors in all pipelines

        :return: None
        """

        all_processors = all_subclasses(BaseProcessor)

        self.logger.info("\n\n Testing all processors \n")
        self.logger.info(f"Found {len(all_processors)} processors")

        for processor in all_processors:
            self.logger.info(f"Testing {processor.__name__}")
            test_processor_class(processor)

    def test_all_pipelines(self):
        """
        Test all pipelines

        :return: None
        """

        pipelines = Pipeline.pipelines.keys()

        for pipeline in pipelines:
            pipe = get_pipeline(
                pipeline,
                selected_configurations="default",
                night=str(datetime.now()).split(" ", maxsplit=1)[0].replace("-", ""),
            )

            test_pipeline(pipe)

            config_list = pipe.all_pipeline_configurations.keys()

            for single_config in config_list:
                processor_list = pipe.set_configuration(single_config)

                for processor in processor_list:
                    test_processor_instance(processor)
