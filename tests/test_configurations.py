"""
Module to check that all combinations of pipelines and configurations are valid
"""

import logging
from datetime import datetime

from mirar.pipelines import Pipeline, get_pipeline
from mirar.processors.base_processor import PrerequisiteError
from mirar.testing import BaseTestCase

logger = logging.getLogger(__name__)


def check_configurations():
    """
    Check each configuration in each pipeline, to ensure that it is valid.
    Will fail if a
    :class:`~mirar.processors.base_processor.PrerequisiteError` is raised.

    :return: None
    """
    pipelines = Pipeline.pipelines.keys()

    for pipeline in pipelines:
        pipe = get_pipeline(
            pipeline,
            selected_configurations="default",
            night=str(datetime.now()).split(" ", maxsplit=1)[0].replace("-", ""),
        )

        logger.info(f"Visualising {pipeline} pipeline")

        config_list = pipe.all_pipeline_configurations.keys()

        for single_config in config_list:
            logger.info(f"Visualising {single_config} configuration")
            try:
                pipe.set_configuration(single_config)
            except PrerequisiteError as exc:
                err = (
                    f"Error for '{pipeline}' pipeline, "
                    f"'{single_config}' configuration"
                )
                logger.error(err)
                raise PrerequisiteError(err) from exc


class TestPipelineConfigurations(BaseTestCase):
    """
    Class to test all pipeline configurations
    """

    def setUp(self):
        self.logger = logging.getLogger(__name__)
        self.logger.setLevel(logging.INFO)

    def test_pipelines(self):
        """
        Test the pipeline configurations

        :return: None
        """
        self.logger.info("\n\n Testing pipeline configurations \n\n")
        check_configurations()
