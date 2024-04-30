"""
Central location for all pipelines. This is where you should add new pipelines.
"""

import importlib
import logging

from mirar.errors import ProcessorError
from mirar.pipelines.base_pipeline import Pipeline

logger = logging.getLogger(__name__)


# Convention: lowercase names


class PipelineConfigError(ProcessorError, KeyError):
    """
    Error raised when a pipeline is not found
    """


def get_pipeline(
    instrument: str, *args, selected_configurations=None, **kwargs
) -> Pipeline:
    """
    Function to get pipeline

    :param instrument: Name of instrument
    :param args: args
    :param selected_configurations: Configurations to use
    :param kwargs: kwargs
    :return: pipeline
    """
    try:
        importlib.import_module(
            f"mirar.pipelines.{instrument.lower()}.{instrument.lower()}_pipeline"
        )
        pipeline = Pipeline.pipelines[instrument.lower()]
        logger.info(f"Found {instrument} pipeline")
    except KeyError as exc:
        err = (
            f"Unrecognised pipeline {instrument}."
            f"Available pipelines are: {Pipeline.pipelines.keys()}"
        )
        logger.error(err)
        raise PipelineConfigError(err) from exc
    except ModuleNotFoundError as exc:
        err = f"Unrecognised pipeline {instrument}."
        logger.error(err)
        raise PipelineConfigError(err) from exc

    return pipeline(selected_configurations=selected_configurations, *args, **kwargs)
