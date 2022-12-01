import logging

from winterdrp.errors import ProcessorError
from winterdrp.pipelines.base_pipeline import Pipeline
from winterdrp.pipelines.summer.summer_pipeline import SummerPipeline
from winterdrp.pipelines.wirc.wirc_pipeline import WircPipeline

logger = logging.getLogger(__name__)


# Convention: lowercase names


class PipelineConfigError(ProcessorError, KeyError):
    pass


def get_pipeline(instrument, selected_configurations=None, *args, **kwargs):

    try:
        pipeline = Pipeline.pipelines[instrument.lower()]
        logger.info(f"Found {instrument} pipeline")
    except KeyError:
        err = f"Unrecognised pipeline {instrument}. Available pipelines are: {Pipeline.pipelines.keys()}"
        logger.error(err)
        raise PipelineConfigError(err)

    return pipeline(selected_configurations=selected_configurations, *args, **kwargs)
