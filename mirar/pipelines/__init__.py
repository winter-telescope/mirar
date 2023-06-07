import logging

from mirar.errors import ProcessorError
from mirar.pipelines.base_pipeline import Pipeline
from mirar.pipelines.sedmv2.sedmv2_pipeline import SEDMv2Pipeline
from mirar.pipelines.summer.summer_pipeline import SummerPipeline
from mirar.pipelines.winter.winter_pipeline import WINTERPipeline
from mirar.pipelines.wirc.wirc_pipeline import WircPipeline

logger = logging.getLogger(__name__)


# Convention: lowercase names


class PipelineConfigError(ProcessorError, KeyError):
    pass


def get_pipeline(instrument, selected_configurations=None, *args, **kwargs):
    try:
        pipeline = Pipeline.pipelines[instrument.lower()]
        logger.info(f"Found {instrument} pipeline")
    except KeyError:
        err = (
            f"Unrecognised pipeline {instrument}. "
            f"Available pipelines are: {Pipeline.pipelines.keys()}"
        )
        logger.error(err)
        raise PipelineConfigError(err)

    return pipeline(selected_configurations=selected_configurations, *args, **kwargs)
