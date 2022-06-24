import logging
from winterdrp.pipelines.base_pipeline import Pipeline
from winterdrp.pipelines.wirc.wirc_pipeline import WircPipeline
#from winterdrp.pipelines.summer.summer_pipeline import SummerPipeline
from winterdrp.pipelines.wirc_imsub.wirc_imsub_pipeline import WircImsubPipeline

logger = logging.getLogger(__name__)

# Convention: lowercase names


def get_pipeline(instrument, configuration, *args, **kwargs):

    try:
        pipeline = Pipeline.pipelines[instrument.lower()]
        logger.info(f"Found {instrument} pipeline")
    except KeyError:
        err = f"Unrecognised pipeline {instrument}. Available pipelines are: {Pipeline.pipelines.keys()}"
        logger.error(err)
        raise KeyError(err)

    return pipeline(pipeline_configuration=configuration, *args, **kwargs)

# def parse_telescope(header):
#
#     telname = None
#
#     try:
#         if header['DETNAM'] == 'Hawaii-II':
#             if header['FPA'] == 'WIRC':
#                 telname = "WIRC"
#     except KeyError:
#         pass
#
#     if telname is None:
#         msg = f"Unable to recognise instrument from fits header."
#         logger.warning(msg)
#         logger.debug(header)
#     else:
#         logger.debug(f'Telescope recognised as {telname}')
#     return telname
#
#
# update_functions = {
#     "WIRC": update_wirc,
#     None: lambda x: x
# }
#
#
# def reformat_raw_data(img):
#
#     header = img[0].header
#
#     return update_functions[parse_telescope(header)](img)
