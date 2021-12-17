import logging
from winterdrp.pipelines.wirc import WircPipeline

logger = logging.getLogger(__name__)

pipelines = {
    "WIRC": WircPipeline
}


def get_pipeline(instrument):

    try:
        pipeline = pipelines[instrument]()
        logger.info(f"Found {instrument} pipeline")
        return pipeline
    except KeyError:
        err = f"Unrecognised pipeline {instrument}. Available pipelines are: {pipelines.keys()}"
        logger.error(err)
        raise KeyError(err)
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
