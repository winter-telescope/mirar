import logging
from winterdrp.pipelines.wirc.wirc_pipeline import WircPipeline

logger = logging.getLogger(__name__)

# Convention: lowercase names

pipelines = {
    "wirc": WircPipeline,
}


def get_pipeline(instrument):

    try:
        pipeline = pipelines[instrument.lower()]
        logger.info(f"Found {instrument} pipeline")
    except KeyError:
        err = f"Unrecognised pipeline {instrument}. Available pipelines are: {pipelines.keys()}"
        logger.error(err)
        raise KeyError(err)

    return pipeline()

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
