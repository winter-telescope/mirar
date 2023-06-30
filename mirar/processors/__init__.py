"""
This module contains the processors that are used to process the raw data
"""
# import logging
from mirar.processors.base_processor import BaseImageProcessor
from mirar.processors.bias import BiasCalibrator
from mirar.processors.dark import DarkCalibrator
from mirar.processors.flat import FlatCalibrator, SkyFlatCalibrator
from mirar.processors.mask import MaskPixelsFromPath
from mirar.processors.utils.image_saver import ImageSaver

#
# logger = logging.getLogger(__name__)
#
#
# def get_processor(processor_name, open_fits, *args, **kwargs):
#
#     try:
#         processor = BaseProcessor.subclasses[processor_name]
#     except KeyError:
#         err = f"Processor type '{processor_name}' not recognised. " \
#               f"The following processors are available:
#               {BaseProcessor.subclasses.keys()}"
#         logger.error(err)
#         raise KeyError(err)
#
#     return processor(open_fits, *args, **kwargs)
