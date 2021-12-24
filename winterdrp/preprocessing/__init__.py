import logging
from winterdrp.preprocessing.bias import BiasCalibrator
from winterdrp.preprocessing.dark import DarkCalibrator
from winterdrp.preprocessing.flat import FlatCalibrator

logger = logging.getLogger(__name__)

processor_map = {
    "bias": BiasCalibrator,
    "dark": DarkCalibrator,
    "flat": FlatCalibrator
}


def get_processor(processor_name, open_fits, *args, **kwargs):

    try:
        processor = processor_map[processor_name]
    except KeyError:
        err = f"Processor type '{processor_name}' not recognised. " \
              f"The following processors are available: {processor_map.keys()}"
        logger.error(err)
        raise KeyError(err)

    return processor(open_fits, *args, **kwargs)
