from mirar.data import Image
from mirar.pipelines.cfht.config import sextractor_astrometry_config


def cfht_anet_sextractor_config_path_generator(image: Image) -> str:
    """
    Generates the sextractor config file path for the winter image
    """
    return sextractor_astrometry_config["config_path"]
