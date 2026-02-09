# from mirar.data import Image
from mirar.pipelines.spring.config import sextractor_astrometry_config


def spring_anet_sextractor_config_path_generator() -> str:
    """
    Generate the path to the ANET SExtractor configuration file for SPRING images
    Parameters
    ----------
    image

    Returns
    -------

    """
    return sextractor_astrometry_config["config_path"]
