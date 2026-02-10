# from mirar.data import Image
from mirar.data import Image
from mirar.pipelines.spring.config import sextractor_astrometry_config


def spring_anet_sextractor_config_path_generator(_image: Image) -> str:
    """
    Generate the path to the ANET SExtractor configuration file for SPRING images
    Parameters
    ----------
    image:Image

    Returns
    -------
    sextractor_config_path:str

    """
    return sextractor_astrometry_config["config_path"]
