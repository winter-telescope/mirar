from mirar.catalog.vizier.ukidssgps import UkidssGPS
from mirar.pipelines.nires.config import sextractor_anet_config


def nires_sextractor_config_path_generator(_) -> str:
    """
    Generates the sextractor config file path for the NIRES image
    """
    return sextractor_anet_config["config_path_boardid_0_2_3_4"]


def nires_astrometric_ref_catalog_generator(_) -> UkidssGPS:
    """
    Generates the astrometric reference catalog for the NIRES image
    """
    return UkidssGPS(min_mag=10, max_mag=22, filter_name="K", search_radius_arcmin=3)
