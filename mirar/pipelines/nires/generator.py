from mirar.catalog.vizier.ukidssgps import UkidssGPS


def nires_astrometric_ref_catalog_generator(_) -> UkidssGPS:
    """
    Generates the astrometric reference catalog for the NIRES image
    """
    return UkidssGPS(min_mag=10, max_mag=22, filter_name="K", search_radius_arcmin=3)
