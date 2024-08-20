"""
Module for querying ZTF using Kowalski
"""

from mirar.catalog.kowalski.base_kowalski_catalog import BaseKowalskiXMatch


class ZTF(BaseKowalskiXMatch):
    """
    ZTF Kowalski catalog
    """

    catalog_name = "ZTF_alerts"
    abbreviation = "ztf"
    projection = {
        "objectId": 1,
        "ra": 1,
        "dec": 1,
    }

    column_names = {
        "objectId": "ztfname",
        "ra": f"{abbreviation}ra",
        "dec": f"{abbreviation}dec",
    }

    column_dtypes = {
        "ztfname": str,
        "ztfra": float,
        "ztfdec": float,
    }

    ra_column_name = f"{abbreviation}ra"
    dec_column_name = f"{abbreviation}dec"
