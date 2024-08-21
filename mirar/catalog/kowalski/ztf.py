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
        "_id": 1,
        "objectId": 1,
        "candidate.ra": 1,
        "candidate.dec": 1,
    }

    column_names = {
        "_id": "ztfid",
        "objectId": "ztfname",
        "candidate.ra": f"{abbreviation}ra",
        "candidate.dec": f"{abbreviation}dec",
    }

    column_dtypes = {
        "ztfid": str,
        "ztfname": str,
        "ztfra": float,
        "ztfdec": float,
    }

    ra_column_name = f"{abbreviation}ra"
    dec_column_name = f"{abbreviation}dec"
