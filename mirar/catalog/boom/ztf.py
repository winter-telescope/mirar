"""
Module for querying ZTF using BOOM
"""

from mirar.catalog.boom.base_boom_catalog import BaseBoomXMatch


class ZTF(BaseBoomXMatch):
    """
    ZTF BOOM catalog
    """

    catalog_name = "ZTF_alerts"
    abbreviation = "ztf"
    projection = {
        "_id": 1,
        "objectId": 1,
        "candidate.ra": 1,
        "candidate.dec": 1,
        "candidate.rb": 1,
        "candidate.drb": 1,
    }

    column_names = {
        "_id": "ztfid",
        "objectId": "ztfname",
        "candidate.ra": f"{abbreviation}ra",
        "candidate.dec": f"{abbreviation}dec",
        "candidate.rb": f"{abbreviation}rb",
        "candidate.drb": f"{abbreviation}drb",
    }

    column_dtypes = {
        "ztfid": str,
        "ztfname": str,
        "ztfra": float,
        "ztfdec": float,
        "ztfrb": float,
        "ztfdrb": float,
    }

    ra_column_name = f"{abbreviation}ra"
    dec_column_name = f"{abbreviation}dec"
