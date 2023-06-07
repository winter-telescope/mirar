"""
Module for querying PS1 using Kowalski
"""
from mirar.catalog.kowalski.base_kowalski_catalog import BaseKowalskiXMatch


class PS1(BaseKowalskiXMatch):
    """
    PS1 Kowalski catalog
    """

    catalog_name = "PS1_DR1"
    abbreviation = "ps"
    projection = {
        "_id": 1,
        "raMean": 1,
        "decMean": 1,
        "gMeanPSFMag": 1,
        "rMeanPSFMag": 1,
        "iMeanPSFMag": 1,
        "zMeanPSFMag": 1,
    }

    column_names = {
        "_id": "psobjectid",
        "raMean": "psra",
        "decMean": "psdec",
        "gMeanPSFMag": "sgmag",
        "rMeanPSFMag": "srmag",
        "iMeanPSFMag": "simag",
        "zMeanPSFMag": "szmag",
    }

    column_dtypes = {
        "psobjectid": float,
        "psra": float,
        "psdec": float,
        "sgmag": float,
        "srmag": float,
        "simag": float,
        "szmag": float,
    }
