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
        "raMean": f"{abbreviation}ra",
        "decMean": f"{abbreviation}dec",
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

    ra_column_name = f"{abbreviation}ra"
    dec_column_name = f"{abbreviation}dec"


class PS1SGSc(BaseKowalskiXMatch):
    """
    PS1 Star/Galaxy Score Kowalski catalog
    """

    catalog_name = "PS1_PSC"
    abbreviation = "sgscore"
    projection = {
        "_id": 1,
        "ra": 1,
        "dec": 1,
        "ps_score": 1,
    }

    column_names = {
        "_id": f"{abbreviation}objid",
        "ra": f"{abbreviation}ra",
        "dec": f"{abbreviation}dec",
        "ps_score": "sgscore",
    }

    column_dtypes = {
        "sgscoreobjid": float,
        "sgscorera": float,
        "sgscoredec": float,
        "sgscore": float,
    }

    ra_column_name = f"{abbreviation}ra"
    dec_column_name = f"{abbreviation}dec"
