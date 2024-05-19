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


class PS1STRM(BaseKowalskiXMatch):
    """
    PS1 STRM Kowalski catalog
    """

    catalog_name = "PS1_STRM"
    abbreviation = "ps1strm"
    projection = {
        "_id": 1,
        "ra": 1,
        "dec": 1,
        "class": 1,
        "z_phot": 1,
        "z_phot_err": 1,
        "prob_QSO": 1,
        "prob_Galaxy": 1,
        "prob_Star": 1,
    }

    column_names = {
        "_id": f"{abbreviation}objid",
        "ra": f"{abbreviation}ra",
        "dec": f"{abbreviation}dec",
        "class": "ps1strmclass",
        "z_phot": "ps1strmzphot",
        "z_phot_err": "ps1strmzphoterr",
        "prob_QSO": "ps1strmprobqso",
        "prob_Galaxy": "ps1strmprobgalaxy",
        "prob_Star": "ps1strmprobstar",
    }

    column_dtypes = {
        "ps1strmobjid": float,
        "ps1strmra": float,
        "ps1strmdec": float,
        "ps1strmclass": str,
        "ps1strmzphot": float,
        "ps1strmzphoterr": float,
        "ps1strmprobqso": float,
        "ps1strmprobgalaxy": float,
        "ps1strmprobstar": float,
    }

    ra_column_name = f"{abbreviation}ra"
    dec_column_name = f"{abbreviation}dec"
