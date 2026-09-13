"""
Module for querying PS1 using BOOM

BOOM merges the separate Kowalski PS1_DR1 (object)/PS1_PSC (star/galaxy
score)/PS1_STRM (photo-z) catalogs into a single "PS1_DR2" catalog, so a
single query (using PS1's own "_id"/ra/dec as the canonical match) can
fetch all three at once, instead of three separate near-sphere queries
that would always match the exact same nearest objects anyway.
"""

from mirar.catalog.boom.base_boom_catalog import BaseBoomXMatch


class PS1(BaseBoomXMatch):
    """
    PS1 BOOM catalog, including star/galaxy score and STRM photo-z
    """

    catalog_name = "PS1_DR2"
    abbreviation = "ps"
    projection = {
        "_id": 1,
        "ra": 1,
        "dec": 1,
        "gMeanPSFMag": 1,
        "rMeanPSFMag": 1,
        "iMeanPSFMag": 1,
        "zMeanPSFMag": 1,
        "ps_score": 1,
        "strm_uid": 1,
        "strm_class": 1,
        "strm_z_phot": 1,
        "strm_z_phot_err": 1,
        "strm_prob_qso": 1,
        "strm_prob_galaxy": 1,
        "strm_prob_star": 1,
    }

    column_names = {
        "_id": "psobjectid",
        "ra": f"{abbreviation}ra",
        "dec": f"{abbreviation}dec",
        "gMeanPSFMag": "sgmag",
        "rMeanPSFMag": "srmag",
        "iMeanPSFMag": "simag",
        "zMeanPSFMag": "szmag",
        "ps_score": "sgscore",
        "strm_uid": "ps1strmobjid",
        "strm_class": "ps1strmclass",
        "strm_z_phot": "ps1strmzphot",
        "strm_z_phot_err": "ps1strmzphoterr",
        "strm_prob_qso": "ps1strmprobqso",
        "strm_prob_galaxy": "ps1strmprobgalaxy",
        "strm_prob_star": "ps1strmprobstar",
    }

    column_dtypes = {
        "psobjectid": float,
        "psra": float,
        "psdec": float,
        "sgmag": float,
        "srmag": float,
        "simag": float,
        "szmag": float,
        "sgscore": float,
        "ps1strmobjid": float,
        "ps1strmclass": str,
        "ps1strmzphot": float,
        "ps1strmzphoterr": float,
        "ps1strmprobqso": float,
        "ps1strmprobgalaxy": float,
        "ps1strmprobstar": float,
    }

    ra_column_name = f"{abbreviation}ra"
    dec_column_name = f"{abbreviation}dec"
