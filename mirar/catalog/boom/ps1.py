"""
Module for querying PS1 using BOOM

BOOM merges the separate Kowalski PS1_DR1 (object)/PS1_PSC (star/galaxy
score)/PS1_STRM (photo-z) catalogs into a single "PS1_DR2" catalog. The
three classes below therefore all query the same BOOM catalog, each with
its own projection, to keep the same external interface (and downstream
column names) as before.
"""

from mirar.catalog.boom.base_boom_catalog import BaseBoomXMatch


class PS1(BaseBoomXMatch):
    """
    PS1 BOOM catalog
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
    }

    column_names = {
        "_id": "psobjectid",
        "ra": f"{abbreviation}ra",
        "dec": f"{abbreviation}dec",
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


class PS1SGSc(BaseBoomXMatch):
    """
    PS1 Star/Galaxy Score BOOM catalog
    """

    catalog_name = "PS1_DR2"
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


class PS1STRM(BaseBoomXMatch):
    """
    PS1 STRM BOOM catalog
    """

    catalog_name = "PS1_DR2"
    abbreviation = "ps1strm"
    projection = {
        # BOOM always returns "_id" unless explicitly excluded; this class
        # uses "strm_uid" (the STRM-specific id) instead, so "_id" must be
        # excluded or it ends up in results with no column_names entry.
        "_id": 0,
        "strm_uid": 1,
        "ra": 1,
        "dec": 1,
        "strm_class": 1,
        "strm_z_phot": 1,
        "strm_z_phot_err": 1,
        "strm_prob_qso": 1,
        "strm_prob_galaxy": 1,
        "strm_prob_star": 1,
    }

    column_names = {
        "strm_uid": f"{abbreviation}objid",
        "ra": f"{abbreviation}ra",
        "dec": f"{abbreviation}dec",
        "strm_class": "ps1strmclass",
        "strm_z_phot": "ps1strmzphot",
        "strm_z_phot_err": "ps1strmzphoterr",
        "strm_prob_qso": "ps1strmprobqso",
        "strm_prob_galaxy": "ps1strmprobgalaxy",
        "strm_prob_star": "ps1strmprobstar",
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
