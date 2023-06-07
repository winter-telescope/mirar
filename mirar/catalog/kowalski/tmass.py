"""
Module for querying 2MASS using Kowalski
"""
from mirar.catalog.kowalski.base_kowalski_catalog import BaseKowalskiXMatch


class TMASS(BaseKowalskiXMatch):
    """
    2MASS Kowalski catalog
    """

    catalog_name = "2MASS_PSC"
    abbreviation = "tm"
    projection = {
        "_id": 0,
        "designation": 1,
        "ra": 1,
        "decl": 1,
        "j_m": 1,
        "j_msigcom": 1,
        "h_m": 1,
        "h_cmsigcom": 1,
        "k_m": 1,
        "k_cmsigcom": 1,
        "ph_qual": 1,
    }

    column_names = {
        "ra": "tmra",
        "decl": "tmdec",
        "j_m": "tmjmag",
        "h_m": "tmhmag",
        "k_m": "tmkmag",
        "j_msigcom": "tmjmagerr",
        "h_cmsigcom": "tmhmagerr",
        "k_cmsigcom": "tmkmagerr",
        "designation": "tmobjectid",
        "ph_qual": "tmph_qual",
    }

    column_dtypes = {
        "tmra": float,
        "tmdec": float,
        "tmjmag": float,
        "tmhmag": float,
        "tmkmag": float,
        "tmjmagerr": float,
        "tmhmagerr": float,
        "tmkmagerr": float,
        "tmobjectid": str,
        "tmph_qual": str,
    }
