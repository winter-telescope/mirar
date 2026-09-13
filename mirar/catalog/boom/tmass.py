"""
Module for querying 2MASS using BOOM
"""

from mirar.catalog.base.base_gaia import offsets_2mass
from mirar.catalog.boom.base_boom_catalog import BaseBoomXMatch


class TMASS(BaseBoomXMatch):
    """
    2MASS BOOM catalog
    """

    catalog_name = "2MASS_PSC"
    abbreviation = "tm"
    projection = {
        "_id": 1,
        "ra": 1,
        "dec": 1,
        "j_m": 1,
        "j_msigcom": 1,
        "h_m": 1,
        "h_msigcom": 1,
        "k_m": 1,
        "k_msigcom": 1,
        "ph_qual": 1,
    }

    column_names = {
        "ra": f"{abbreviation}ra",
        "dec": f"{abbreviation}dec",
        "j_m": "tmjmag",
        "h_m": "tmhmag",
        "k_m": "tmkmag",
        "j_msigcom": "tmjmagerr",
        "h_msigcom": "tmhmagerr",
        "k_msigcom": "tmkmagerr",
        "_id": "tmobjectid",
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

    ra_column_name = f"{abbreviation}ra"
    dec_column_name = f"{abbreviation}dec"

    @staticmethod
    def update_data(data: dict) -> dict:
        """
        For a given catalog, update the data with any extra information

        :param data: BOOM data
        :return: updated data
        """

        new = {}

        for key, val in data.items():
            new_list = []

            for vega_dict in val:
                ab_dict = dict(vega_dict)

                for tmass_filter, offset in offsets_2mass.items():
                    tmass_key = f"{tmass_filter}_m"
                    vega_mag = ab_dict[tmass_key]
                    if vega_mag is not None:
                        ab_dict[tmass_key] = vega_mag + offset

                new_list.append(ab_dict)

            new[key] = new_list

        return data
