"""
Module for querying Gaia using Kowalski
"""

import astropy.table
import numpy as np
import pandas as pd
from astropy.table import Table

from mirar.catalog.base.base_gaia import BaseGaia2Mass
from mirar.catalog.kowalski.base_kowalski_catalog import BaseKowalskiXMatch


class Gaia2MassKowalski(BaseKowalskiXMatch, BaseGaia2Mass):
    """
    Gaia Kowalski catalog
    """

    catalog_name = "Gaia_DR2_2MASS_best_neighbour"
    abbreviation = "gaia2mass"
    projection = {
        "_id": 1,
        "number_of_neighbours": 1,
        "number_of_mates": 1,
        "best_neighbour_multiplicity": 1,
        "ra": 1,
        "dec": 1,
        "err_maj": 1,
        "err_min": 1,
        "j_m": 1,
        "j_msigcom": 1,
        "ph_qual": 1,
    }

    def __init__(
        self,
        *args,
        **kwargs,
    ):
        super().__init__(*args, **kwargs)
        self.search_radius_arcsec = self.search_radius_arcmin * 60.0

    @property
    def kowalski_filter(self) -> dict:
        """
        Filter for Kowalski query

        :return: filter
        """
        return {
            "number_of_mates": {"$eq": 0},
            "number_of_neighbours": {"$eq": 1},
            f"{self.filter_name.lower()}_m": {"$gt": self.min_mag, "$lt": self.max_mag},
            f"{self.filter_name.lower()}_msigcom": {"$lt": 1.086 / self.snr_threshold},
        }

    @property
    def column_names(self) -> dict:
        return {
            "_id": f"{self.abbreviation}_objectid",
            "number_of_neighbours": f"{self.abbreviation}_number_of_neighbours",
            "number_of_mates": f"{self.abbreviation}_number_of_mates",
            "best_neighbour_multiplicity": f"{self.abbreviation}_best_neighbour_multiplicity",
            "ra": "ra",
            "dec": "dec",
            "err_maj": f"{self.abbreviation}_err_maj",
            "err_min": f"{self.abbreviation}_err_min",
            "j_m": f"{self.abbreviation}_j_m",
            "j_msigcom": f"{self.abbreviation}_j_msigcom",
            "ph_qual": f"{self.abbreviation}_ph_qual",
        }

    @property
    def column_dtypes(self) -> dict:
        return {
            f"{self.abbreviation}_objectid": float,
            f"{self.abbreviation}_number_of_neighbours": float,
            f"{self.abbreviation}_number_of_mates": float,
            f"{self.abbreviation}_best_neighbour_multiplicity": float,
            f"{self.abbreviation}_ra": float,
            f"{self.abbreviation}_dec": float,
            f"{self.abbreviation}_err_maj": float,
            f"{self.abbreviation}_err_min": float,
            f"{self.abbreviation}_j_m": float,
            f"{self.abbreviation}_j_msigcom": float,
            f"{self.abbreviation}_ph_qual": str,
        }

    @property
    def ra_column_name(self) -> str:
        return f"{self.abbreviation}_ra"

    @property
    def dec_column_name(self) -> str:
        return f"{self.abbreviation}_dec"

    def get_source_table(self, ra_deg: float, dec_deg: float) -> astropy.table.Table:
        query_coords = {0: [ra_deg, dec_deg]}

        response = BaseKowalskiXMatch.query(self, coords=query_coords)

        base_list = pd.DataFrame(response["0"])

        err = np.mean([base_list["err_maj"], base_list["err_min"]], axis=0)

        src_list = Table.from_pandas(base_list)
        src_list["ra_errdeg"] = err / 3.6e6
        src_list["dec_errdeg"] = err / 3.6e6
        src_list["FLAGS"] = 0
        src_list["magnitude"] = src_list[f"{self.filter_name.lower()}_m"]
        src_list["magnitude_err"] = src_list[f"{self.filter_name.lower()}_msigcom"]

        return src_list
