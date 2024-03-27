"""
Module for querying Gaia using Kowalski
"""

from mirar.catalog.kowalski.base_kowalski_catalog import BaseKowalskiXMatch


class Gaia(BaseKowalskiXMatch):
    """
    Gaia Kowalski catalog
    """

    catalog_name = "Gaia_EDR3"
    abbreviation = "gaia"
    projection = {
        "_id": 1,
        "ra": 1,
        "dec": 1,
        "parallax": 1,
        "parallax_error": 1,
        "parallax_over_error": 1,
        "ruwe": 1,
    }

    @property
    def column_names(self) -> dict:
        return {
            "_id": f"{self.abbreviation}_objectid",
            "ra": f"{self.abbreviation}_ra",
            "dec": f"{self.abbreviation}_dec",
            "parallax": f"{self.abbreviation}_parallax",
            "parallax_error": f"{self.abbreviation}_parallax_error",
            "parallax_over_error": f"{self.abbreviation}_parallax_over_error",
            "ruwe": f"{self.abbreviation}_ruwe",
        }

    @property
    def column_dtypes(self) -> dict:
        return {
            f"{self.abbreviation}_objectid": float,
            f"{self.abbreviation}_ra": float,
            f"{self.abbreviation}_dec": float,
            f"{self.abbreviation}_parallax": float,
            f"{self.abbreviation}_parallax_error": float,
            f"{self.abbreviation}_parallax_over_error": float,
            f"{self.abbreviation}_ruwe": float,
        }

    @property
    def ra_column_name(self) -> str:
        return f"{self.abbreviation}_ra"

    @property
    def dec_column_name(self) -> str:
        return f"{self.abbreviation}_dec"


class GaiaBright(Gaia):
    """
    Gaia Bright Kowalski catalog (Mg < 14)
    """

    abbreviation = "gaiabright"

    @property
    def kowalski_filter(self) -> dict:
        """
        Filter for Kowalski query

        :return: filter
        """
        return {"phot_g_mean_mag": {"$lt": 14}}
