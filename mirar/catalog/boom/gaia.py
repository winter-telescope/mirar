"""
Module for querying Gaia using BOOM
"""

from mirar.catalog.boom.base_boom_catalog import BaseBoomXMatch


class Gaia(BaseBoomXMatch):
    """
    Gaia BOOM catalog
    """

    catalog_name = "Gaia_DR3"
    abbreviation = "gaia"
    projection = {
        "_id": 1,
        "ra": 1,
        "dec": 1,
        "parallax": 1,
        "parallax_error": 1,
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

    @staticmethod
    def update_data(data: dict) -> dict:
        """
        Compute parallax_over_error, which Gaia_DR3 doesn't provide.

        :param data: BOOM data
        :return: updated data
        """
        new = {}
        for key, matches in data.items():
            new_matches = []
            for match in matches:
                new_match = dict(match)
                parallax = new_match.get("parallax")
                parallax_error = new_match.get("parallax_error")
                if (
                    parallax is not None
                    and parallax_error is not None
                    and parallax_error != 0
                ):
                    new_match["parallax_over_error"] = parallax / parallax_error
                else:
                    new_match["parallax_over_error"] = None
                new_matches.append(new_match)
            new[key] = new_matches
        return new


class GaiaBright(Gaia):
    """
    Gaia Bright BOOM catalog (Mg < 14)
    """

    abbreviation = "gaiabright"

    @property
    def boom_filter(self) -> dict:
        """
        Filter for BOOM query

        :return: filter
        """
        return {"phot_g_mean_mag": {"$lt": 14}}
