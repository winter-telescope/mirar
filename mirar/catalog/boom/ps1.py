"""
Module for querying PS1 using BOOM
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
        "psobjectid": int,
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

    @property
    def boom_filter(self) -> dict:
        """
        Excludes PS1_DR2 docs missing ps_score - absent from Kowalski's real PS1_DR1.

        :return: filter
        """
        return {"ps_score": {"$exists": True}}

    @staticmethod
    def update_data(data: dict) -> dict:
        """
        Replace PS1's -999 "no measurement" sentinel with None.

        :param data: BOOM data
        :return: updated data
        """
        mag_keys = ("gMeanPSFMag", "rMeanPSFMag", "iMeanPSFMag", "zMeanPSFMag")
        new = {}
        for name, matches in data.items():
            new_matches = []
            for match in matches:
                new_match = dict(match)
                for key in mag_keys:
                    if new_match.get(key) == -999:
                        new_match[key] = None
                new_matches.append(new_match)
            new[name] = new_matches
        return new
