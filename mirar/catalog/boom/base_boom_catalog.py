"""
Module containing the base BOOM catalog object
"""

import logging
from abc import ABC
from typing import Optional

import numpy as np
from astropy import units as u
from astropy.coordinates import SkyCoord
from blastwave.errors import BOOMCredentialsError
from blastwave.query.boom import BoomClient

from mirar.catalog.base.base_xmatch_catalog import BaseXMatchCatalog
from mirar.errors import ProcessorError

logger = logging.getLogger(__name__)


class BOOMError(ProcessorError):
    """Error relating to BOOM"""


def get_boom_client() -> BoomClient:
    """
    Get a BoomClient object, using credentials stored in the environment
    (BOOM_API_USER / BOOM_API_PASSWORD, or a .env file)

    :return: BoomClient object
    """
    boom_instance = BoomClient()

    try:
        boom_instance.get_session_headers()
    except BOOMCredentialsError as exc:
        logger.error(str(exc))
        raise BOOMError(str(exc)) from exc

    res = boom_instance.ping()
    if not res.ok:
        err = "Error connecting to BOOM. Are your credentials right?"
        logger.error(err)
        raise BOOMError(err)

    return boom_instance


def flatten_boom_data(matches: list[dict]) -> list[dict]:
    """
    Flatten a BOOM data dict

    :param matches: List of matches
    :return: Flattened list of depth-1 dictionaries
    """

    new = []

    if len(matches) > 0:
        for match in matches:
            new_dict = {}
            if isinstance(match, dict):
                for key, val in match.items():
                    if isinstance(val, dict):
                        for subkey, subval in val.items():
                            new_dict[f"{key}.{subkey}"] = subval
                    else:
                        new_dict[key] = val
            new.append(new_dict)

    return new


class BaseBoomXMatch(BaseXMatchCatalog, ABC):
    """
    Base class for a catalog using BOOM
    """

    @property
    def boom_filter(self) -> Optional[dict]:
        """
        Filter for BOOM query

        :return: filter
        """
        return None

    def __init__(
        self,
        *args,
        boom: Optional[BoomClient] = None,
        max_time_ms: float = 10000,
        **kwargs,
    ):
        super().__init__(*args, **kwargs)
        self.max_time_ms = max_time_ms
        self.boom = boom

    def near_query_boom(self, coords: dict) -> dict:
        """
        Performs a batched BOOM cone search around each coordinate

        BOOM's cone_search does not sort matches by distance before
        applying `limit` - confirmed empirically that a `limit=3` query
        can return matches over 4x further away than others it excluded
        within the same radius. `limit` is omitted here (BOOM then
        returns every match within `radius`, which is itself small - a
        crossmatch search radius, not a survey area), and matches are
        sorted by true angular separation and truncated to num_sources
        ourselves instead.

        :param coords: dict of {name: [ra, dec]}
        :return: dict of {name: matches}
        """
        payload = {
            "catalog_name": self.catalog_name,
            "object_coordinates": coords,
            "radius": self.search_radius_arcsec,
            "unit": "Arcseconds",
            "projection": self.projection,
            "max_time_ms": int(self.max_time_ms),
        }
        if self.boom_filter is not None:
            payload["filter"] = self.boom_filter

        logger.debug(f"BOOM client is {self.boom}")
        res = self.boom.api("post", "queries/cone_search", data=payload)

        if not res.ok:
            err = f"BOOM cone search on '{self.catalog_name}' failed: {res.text}"
            logger.error(err)
            res.raise_for_status()

        data = res.json()["data"]
        flattened = {key: flatten_boom_data(matches) for key, matches in data.items()}
        return self._sort_by_separation(flattened, coords)

    def _sort_by_separation(self, flattened: dict, coords: dict) -> dict:
        """
        Sort each coordinate's matches by true angular separation and
        truncate to num_sources.

        :param flattened: dict of {name: matches}, as returned by BOOM
        :param coords: dict of {name: [ra, dec]}, the queried coordinates
        :return: dict of {name: matches}, nearest-first, truncated
        """
        raw_ra_key = next(
            key for key, val in self.column_names.items() if val == self.ra_column_name
        )
        raw_dec_key = next(
            key for key, val in self.column_names.items() if val == self.dec_column_name
        )

        sorted_data = {}
        for name, matches in flattened.items():
            if len(matches) == 0:
                sorted_data[name] = matches
                continue

            query_coord = SkyCoord(
                ra=coords[name][0] * u.deg, dec=coords[name][1] * u.deg
            )
            match_coords = SkyCoord(
                ra=[match[raw_ra_key] for match in matches] * u.deg,
                dec=[match[raw_dec_key] for match in matches] * u.deg,
            )
            order = np.argsort(query_coord.separation(match_coords))
            sorted_data[name] = [matches[i] for i in order[: self.num_sources]]

        return sorted_data

    def query(self, coords) -> dict:
        """
        Uses a BOOM client to query for sources around coords

        :param coords: ra/dec
        :return: crossmatch sources
        """
        if self.boom is None:
            self.boom = get_boom_client()
        logger.debug("Querying BOOM")
        data = self.near_query_boom(coords)
        data = self.update_data(data)
        return data

    @staticmethod
    def update_data(data: dict) -> dict:
        """
        For a given catalog, update the data with any extra information

        :param data: BOOM data
        :return: updated data
        """
        return data
