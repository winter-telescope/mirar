"""
Module containing the base BOOM catalog object
"""

import logging
from abc import ABC
from typing import Optional

from blastwave.errors import BOOMCredentialsError
from blastwave.query.boom import BOOMClient

from mirar.catalog.base.base_xmatch_catalog import BaseXMatchCatalog
from mirar.errors import ProcessorError

logger = logging.getLogger(__name__)


class BOOMError(ProcessorError):
    """Error relating to BOOM"""


def get_boom_client() -> BOOMClient:
    """
    Get a BOOMClient object, using credentials stored in the environment
    (BOOM_API_USER / BOOM_API_PASSWORD, or a .env file)

    :return: BOOMClient object
    """
    boom_instance = BOOMClient()

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
        boom: Optional[BOOMClient] = None,
        **kwargs,
    ):
        super().__init__(*args, **kwargs)
        self.boom = boom

    def near_query_boom(self, coords: dict) -> dict:
        """
        Performs a BOOM cone search around each coordinate using
        BOOMClient.cone_search().

        :param coords: dict of {name: [ra, dec]}
        :return: dict of {name: matches}
        """
        results = {}
        for name, (ra, dec) in coords.items():
            raw_matches = self.boom.cone_search(
                ra=ra,
                dec=dec,
                radius_arcsec=self.search_radius_arcsec,
                catalog=self.catalog_name,
                limit=self.num_sources,
                filter_query=self.boom_filter,
                projection=self.projection,
            )
            results[name] = flatten_boom_data(raw_matches)

        return results

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
