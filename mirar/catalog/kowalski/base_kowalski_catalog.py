"""
Module containing the base Kowalski catalog object
"""
import logging
import os
from abc import ABC
from typing import Optional

from penquins import Kowalski

from mirar.catalog.base_catalog import BaseXMatchCatalog
from mirar.errors import ProcessorError

logger = logging.getLogger(__name__)


class KowalskiError(ProcessorError):
    """Error relating to Kowalski"""


PROTOCOL, HOST, PORT = "https", "kowalski.caltech.edu", 443


def get_kowalski() -> Kowalski:
    """
    Get a Kowalski object, using credentials stored in the environment

    :return: Kowalski object
    """

    token_kowalski = os.environ.get("KOWALSKI_TOKEN")

    if token_kowalski is not None:
        logger.debug("Using kowalski token")

        kowalski_instance = Kowalski(
            token=token_kowalski, protocol=PROTOCOL, host=HOST, port=PORT
        )

    else:
        username_kowalski = os.environ.get("KOWALSKI_USER")
        password_kowalski = os.environ.get("KOWALSKI_PWD")

        if username_kowalski is None:
            err = (
                "Kowalski username not provided, "
                "please run export KOWALSKI_USER=<user>"
            )
            logger.error(err)
            raise KowalskiError(err)

        if password_kowalski is None:
            err = (
                "Kowalski password not provided, "
                "please run export KOWALSKI_PWD=<user>"
            )
            logger.error(err)
            raise KowalskiError(err)

        kowalski_instance = Kowalski(
            username=username_kowalski,
            password=password_kowalski,
            protocol=PROTOCOL,
            host=HOST,
            port=PORT,
        )

    if not kowalski_instance.ping():
        err = "Error connecting to Kowalski. Are your credentials right?"
        logger.error(err)
        raise KowalskiError(err)

    return kowalski_instance


class BaseKowalskiXMatch(BaseXMatchCatalog, ABC):
    """
    Base class for a catalog using Kowalski
    """

    def __init__(
        self,
        *args,
        kowalski: Optional[Kowalski] = None,
        max_time_ms: float = 10000,
        **kwargs,
    ):
        super().__init__(*args, **kwargs)
        self.max_time_ms = max_time_ms
        self.kowalski = kowalski

    def near_query_kowalski(self, coords: dict) -> dict:
        """
        Performs a Kowalski query around coords

        :param coords: ra/dec
        :return: crossmatch dict
        """
        query = {
            "query_type": "near",
            "query": {
                "max_distance": self.search_radius_arcsec,
                "distance_units": "arcsec",
                "radec": coords,
                "catalogs": {
                    f"{self.catalog_name}": {
                        "filter": {},
                        "projection": self.projection,
                    }
                },
            },
            "kwargs": {
                "max_time_ms": self.max_time_ms,
                "limit": self.num_sources,
            },
        }
        logger.debug(f"Kowalski is {self.kowalski}")
        response = self.kowalski.query(query=query)
        data = response.get("default").get("data")

        return data[self.catalog_name]

    def query(self, coords) -> dict:
        """
        Uses a Kowalski object to query for sources around coords

        :param coords: ra/dec
        :return: crossmatch sources
        """
        if self.kowalski is None:
            self.kowalski = get_kowalski()
        logger.debug("Querying kowalski")
        data = self.near_query_kowalski(coords)
        return data
