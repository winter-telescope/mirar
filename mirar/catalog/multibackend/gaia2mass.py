"""
Composite catalog for Gaia 2Mass
"""

import logging
from typing import Type

from mirar.catalog.base.base_catalog import BaseCatalog, BaseMultiBackendCatalog
from mirar.catalog.kowalski import get_kowalski
from mirar.catalog.kowalski.gaia2mass import Gaia2MassKowalski
from mirar.catalog.tap.gaia2mass import Gaia, Gaia2MassTAP
from mirar.catalog.vizier.gaia2mass import Gaia2MassVizier

logger = logging.getLogger(__name__)


class Gaia2Mass(BaseMultiBackendCatalog):
    """
    Composite catalog for Gaia 2Mass
    """

    abbreviation = "tmass"

    @staticmethod
    def set_backend(backend: str | None) -> Type[BaseCatalog]:

        if backend is None:
            k = get_kowalski()
            if k.ping():
                backend = "kowalski"

        if backend is None:
            # pylint: disable=protected-access,no-member
            if Gaia._TapPlus__getconnhandler().get_response_status() == 200:
                # Gaia goes down sometimes
                # Response status 0 means it's down, 200 when up and working
                backend = "gaia_tap"
            else:
                logger.warning("Gaia TAP service is down, cannot use default backend.")

        logger.debug(f"Backend for Gaia2Mass: {backend}")

        if backend == "kowalski":
            return Gaia2MassKowalski
        if backend == "gaia_tap":
            return Gaia2MassTAP
        if backend == "vizier":
            return Gaia2MassVizier

        raise NotImplementedError(f"Backend '{backend}' not implemented for Gaia2Mass")
