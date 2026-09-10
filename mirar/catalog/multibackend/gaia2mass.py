"""
Composite catalog for Gaia 2Mass
"""

import logging
from typing import Type

from mirar.catalog.base.base_catalog import BaseCatalog, BaseMultiBackendCatalog
from mirar.catalog.tap.gaia2mass import Gaia2MassARI, Gaia2MassTAP
from mirar.catalog.vizier.gaia2mass import Gaia2MassVizier

logger = logging.getLogger(__name__)

DEFAULT_GAIA2MASS_BACKEND = "vizier"


class Gaia2Mass(BaseMultiBackendCatalog):
    """
    Composite catalog for Gaia 2Mass
    """

    abbreviation = "tmass"

    @staticmethod
    def set_backend(backend: str | None) -> Type[BaseCatalog]:

        if backend is None:
            backend = DEFAULT_GAIA2MASS_BACKEND

        logger.debug(f"Backend for Gaia2Mass: {backend}")

        if backend == "gaia_ari":
            return Gaia2MassARI
        if backend == "vizier":
            return Gaia2MassVizier
        if backend == "gaia_tap":
            return Gaia2MassTAP

        raise NotImplementedError(f"Backend '{backend}' not implemented for Gaia2Mass")
