"""
Composite catalog for Gaia 2Mass
"""

import logging
from typing import Type

from mirar.catalog.base.base_catalog import BaseCatalog, BaseMultiBackendCatalog
from mirar.catalog.tap.gaia2mass import Gaia, Gaia2MassARI, Gaia2MassTAP, gaia_ari
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
            backend = "vizier"

        if backend is None:

            # Check server is alive
            cmd = (
                "SELECT table_name from tap_schema.tables "
                "WHERE table_name = 'gaiadr3.gaia_source'"
            )
            job = gaia_ari.launch_job(cmd, dump_to_file=False)
            job.get_results()

            # pylint: disable=protected-access,no-member
            if gaia_ari._TapPlus__getconnhandler().get_response_status() == 200:
                # Gaia ARI also goes down sometimes
                # Response status 0 means it's down, 200 when up and working
                backend = "gaia_ari"

        if backend is None:
            # pylint: disable=protected-access,no-member
            if Gaia._TapPlus__getconnhandler().get_response_status() == 200:
                # Gaia goes down sometimes
                # Response status 0 means it's down, 200 when up and working
                backend = "gaia_tap"

        logger.debug(f"Backend for Gaia2Mass: {backend}")

        if backend == "gaia_ari":
            return Gaia2MassARI
        if backend == "vizier":
            return Gaia2MassVizier
        if backend == "gaia_tap":
            return Gaia2MassTAP

        raise NotImplementedError(f"Backend '{backend}' not implemented for Gaia2Mass")
