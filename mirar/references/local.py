"""
Module for getting reference images stored locally
"""
import logging
from pathlib import Path
from typing import Type

from astropy.io import fits

from mirar.data import Image
from mirar.paths import LATEST_WEIGHT_SAVE_KEY
from mirar.processors.sqldatabase.base_model import BaseDB
from mirar.references.base_reference_generator import BaseReferenceGenerator

logger = logging.getLogger(__name__)


class RefFromPath(BaseReferenceGenerator):
    """
    Get locally saved ref from path
    """

    abbreviation = "local"

    def __init__(
        self,
        path: str,
        filter_name: str,
        write_to_db: bool = False,
        db_table: Type[BaseDB] = None,
        duplicate_protocol: str = "replace",
        q3c_bool: bool = True,
    ):
        super().__init__(
            filter_name, write_to_db, db_table, duplicate_protocol, q3c_bool
        )
        self.path = path

    def get_reference(self, image: Image) -> (fits.PrimaryHDU, fits.PrimaryHDU):
        with fits.open(self.path) as hdul:
            ref_hdu = hdul[0].copy()
            weight_path = Path(hdul[0].header[LATEST_WEIGHT_SAVE_KEY])

        ref_weight_hdu = None
        if weight_path.exists():
            with fits.open(weight_path) as hdul:
                ref_weight_hdu = hdul[0].copy()

        return ref_hdu, ref_weight_hdu
