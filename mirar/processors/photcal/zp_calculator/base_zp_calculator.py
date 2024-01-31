"""
Module for calculating zero points
"""

import logging

from astropy.table import Table

from mirar.data import Image

logger = logging.getLogger(__name__)


class BaseZeroPointCalculator:
    """
    Base class for actually calculating zero point calculators
    """

    def calculate_zeropoint(
        self,
        image: Image,
        matched_ref_cat: Table,
        matched_img_cat: Table,
        colnames: list[str],
    ) -> Image:
        """
        Function to calculate zero point from two catalogs and add the zeropoint
        information to the image header
        Args:
            matched_ref_cat: Reference catalog table
            matched_img_cat: Catalog of sources
            image: Image object
            colnames: List of column names from the image catalog to use for
            calculating zero point. The reference catalog is assumed to have the
            "magnitude" column, as it comes from mirar.catalog.base_catalog.BaseCatalog
        Returns:
            Image object with zero point information added to the header.
        """
        raise NotImplementedError
