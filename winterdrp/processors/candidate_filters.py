import logging

import astropy.table
import pandas as pd

from winterdrp.processors.base_processor import BaseDataframeProcessor
import numpy as np
from astropy.io import fits
from collections.abc import Callable
from winterdrp.catalog.base_catalog import BaseCatalog
from pandas import DataFrame
from winterdrp.processors.astromatic.sextractor.sextractor import Sextractor
from winterdrp.processors.astromatic.sextractor.sourceextractor import run_sextractor_dual
from winterdrp.utils.ldac_tools import get_table_from_ldac
from winterdrp.paths import get_output_dir
import io
import gzip
import os

logger = logging.getLogger(__name__)


class EdgeCandidatesMask(BaseDataframeProcessor):

    def __init__(self,
                 edge_boundary_size: float,
                 image_xsize_column_key: str = 'X_SHAPE',
                 image_ysize_column_key: str = 'Y_SHAPE',
                 x_column_key: str = 'X_IMAGE',
                 y_column_key: str = 'Y_IMAGE',
                 *args,
                 **kwargs):
        super(EdgeCandidatesMask, self).__init__(*args, **kwargs)

        self.edge_boundary_size = edge_boundary_size
        self.image_xsize_column_key = image_xsize_column_key
        self.image_ysize_column_key = image_ysize_column_key
        self.x_column_key = x_column_key
        self.y_column_key = y_column_key

    def _apply_to_images(
            self,
            tables: list[pd.DataFrame]
    ) -> list[pd.DataFrame]:
        all_masked_tables = []

        for table in tables:
            x_coords = table[self.x_column_key]
            y_coords = table[self.y_column_key]
            image_xsize = table[self.image_xsize_column_key]
            image_ysize = table[self.image_ysize_column_key]

            logger.info(f'Applying edge-filter to {len(table)} candidates')
            near_x_edge_mask = (x_coords < self.edge_boundary_size) | (
                    x_coords > image_xsize - self.edge_boundary_size)
            near_y_edge_mask = (y_coords < self.edge_boundary_size) | (
                    y_coords > image_ysize - self.edge_boundary_size)

            near_edge_mask = near_y_edge_mask | near_x_edge_mask
            masked_table = table[np.invert(near_edge_mask)]

            logger.info(f'Edge-filter passed {len(masked_table)} candidates')
            all_masked_tables.append(masked_table)

        return all_masked_tables


class BrightStarCandidateMask(BaseDataframeProcessor):

    def __init__(self,
                 ref_catalog_generator: Callable[[astropy.io.fits.Header], BaseCatalog],
                 bright_star_thresh_mag: float = 14,
                 *args,
                 **kwargs):
        super(BrightStarCandidateMask, self).__init__(*args, **kwargs)
        self.ref_catalog_generator = ref_catalog_generator
        self.bright_star_thresh_mag = bright_star_thresh_mag

    def _apply_to_images(
            self,
            tables: list[DataFrame]
    ) -> list[DataFrame]:
        pass
