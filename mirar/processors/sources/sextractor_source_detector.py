"""
Module to detect candidates in an image
"""

import logging
from pathlib import Path

import astropy
import numpy as np
import pandas as pd
from astropy import units as u
from astropy.coordinates import SkyCoord, match_coordinates_sky

from mirar.data import Image, ImageBatch, SourceBatch, SourceTable
from mirar.paths import (
    BASE_NAME_KEY,
    CAND_DEC_KEY,
    CAND_RA_KEY,
    SEXTRACTOR_HEADER_KEY,
    XPOS_KEY,
    YPOS_KEY,
    get_output_dir,
)
from mirar.processors.astromatic import Sextractor
from mirar.processors.base_processor import BaseSourceGenerator, PrerequisiteError
from mirar.utils.ldac_tools import get_table_from_ldac

logger = logging.getLogger(__name__)


def generate_candidates_table(
    image: Image,
    sextractor_catalog_path: str | Path,
    target_only: bool = True,
) -> pd.DataFrame:
    """
    Generate a candidates table from a sextractor catalog
    :param sextractor_catalog_path: Path to the sextractor catalog
    :return: Candidates table
    """
    det_srcs = get_table_from_ldac(sextractor_catalog_path)
    logger.debug(f"Found {len(det_srcs)} sources in image.")

    multi_col_mask = [det_srcs.dtype[i].shape != () for i in range(len(det_srcs.dtype))]
    if np.sum(multi_col_mask) != 0:
        logger.warning(
            "Sextractor catalog contains multi-dimensional columns, "
            "removing them before converting to pandas dataframe"
        )
        to_remove = np.array(det_srcs.colnames)[multi_col_mask]
        det_srcs.remove_columns(to_remove)

    if target_only:
        logger.debug(
            f"Isolating target from sextractor catalog, "
            f"at position {image[CAND_RA_KEY]}, {image[CAND_DEC_KEY]}."
        )
        det_srcs = isolate_target(image, det_srcs)

    det_srcs = det_srcs.to_pandas()

    ydims, xdims = image.get_data().shape
    det_srcs["NAXIS1"] = xdims
    det_srcs["NAXIS2"] = ydims
    det_srcs[XPOS_KEY] = det_srcs["X_IMAGE"] - 1
    det_srcs[YPOS_KEY] = det_srcs["Y_IMAGE"] - 1
    det_srcs[CAND_RA_KEY] = det_srcs["ALPHAWIN_J2000"]
    det_srcs[CAND_DEC_KEY] = det_srcs["DELTAWIN_J2000"]
    det_srcs["fwhm"] = det_srcs["FWHM_IMAGE"]
    det_srcs["elong"] = det_srcs["ELONGATION"]

    return det_srcs


def isolate_target(
    image: Image,
    sextractor_catalog: astropy.table.Table,
) -> astropy.table.Table:
    """
    Args:
        image: Image object containing the target source coordinates
        sextractor_catalog: sextractor catalog as an astropy Table
    Returns: Table with len=1, the target source from the sextractor catalog
    """
    cat_coords = SkyCoord(
        ra=sextractor_catalog["ALPHAWIN_J2000"],
        dec=sextractor_catalog["DELTAWIN_J2000"],
        unit=(u.deg, u.deg),
    )
    targ_coords = SkyCoord(
        ra=image[CAND_RA_KEY], dec=image[CAND_DEC_KEY], unit=(u.deg, u.deg)
    )

    idx, _, __ = match_coordinates_sky(targ_coords, cat_coords)
    matched_sextractor_catalog = sextractor_catalog[idx]

    logger.debug(
        f"Found nearest neighbor source at "
        f"{matched_sextractor_catalog['ALPHAWIN_J2000']}, "
        f"{matched_sextractor_catalog['DELTAWIN_J2000']}"
    )
    return astropy.table.Table(matched_sextractor_catalog)


class SextractorSourceDetector(BaseSourceGenerator):
    """
    Class that retrieves a sextractor catalog and saves all sources to a sourcetable
    """

    base_key = "DETSOURC"

    def __init__(
        self,
        output_sub_dir: str = "sources",
        target_only: bool = True,
    ):
        """
        :param output_sub_dir: subdirectory to output files
        :param target_only: whether the returned sourcetable should contain the
        target source only (vs. all sources in the image)
        """
        super().__init__()
        self.output_sub_dir = output_sub_dir
        self.target_only = target_only

    def __str__(self) -> str:
        return "Retrieves a sextractor catalog and converts it to a sourcetable"

    def get_sub_output_dir(self) -> Path:
        """
        Returns: output sub-directory
        """
        return get_output_dir(self.output_sub_dir, self.night_sub_dir)

    def _apply_to_images(
        self,
        batch: ImageBatch,
    ) -> SourceBatch:
        all_sources = SourceBatch()
        for image in batch:
            srcs_table = generate_candidates_table(
                image=image,
                sextractor_catalog_path=image[SEXTRACTOR_HEADER_KEY],
                target_only=self.target_only,
            )

            if len(srcs_table) > 0:
                x_shape, y_shape = image.get_data().shape
                srcs_table["X_SHAPE"] = x_shape
                srcs_table["Y_SHAPE"] = y_shape

            metadata = self.get_metadata(image)

            if len(srcs_table) == 0:
                msg = f"No sources found in image {image[BASE_NAME_KEY]}"
                logger.warning(msg)

            else:
                msg = f"Found {len(srcs_table)} sources in image {image[BASE_NAME_KEY]}"
                logger.debug(msg)
                all_sources.append(SourceTable(srcs_table, metadata=metadata))

        return all_sources

    def check_prerequisites(self):
        check = np.sum([isinstance(x, Sextractor) for x in self.preceding_steps])
        if check == 0:
            raise PrerequisiteError(
                "Sextractor must be run before SextractorSourceDetector"
            )
