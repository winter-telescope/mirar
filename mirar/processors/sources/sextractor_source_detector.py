"""
Module to detect candidates in an image
"""

import logging
import os
from pathlib import Path

import pandas as pd

from mirar.data import Image, ImageBatch, SourceBatch, SourceTable
from mirar.paths import (
    BASE_NAME_KEY,
    CAND_DEC_KEY,
    CAND_RA_KEY,
    LATEST_WEIGHT_SAVE_KEY,
    XPOS_KEY,
    YPOS_KEY,
    get_output_dir,
)
from mirar.processors.astromatic.sextractor.sourceextractor import run_sextractor_single
from mirar.processors.base_processor import BaseSourceGenerator
from mirar.utils.ldac_tools import get_table_from_ldac

logger = logging.getLogger(__name__)


def generate_candidates_table(
    image: Image,
    sex_catalog_path: str | Path,
) -> pd.DataFrame:
    """
    Generate a candidates table from a sextractor catalog
    :param sex_catalog_path: Path to the sextractor catalog
    :return: Candidates table
    """
    det_srcs = get_table_from_ldac(sex_catalog_path)
    det_srcs = det_srcs.to_pandas()

    logger.debug(f"Found {len(det_srcs)} sources in image.")

    ydims, xdims = image.get_data().shape
    det_srcs["NAXIS1"] = xdims
    det_srcs["NAXIS2"] = ydims
    det_srcs[XPOS_KEY] = det_srcs["X_IMAGE"] - 1
    det_srcs[YPOS_KEY] = det_srcs["Y_IMAGE"] - 1
    xpeaks, ypeaks = det_srcs["XPEAK_IMAGE"] - 1, det_srcs["YPEAK_IMAGE"] - 1
    det_srcs["xpeak"] = xpeaks
    det_srcs["ypeak"] = ypeaks
    det_srcs[CAND_RA_KEY] = det_srcs["ALPHA_J2000"]
    det_srcs[CAND_DEC_KEY] = det_srcs["DELTA_J2000"]
    det_srcs["fwhm"] = det_srcs["FWHM_IMAGE"]
    det_srcs["aimage"] = det_srcs["A_IMAGE"]
    det_srcs["bimage"] = det_srcs["B_IMAGE"]
    det_srcs["aimagerat"] = det_srcs["aimage"] / det_srcs["fwhm"]
    det_srcs["bimagerat"] = det_srcs["bimage"] / det_srcs["fwhm"]
    det_srcs["elong"] = det_srcs["ELONGATION"]

    return det_srcs


class SextractorSourceDetector(BaseSourceGenerator):
    """
    Class to detect sources by running sourceextractor on an image,
    then saves them all to a sourcetable
    """

    base_key = "DETSOURC"

    def __init__(  # pylint: disable=too-many-arguments
        self,
        cand_det_sextractor_config: str,
        cand_det_sextractor_filter: str,
        cand_det_sextractor_nnw: str,
        cand_det_sextractor_params: str,
        output_sub_dir: str = "sources",
    ):
        super().__init__()
        self.output_sub_dir = output_sub_dir
        self.cand_det_sextractor_config = cand_det_sextractor_config
        self.cand_det_sextractor_filter = cand_det_sextractor_filter
        self.cand_det_sextractor_nnw = cand_det_sextractor_nnw
        self.cand_det_sextractor_params = cand_det_sextractor_params

    def __str__(self) -> str:
        return (
            "Extracts detected sources from images, "
            "and converts them to a pandas dataframe"
        )

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
            image_path = os.path.join(self.get_sub_output_dir(), image[BASE_NAME_KEY])
            weight = image[LATEST_WEIGHT_SAVE_KEY]

            cands_catalog_name = image_path.replace(".fits", ".dets")
            cands_catalog_name, _ = run_sextractor_single(
                img=image_path,
                output_dir=self.get_sub_output_dir(),
                catalog_name=cands_catalog_name,
                config=self.cand_det_sextractor_config,
                parameters_name=self.cand_det_sextractor_params,
                filter_name=self.cand_det_sextractor_filter,
                starnnw_name=self.cand_det_sextractor_nnw,
                weight_image=weight,
                gain=1.0,
            )

            srcs_table = generate_candidates_table(
                image=image,
                sex_catalog_path=cands_catalog_name,
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
