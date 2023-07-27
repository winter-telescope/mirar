"""
Module to detect candidates in an image
"""

import logging
import os
from pathlib import Path

import numpy as np
import pandas as pd
from astropy.io import fits

from mirar.data import Image, ImageBatch, SourceBatch, SourceTable
from mirar.data.utils import encode_img
from mirar.paths import (
    BASE_NAME_KEY,
    CAND_DEC_KEY,
    CAND_RA_KEY,
    LATEST_SAVE_KEY,
    NORM_PSFEX_KEY,
    REF_IMG_KEY,
    UNC_IMG_KEY,
    XPOS_KEY,
    YPOS_KEY,
    ZP_KEY,
    get_output_dir,
)
from mirar.processors.astromatic.sextractor.sourceextractor import run_sextractor_dual
from mirar.processors.base_processor import BaseSourceGenerator
from mirar.processors.photometry.utils import make_cutouts
from mirar.utils.ldac_tools import get_table_from_ldac

logger = logging.getLogger(__name__)


class SourceDetector(BaseSourceGenerator):
    """
    Class to detect candidates by running sourceextractor on an image
    """

    base_key = "DETCANDS"

    def __init__(
        self,
        cand_det_sextractor_config: str,
        cand_det_sextractor_filter: str,
        cand_det_sextractor_nnw: str,
        cand_det_sextractor_params: str,
        output_sub_dir: str = "candidates",
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
        Get output sub-directory
        Returns:

        """
        return get_output_dir(self.output_sub_dir, self.night_sub_dir)

    def generate_candidates_table(
        self,
        diff: Image,
        scorr_catalog_path: str | Path,
        sci_resamp_image_path: str | Path,
        ref_resamp_image_path: str | Path,
        diff_scorr_path: str | Path,
        diff_psf_path: str | Path,
        diff_unc_path: str | Path,
    ) -> pd.DataFrame:
        """
        Generate a candidates table from a difference image

        :param scorr_catalog_path: Path to the scorr catalog
        :param sci_resamp_image_path:
        :param ref_resamp_image_path:
        :param diff_path:
        :param diff_scorr_path:
        :param diff_psf_path:
        :param diff_unc_path:
        :return:
        """
        det_srcs = get_table_from_ldac(scorr_catalog_path)
        det_srcs = det_srcs.to_pandas()

        diff_path = diff[LATEST_SAVE_KEY]

        if len(det_srcs) == 0:
            return pd.DataFrame()
        logger.debug(f"Found {len(det_srcs)} candidates in image.")
        det_srcs[XPOS_KEY] = det_srcs["X_IMAGE"] - 1
        det_srcs[YPOS_KEY] = det_srcs["Y_IMAGE"] - 1

        scorr_data = fits.getdata(diff_scorr_path)
        xpeaks, ypeaks = det_srcs["XPEAK_IMAGE"] - 1, det_srcs["YPEAK_IMAGE"] - 1
        scorr_peaks = scorr_data[ypeaks, xpeaks]
        det_srcs["xpeak"] = xpeaks
        det_srcs["ypeak"] = ypeaks
        det_srcs["scorr"] = scorr_peaks

        cutout_size_display = 40

        display_sci_ims = []
        display_ref_ims = []
        display_diff_ims = []

        for _, row in det_srcs.iterrows():
            xpeak, ypeak = int(row["xpeak"]), int(row["ypeak"])

            display_sci_cutout, display_ref_cutout, display_diff_cutout = make_cutouts(
                [sci_resamp_image_path, ref_resamp_image_path, diff_path],
                (xpeak, ypeak),
                cutout_size_display,
            )

            display_sci_bit = encode_img(display_sci_cutout)
            display_ref_bit = encode_img(display_ref_cutout)
            display_diff_bit = encode_img(display_diff_cutout)

            display_sci_ims.append(display_sci_bit)
            display_ref_ims.append(display_ref_bit)
            display_diff_ims.append(display_diff_bit)

        det_srcs["cutoutScience"] = display_sci_ims
        det_srcs["cutoutTemplate"] = display_ref_ims
        det_srcs["cutoutDifference"] = display_diff_ims

        det_srcs[ZP_KEY] = diff[ZP_KEY]
        det_srcs[LATEST_SAVE_KEY] = diff[LATEST_SAVE_KEY]
        det_srcs["magzpsci"] = diff[ZP_KEY]
        diff_zp_unc = diff[f"{ZP_KEY}_STD"]
        det_srcs["magzpsciunc"] = diff_zp_unc
        det_srcs["diffimname"] = diff[LATEST_SAVE_KEY]
        det_srcs["sciimname"] = sci_resamp_image_path
        det_srcs["refimname"] = ref_resamp_image_path
        det_srcs[NORM_PSFEX_KEY] = diff_psf_path
        det_srcs[UNC_IMG_KEY] = diff_unc_path
        det_srcs[CAND_RA_KEY] = det_srcs["ALPHA_J2000"]
        det_srcs[CAND_DEC_KEY] = det_srcs["DELTA_J2000"]
        det_srcs["fwhm"] = det_srcs["FWHM_IMAGE"]
        det_srcs["aimage"] = det_srcs["A_IMAGE"]
        det_srcs["bimage"] = det_srcs["B_IMAGE"]
        det_srcs["aimagerat"] = det_srcs["aimage"] / det_srcs["fwhm"]
        det_srcs["bimagerat"] = det_srcs["bimage"] / det_srcs["fwhm"]
        det_srcs["elong"] = det_srcs["ELONGATION"]

        det_srcs["jd"] = fits.getval(sci_resamp_image_path, "MJD-OBS") + 2400000.5
        det_srcs["exptime"] = diff["EXPTIME"]
        det_srcs["field"] = fits.getval(sci_resamp_image_path, "FIELDID")
        det_srcs["programpi"] = fits.getval(sci_resamp_image_path, "PROGPI")
        det_srcs["programid"] = fits.getval(sci_resamp_image_path, "PROGID")
        det_srcs["fid"] = fits.getval(sci_resamp_image_path, "FILTERID")
        det_srcs["candid"] = np.array(
            det_srcs["jd"] * 100, dtype=int
        ) * 10000 + np.arange(len(det_srcs))
        det_srcs["diffmaglim"] = diff["DIFFMLIM"]
        det_srcs["isdiffpos"] = 1

        logger.debug(det_srcs["diffmaglim"])

        return det_srcs

    def _apply_to_images(
        self,
        batch: ImageBatch,
    ) -> SourceBatch:
        all_cands = SourceBatch()
        for image in batch:
            scorr_image_path = os.path.join(self.get_sub_output_dir(), image["DIFFSCR"])
            diff_image_path = os.path.join(self.get_sub_output_dir(), image["DIFFIMG"])
            diff_psf_path = os.path.join(self.get_sub_output_dir(), image["DIFFPSF"])
            diff_unc_path = os.path.join(self.get_sub_output_dir(), image["DIFFUNC"])

            scorr_mask_path = os.path.join(self.get_sub_output_dir(), image["SCORMASK"])
            cands_catalog_name = diff_image_path.replace(".fits", ".dets")
            cands_catalog_name, _ = run_sextractor_dual(
                det_image=scorr_image_path,
                measure_image=diff_image_path,
                output_dir=self.get_sub_output_dir(),
                catalog_name=cands_catalog_name,
                config=self.cand_det_sextractor_config,
                parameters_name=self.cand_det_sextractor_params,
                filter_name=self.cand_det_sextractor_filter,
                starnnw_name=self.cand_det_sextractor_nnw,
                weight_image=scorr_mask_path,
                gain=1.0,
            )

            sci_image_path = os.path.join(
                self.get_sub_output_dir(), image[BASE_NAME_KEY]
            )
            ref_image_path = os.path.join(self.get_sub_output_dir(), image[REF_IMG_KEY])
            srcs_table = self.generate_candidates_table(
                diff=image,
                scorr_catalog_path=cands_catalog_name,
                sci_resamp_image_path=sci_image_path,
                ref_resamp_image_path=ref_image_path,
                diff_scorr_path=scorr_image_path,
                diff_psf_path=diff_psf_path,
                diff_unc_path=diff_unc_path,
            )

            if len(srcs_table) > 0:
                x_shape, y_shape = image.get_data().shape
                srcs_table["X_SHAPE"] = x_shape
                srcs_table["Y_SHAPE"] = y_shape

            metadata = {}

            for key in image.keys():
                if key != "COMMENT":
                    metadata[key] = image[key]

            if len(srcs_table) == 0:
                msg = f"No sources found in image {image[BASE_NAME_KEY]}"
                logger.warning(msg)

            else:
                msg = f"Found {len(srcs_table)} sources in image {image[BASE_NAME_KEY]}"
                logger.debug(msg)
                all_cands.append(SourceTable(srcs_table, metadata=metadata))

        return all_cands
