"""
Module to detect candidates in an image
"""
import gzip
import io
import logging
import os
from pathlib import Path

import numpy as np
import pandas as pd
from astropy.io import fits

from winterdrp.data import ImageBatch, SourceBatch, SourceTable
from winterdrp.paths import BASE_NAME_KEY, REF_IMG_KEY, core_fields, get_output_dir
from winterdrp.processors.astromatic.sextractor.sourceextractor import (
    run_sextractor_dual,
)
from winterdrp.processors.base_processor import BaseCandidateGenerator
from winterdrp.utils.ldac_tools import get_table_from_ldac

# TODO : Move photometry to its own thing like catalogs, user can choose
# whichever way they want to do photometry
logger = logging.getLogger(__name__)


class DetectCandidates(BaseCandidateGenerator):
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
        return f"<Extracts detected sources from images, and converts them to a pandas dataframe>"

    def get_sub_output_dir(self):
        return Path(get_output_dir(self.output_sub_dir, self.night_sub_dir))

    def get_path(self, name: str) -> Path:
        return self.get_sub_output_dir().joinpath(name)

    def make_alert_cutouts(self, image_path: str, position, half_size):
        data = self.open_fits(image_path).get_data()

        y_image_size, x_image_size = np.shape(data)
        x, y = position
        # logger.debug(f'{x},{y},{np.shape(data)}')
        if y < half_size:
            cutout = data[0 : y + half_size + 1, x - half_size : x + half_size + 1]
            n_pix = half_size - y
            cutout = np.pad(cutout, ((n_pix, 0), (0, 0)), "constant")

        elif y + half_size + 1 > y_image_size:
            cutout = data[
                y - half_size : y_image_size, x - half_size : x + half_size + 1
            ]
            n_pix = (half_size + y + 1) - y_image_size
            cutout = np.pad(cutout, ((0, n_pix), (0, 0)), "constant")

        elif x < half_size:
            cutout = data[y - half_size : y + half_size + 1, 0 : x + half_size + 1]
            n_pix = half_size - x
            cutout = np.pad(cutout, ((0, 0), (n_pix, 0)), "constant")
        elif x + half_size > x_image_size:
            cutout = data[
                y - half_size : y + half_size + 1, x - half_size : x_image_size
            ]
            n_pix = (half_size + x + 1) - x_image_size
            cutout = np.pad(cutout, ((0, 0), (0, n_pix)), "constant")
        else:
            cutout = data[
                y - half_size : y + half_size + 1, x - half_size : x + half_size + 1
            ]
        return cutout

    @staticmethod
    def makebitims(image):
        ######################################################
        # make bit images of the cutouts for the marshal
        #
        # Inputs:
        # image: input image cutout
        #
        # Returns:
        # buf2: a gzipped fits file of the cutout image as
        #  a BytesIO object
        ######################################################

        # open buffer and store image in memory
        buf = io.BytesIO()
        buf2 = io.BytesIO()
        fits.writeto(buf, image)
        with gzip.open(buf2, "wb") as fz:
            fz.write(buf.getvalue())

        return buf2

    def generate_candidates_table(
        self,
        scorr_catalog_name: str | Path,
        sci_resamp_path: str | Path,
        ref_resamp_path: str | Path,
        diff_path: str | Path,
        diff_scorr_path: str | Path,
        diff_psf_path: str | Path,
        diff_unc_path: str | Path,
    ) -> pd.DataFrame:
        det_src_table = get_table_from_ldac(scorr_catalog_name)

        if len(det_src_table) == 0:
            return pd.DataFrame()

        logger.info(f"Found {len(det_src_table)} candidates in image {diff_path}.")
        det_src_table["xpos"] = det_src_table["X_IMAGE"] - 1
        det_src_table["ypos"] = det_src_table["Y_IMAGE"] - 1

        diff_image = self.open_fits(diff_path)
        sci_resamp_image = self.open_fits(sci_resamp_path)

        scorr_data = self.open_fits(diff_scorr_path).get_data()

        xpeaks, ypeaks = (
            det_src_table["XPEAK_IMAGE"] - 1,
            det_src_table["YPEAK_IMAGE"] - 1,
        )
        scorr_peaks = scorr_data[ypeaks, xpeaks]
        det_src_table["xpeak"] = xpeaks
        det_src_table["ypeak"] = ypeaks
        det_src_table["scorr"] = scorr_peaks

        cutout_size_display = 40

        display_sci_ims = []
        display_ref_ims = []
        display_diff_ims = []

        for ind, src in enumerate(det_src_table):
            xpeak, ypeak = int(xpeaks[ind]), int(ypeaks[ind])

            display_sci_cutout = self.make_alert_cutouts(
                sci_resamp_path, (xpeak, ypeak), cutout_size_display
            )
            display_ref_cutout = self.make_alert_cutouts(
                ref_resamp_path, (xpeak, ypeak), cutout_size_display
            )
            display_diff_cutout = self.make_alert_cutouts(
                diff_path, (xpeak, ypeak), cutout_size_display
            )

            display_sci_bit = self.makebitims(display_sci_cutout.astype(np.float32))
            display_ref_bit = self.makebitims(display_ref_cutout.astype(np.float32))
            display_diff_bit = self.makebitims(display_diff_cutout.astype(np.float32))

            display_sci_ims.append(display_sci_bit)
            display_ref_ims.append(display_ref_bit)
            display_diff_ims.append(display_diff_bit)

        det_src_table["cutoutScience"] = display_sci_ims
        det_src_table["cutoutTemplate"] = display_ref_ims
        det_src_table["cutoutDifference"] = display_diff_ims

        diff_zp = float(diff_image["ZP"])
        det_src_table["magzpsci"] = diff_zp
        diff_zp_unc = float(diff_image["ZP_std"])
        det_src_table["magzpsciunc"] = diff_zp_unc
        det_src_table["diffimname"] = diff_path
        det_src_table["sciimname"] = sci_resamp_path
        det_src_table["refimname"] = ref_resamp_path
        det_src_table["diffpsfname"] = diff_psf_path
        det_src_table["diffuncname"] = diff_unc_path

        det_src_table["ra"] = det_src_table["ALPHA_J2000"]
        det_src_table["dec"] = det_src_table["DELTA_J2000"]
        det_src_table["fwhm"] = det_src_table["FWHM_IMAGE"]
        det_src_table["aimage"] = det_src_table["A_IMAGE"]
        det_src_table["bimage"] = det_src_table["B_IMAGE"]
        det_src_table["aimagerat"] = det_src_table["aimage"] / det_src_table["fwhm"]
        det_src_table["bimagerat"] = det_src_table["bimage"] / det_src_table["fwhm"]
        det_src_table["elong"] = det_src_table["ELONGATION"]

        det_src_table["jd"] = sci_resamp_image["MJD-OBS"] + 2400000.5
        det_src_table["exptime"] = diff_image["EXPTIME"]
        det_src_table["field"] = sci_resamp_image["FIELDID"]
        det_src_table["programpi"] = sci_resamp_image["PROGPI"]
        det_src_table["programid"] = sci_resamp_image["PROGID"]
        det_src_table["fid"] = sci_resamp_image["FILTERID"]
        det_src_table["candid"] = np.array(
            det_src_table["jd"] * 100, dtype=int
        ) * 10000 + np.arange(len(det_src_table))
        det_src_table["diffmaglim"] = diff_image["DIFFMLIM"]
        det_src_table["isdiffpos"] = 1
        det_src_table = det_src_table.to_pandas()
        logger.info(det_src_table["diffmaglim"])

        return det_src_table

    def _apply_to_images(
        self,
        batch: ImageBatch,
    ) -> SourceBatch:
        all_cands = SourceBatch()
        for image in batch:

            scorr_image_path = self.get_path(image["DIFFSCR"])
            diff_image_path = self.get_path(image["DIFFIMG"])
            diff_psf_path = self.get_path(image["DIFFPSF"])
            diff_unc_path = self.get_path(image["DIFFUNC"])

            scorr_mask_path = self.get_path(image["SCORMASK"])
            cands_catalog_name = str(diff_image_path).replace(".fits", ".dets")

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

            sci_image_path = self.get_path(image[BASE_NAME_KEY])
            ref_image_path = self.get_path(image[REF_IMG_KEY])

            cands_table = self.generate_candidates_table(
                scorr_catalog_name=cands_catalog_name,
                sci_resamp_path=sci_image_path,
                ref_resamp_path=ref_image_path,
                diff_path=diff_image_path,
                diff_scorr_path=scorr_image_path,
                diff_psf_path=diff_psf_path,
                diff_unc_path=diff_unc_path,
            )

            if len(cands_table) > 0:
                x_shape, y_shape = image.get_data().shape

                cands_table["X_SHAPE"] = x_shape
                cands_table["Y_SHAPE"] = y_shape

            metadata = {}

            for key in core_fields:
                metadata[key] = image[key]

            all_cands.append(SourceTable(cands_table, metadata=metadata))

        return all_cands
