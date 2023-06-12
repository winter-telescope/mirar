"""
Module with processors to perform aperture photometry
"""
import shutil
from typing import Optional

import numpy as np

from mirar.data import ImageBatch, SourceBatch
from mirar.paths import (
    APFLUX_PREFIX_KEY,
    APFLUXUNC_PREFIX_KEY,
    APMAG_PREFIX_KEY,
    APMAGUNC_PREFIX_KEY,
    ZP_KEY,
)
from mirar.processors.photometry.base_photometry import (
    AperturePhotometry,
    BaseCandidatePhotometry,
    BaseImagePhotometry,
)


class CandidateAperturePhotometry(BaseCandidatePhotometry):
    """
    Processor to run aperture photometry on all candidates in candidate table
    """

    base_key = "APERPHOTDF"

    def __init__(
        self,
        *args,
        aper_diameters: float | list[float] = 10.0,
        bkg_in_diameters: float | list[float] = 25.0,
        bkg_out_diameters: float | list[float] = 40.0,
        zp_colname: str = ZP_KEY,
        col_suffix_list: str | list[str] = None,
        **kwargs,
    ):
        super().__init__(*args, **kwargs)

        self.aperture_photometer = AperturePhotometry(
            aper_diameters=aper_diameters,
            bkg_in_diameters=bkg_in_diameters,
            bkg_out_diameters=bkg_out_diameters,
        )
        self.col_suffix_list = col_suffix_list
        if self.col_suffix_list is None:
            self.col_suffix_list = self.aperture_photometer.aper_diameters
        self.zp_colname = zp_colname

    def _apply_to_candidates(
        self,
        batch: SourceBatch,
    ) -> SourceBatch:
        for source_table in batch:
            candidate_table = source_table.get_data()
            all_fluxes, all_fluxuncs = [], []
            for cand_ind in range(len(candidate_table)):
                row = candidate_table.iloc[cand_ind]

                image_cutout, unc_image_cutout = self.generate_cutouts(row)

                fluxes, fluxuncs = self.aperture_photometer.perform_photometry(
                    image_cutout=image_cutout, unc_image_cutout=unc_image_cutout
                )
                all_fluxes.append(fluxes)
                all_fluxuncs.append(fluxuncs)
            all_fluxes = np.array(all_fluxes).T
            all_fluxuncs = np.array(all_fluxuncs).T

            for ind, suffix in enumerate(self.col_suffix_list):
                flux, fluxunc = all_fluxes[ind], all_fluxuncs[ind]
                candidate_table[f"{APFLUX_PREFIX_KEY}{suffix}"] = flux
                candidate_table[f"{APFLUXUNC_PREFIX_KEY}{suffix}"] = fluxunc
                candidate_table[f"{APMAG_PREFIX_KEY}{suffix}"] = np.array(
                    candidate_table[self.zp_colname], dtype=float
                ) - 2.5 * np.log10(flux)
                candidate_table[f"{APMAGUNC_PREFIX_KEY}{suffix}"] = (
                    1.086 * fluxunc / flux
                )

            source_table.set_data(candidate_table)

        if self.photometry_out_temp_dir is not None:
            shutil.rmtree(self.photometry_out_temp_dir)
        return batch


class ImageAperturePhotometry(BaseImagePhotometry):
    """
    Processor to run aperture photometry at the RA/Dec specified in the header
    """

    base_key = "APERPHOTIM"

    def __init__(
        self,
        *args,
        aper_diameters: float | list[float] = 10.0,
        bkg_in_diameters: float | list[float] = 25.0,
        bkg_out_diameters: float | list[float] = 40.0,
        zp_colname: str = ZP_KEY,
        col_suffix_list: Optional[list[str]] = None,
        **kwargs,
    ):
        super().__init__(*args, **kwargs)

        self.aperture_photometer = AperturePhotometry(
            aper_diameters=aper_diameters,
            bkg_in_diameters=bkg_in_diameters,
            bkg_out_diameters=bkg_out_diameters,
        )
        self.col_suffix_list = col_suffix_list
        if self.col_suffix_list is None:
            self.col_suffix_list = self.aperture_photometer.aper_diameters

        self.zp_colname = zp_colname

    def _apply_to_images(
        self,
        batch: ImageBatch,
    ) -> ImageBatch:
        for image in batch:
            image_cutout, unc_image_cutout = self.generate_cutouts(image)

            fluxes, fluxuncs = self.aperture_photometer.perform_photometry(
                image_cutout, unc_image_cutout
            )

            for ind, flux in enumerate(fluxes):
                fluxunc = fluxuncs[ind]
                suffix = self.col_suffix_list[ind]
                image[f"fluxap{suffix}"] = flux
                image[f"fluxunc{suffix}"] = fluxunc
                image[f"magap{suffix}"] = float(
                    image[self.zp_colname]
                ) - 2.5 * np.log10(flux)
                image[f"magerrap{suffix}"] = 1.086 * fluxunc / flux

        if self.photometry_out_temp_dir is not None:
            shutil.rmtree(self.photometry_out_temp_dir)
        return batch
