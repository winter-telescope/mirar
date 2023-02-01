from typing import Optional

import numpy as np
from winterdrp.data import ImageBatch, SourceBatch
from winterdrp.paths import ZP_KEY
from winterdrp.processors.photometry.base_photometry import (
    AperturePhotometry,
    BaseCandidatePhotometry,
    BaseImagePhotometry,
)
from winterdrp.processors.photometry.utils import make_cutouts


class CandidateAperturePhotometry(BaseCandidatePhotometry):
    base_key = "APERPHOTDF"
    def __init__(
        self,
        aper_diameters: float | list[float] = 10.0,
        bkg_in_diameters: float | list[float] = 25.0,
        bkg_out_diameters: float | list[float] = 40.0,
        col_suffix_list: Optional[list[str]] = None,
        zp_colname="magzpsci",
        *args,
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
                imagename, unc_imagename = self.get_filenames(row)
                x, y = self.get_physical_coordinates(row)

                image_cutout, unc_image_cutout = make_cutouts(
                    image_paths=[imagename, unc_imagename],
                    position=(x, y),
                    half_size=self.phot_cutout_size,
                )

                fluxes, fluxuncs = self.aperture_photometer.perform_photometry(
                    image_cutout=image_cutout, unc_image_cutout=unc_image_cutout
                )
                all_fluxes.append(fluxes)
                all_fluxuncs.append(fluxuncs)
            all_fluxes = np.array(all_fluxes).T
            all_fluxuncs = np.array(all_fluxuncs).T

            for ind, suffix in enumerate(self.col_suffix_list):
                flux, fluxunc = all_fluxes[ind], all_fluxuncs[ind]
                candidate_table[f"fluxap{suffix}"] = flux
                candidate_table[f"fluxuncap{suffix}"] = fluxunc
                candidate_table[f"magap{suffix}"] = candidate_table[
                    self.zp_colname
                ] - 2.5 * np.log10(flux)
                candidate_table[f"sigmagap{suffix}"] = 1.086 * fluxunc / flux

            source_table.set_data(candidate_table)
        return batch


class ImageAperturePhotometry(BaseImagePhotometry):
    def __init__(
        self,
        aper_diameters: float | list[float] = 10.0,
        bkg_in_diameters: float | list[float] = 25.0,
        bkg_out_diameters: float | list[float] = 40.0,
        col_suffix_list: Optional[list[str]] = None,
    ):
        super().__init__()

        self.aperture_photometer = AperturePhotometry(
            aper_diameters=aper_diameters,
            bkg_in_diameters=bkg_in_diameters,
            bkg_out_diameters=bkg_out_diameters,
        )
        self.col_suffix_list = col_suffix_list
        if self.col_suffix_list is None:
            self.col_suffix_list = self.aperture_photometer.aper_diameters

    def _apply_to_images(
        self,
        batch: ImageBatch,
    ) -> ImageBatch:
        for image in batch:
            imagename, unc_imagename = self.get_filenames(image)
            x, y = self.get_physical_coordinates(image)
            image_cutout, unc_image_cutout = make_cutouts(
                image_paths=[imagename, unc_imagename],
                position=(x, y),
                half_size=self.phot_cutout_size,
            )

            fluxes, fluxuncs = self.aperture_photometer.perform_photometry(
                image_cutout, unc_image_cutout
            )

            for ind in range(len(fluxes)):
                flux, fluxunc = fluxes[ind], fluxuncs[ind]
                suffix = self.col_suffix_list[ind]
                image[f"fluxap{suffix}"] = flux
                image[f"fluxunc{suffix}"] = fluxunc
                image[f"magap{suffix}"] = image[ZP_KEY] - 2.5 * np.log10(flux)
                image[f"magerrap{suffix}"] = 1.086 * fluxunc / flux

        return batch
