"""
Module with processors to perform aperture photometry
"""

import numpy as np

from mirar.data import SourceBatch
from mirar.paths import (
    APFLUX_PREFIX_KEY,
    APFLUXUNC_PREFIX_KEY,
    APMAG_PREFIX_KEY,
    APMAGUNC_PREFIX_KEY,
    get_output_dir,
)
from mirar.processors.photometry.base_photometry import BasePhotometryProcessor
from mirar.processors.photometry.utils import aper_photometry, get_mags_from_fluxes


class AperturePhotometry(BasePhotometryProcessor):
    """
    Processor to run aperture photometry on all candidates in candidate table
    """

    base_key = "APERPHOT"

    def __init__(
        self,
        *args,
        aper_diameters: float | list[float] = 10.0,
        bkg_in_diameters: float | list[float] = 25.0,
        bkg_out_diameters: float | list[float] = 40.0,
        col_suffix_list: str | list[str] = None,
        **kwargs,
    ):
        super().__init__(*args, **kwargs)

        if not isinstance(aper_diameters, list):
            aper_diameters = [aper_diameters]
        if not isinstance(bkg_in_diameters, list):
            bkg_in_diameters = [bkg_in_diameters]
        if not isinstance(bkg_out_diameters, list):
            bkg_out_diameters = [bkg_out_diameters]

        self.aper_diameters = aper_diameters
        self.bkg_in_diameters = bkg_in_diameters
        self.bkg_out_diameters = bkg_out_diameters

        self.col_suffix_list = col_suffix_list
        if self.col_suffix_list is None:
            self.col_suffix_list = self.aper_diameters

    def perform_photometry(
        self, image_cutout: np.array, unc_image_cutout: np.array
    ) -> tuple[list[float], list[float]]:
        fluxes, fluxuncs = [], []
        for ind, aper_diam in enumerate(self.aper_diameters):
            flux, fluxunc = aper_photometry(
                image_cutout,
                unc_image_cutout,
                aper_diam,
                self.bkg_in_diameters[ind],
                self.bkg_out_diameters[ind],
            )
            fluxes.append(flux)
            fluxuncs.append(fluxunc)
        return fluxes, fluxuncs

    def _apply_to_sources(
        self,
        batch: SourceBatch,
    ) -> SourceBatch:
        for source_table in batch:
            candidate_table = source_table.get_data()

            metadata = source_table.get_metadata()

            all_fluxes, all_fluxuncs = [], []
            temp_imagename, temp_unc_imagename = self.save_temp_image_uncimage(metadata)

            for cand_ind in range(len(candidate_table)):
                row = candidate_table.iloc[cand_ind]

                image_cutout, unc_image_cutout = self.generate_cutouts(
                    imagename=temp_imagename,
                    unc_imagename=temp_unc_imagename,
                    data_item=row,
                )

                fluxes, fluxuncs = self.perform_photometry(
                    image_cutout=image_cutout, unc_image_cutout=unc_image_cutout
                )
                all_fluxes.append(fluxes)
                all_fluxuncs.append(fluxuncs)

                if self.save_cutouts:
                    image_cutout_path = get_output_dir(
                        self.temp_output_sub_dir, self.night_sub_dir
                    ).joinpath(f"image_cutout_{cand_ind}.dat")
                    np.savetxt(X=image_cutout, fname=image_cutout_path)
                    unc_image_cutout_path = get_output_dir(
                        self.temp_output_sub_dir, self.night_sub_dir
                    ).joinpath(f"unc_image_cutout_{cand_ind}.dat")
                    np.savetxt(X=unc_image_cutout, fname=unc_image_cutout_path)

            all_fluxes = np.array(all_fluxes).T
            all_fluxuncs = np.array(all_fluxuncs).T

            for ind, suffix in enumerate(self.col_suffix_list):
                flux, fluxunc = all_fluxes[ind], all_fluxuncs[ind]
                candidate_table[f"{APFLUX_PREFIX_KEY}{suffix}"] = flux
                candidate_table[f"{APFLUXUNC_PREFIX_KEY}{suffix}"] = fluxunc

                magnitudes, magnitudes_unc = get_mags_from_fluxes(
                    flux_list=np.array(flux, dtype=float),
                    fluxunc_list=np.array(fluxunc, dtype=float),
                    zeropoint=float(source_table[self.zp_key]),
                    zeropoint_unc=float(source_table[self.zp_std_key]),
                )
                candidate_table[f"{APMAG_PREFIX_KEY}{suffix}"] = magnitudes
                candidate_table[f"{APMAGUNC_PREFIX_KEY}{suffix}"] = magnitudes_unc

            temp_imagename.unlink()
            temp_unc_imagename.unlink()
            source_table.set_data(candidate_table)

        return batch
