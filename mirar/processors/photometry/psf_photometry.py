"""
Module with processors to perform point-spread-function photometry
"""
import logging
from pathlib import Path

import numpy as np

from mirar.data import SourceBatch
from mirar.paths import MAG_PSF_KEY, MAGERR_PSF_KEY, PSF_FLUX_KEY, PSF_FLUXUNC_KEY
from mirar.processors.astromatic.psfex import PSFex
from mirar.processors.base_processor import PrerequisiteError
from mirar.processors.photometry.base_photometry import BasePhotometryProcessor
from mirar.processors.photometry.utils import (
    get_mags_from_fluxes,
    make_psf_shifted_array,
    psf_photometry,
)

logger = logging.getLogger(__name__)


class SourcePSFPhotometry(BasePhotometryProcessor):
    """
    Processor to run PSF photometry on a source table
    """

    base_key = "PSFPHOT"

    def perform_photometry(
        self, image_cutout: np.array, unc_image_cutout: np.array, psf_filename: Path
    ) -> tuple[float, float, float, float, float]:
        psfmodels = make_psf_shifted_array(
            psf_filename=psf_filename,
            cutout_size_psf_phot=int(image_cutout.shape[0] / 2),
        )

        flux, fluxunc, minchi2, xshift, yshift = psf_photometry(
            image_cutout, unc_image_cutout, psfmodels
        )
        return flux, fluxunc, minchi2, xshift, yshift

    def get_psf_filename(self, row):
        """
        Function to get the name of psf file
        Args:
            row: row of a pandas Dataframe

        Returns:

        """
        psf_filename = row[self.psf_file_key]
        return psf_filename

    def _apply_to_sources(
        self,
        batch: SourceBatch,
    ) -> SourceBatch:
        for source_table in batch:
            candidate_table = source_table.get_data()

            metadata = source_table.get_metadata()

            fluxes, fluxuncs, minchi2s, xshifts, yshifts = [], [], [], [], []

            for _, row in candidate_table.iterrows():
                image_cutout, unc_image_cutout = self.generate_cutouts(row, metadata)
                psf_filename = self.get_psf_filename(row)
                (
                    flux,
                    fluxunc,
                    minchi2,
                    xshift,
                    yshift,
                ) = self.perform_photometry(
                    image_cutout, unc_image_cutout, psf_filename=psf_filename
                )
                fluxes.append(flux)
                fluxuncs.append(fluxunc)
                minchi2s.append(minchi2)
                xshifts.append(xshift)
                yshifts.append(yshift)

            candidate_table[PSF_FLUX_KEY] = fluxes
            candidate_table[PSF_FLUXUNC_KEY] = fluxuncs
            candidate_table["chipsf"] = minchi2s
            candidate_table["xshift"] = xshifts
            candidate_table["yshift"] = yshifts

            magnitudes, magnitudes_unc = get_mags_from_fluxes(
                flux_list=fluxes,
                fluxunc_list=fluxuncs,
                zeropoint_list=np.array(candidate_table[self.zp_key], dtype=float),
                zeropoint_unc_list=np.array(
                    candidate_table[self.zp_std_key], dtype=float
                ),
            )

            candidate_table[MAG_PSF_KEY] = magnitudes
            candidate_table[MAGERR_PSF_KEY] = magnitudes_unc

            source_table.set_data(candidate_table)

        return batch

    def check_prerequisites(
        self,
    ):
        mask = [isinstance(x, PSFex) for x in self.preceding_steps]
        if np.sum(mask) < 1:
            err = (
                f"{self.__module__} requires {PSFex} as a prerequisite. "
                f"However, the following steps were found: {self.preceding_steps}."
            )
            logger.error(err)
            raise PrerequisiteError(err)
