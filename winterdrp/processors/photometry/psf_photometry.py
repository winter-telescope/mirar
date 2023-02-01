"""
Module with processors to perform point-spread-function photometry
"""
import numpy as np

from winterdrp.data import Image, ImageBatch, SourceBatch
from winterdrp.paths import (
    MAG_PSF_KEY,
    MAGERR_PSF_KEY,
    NORM_PSFEX_KEY,
    PSF_FLUX_KEY,
    PSF_FLUXUNC_KEY,
    ZP_KEY,
)
from winterdrp.processors.photometry.base_photometry import (
    BaseCandidatePhotometry,
    BaseImagePhotometry,
    PSFPhotometry,
)


class CandidatePSFPhotometry(BaseCandidatePhotometry):
    """
    Processor to run PSF photometry on all candidates in candidate table
    """
    def __init__(self, zp_colname="magzpsci"):
        super().__init__()
        self.zp_colname = zp_colname

    def get_psf_filename(self, row):
        psf_filename = row[self.psf_file_colname]
        return psf_filename

    def _apply_to_candidates(
        self,
        batch: SourceBatch,
    ) -> SourceBatch:

        for source_table in batch:
            candidate_table = source_table.get_data()

            fluxes, fluxuncs, minchi2s, xshifts, yshifts = [], [], [], [], []

            for ind in range(len(candidate_table)):
                row = candidate_table.iloc[ind]

                image_cutout, unc_image_cutout = self.generate_cutouts(row)
                psf_filename = self.get_psf_filename(row)
                psf_photometer = PSFPhotometry(psf_filename=psf_filename)
                (
                    flux,
                    fluxunc,
                    minchi2,
                    xshift,
                    yshift,
                ) = psf_photometer.perform_photometry(image_cutout, unc_image_cutout)
                fluxes.append(flux)
                fluxuncs.append(fluxunc)
                minchi2s.append(minchi2)
                xshifts.append(xshift)
                yshifts.append(yshift)

            candidate_table["psf_flux"] = fluxes
            candidate_table["psf_fluxunc"] = fluxuncs
            candidate_table["chipsf"] = minchi2s
            candidate_table["xshift"] = xshifts
            candidate_table["yshift"] = yshifts
            candidate_table["magpsf"] = candidate_table[
                self.zp_colname
            ] - 2.5 * np.log10(candidate_table["psf_flux"])
            candidate_table["sigmapsf"] = (
                1.086 * candidate_table["psf_fluxunc"] / candidate_table["psf_flux"]
            )

            source_table.set_data(candidate_table)
        return batch


class ImagePSFPhotometry(BaseImagePhotometry):
    """
    Processor to run PSF photometry at the RA/Dec specified in the header
    """

    def get_psf_filename(self, image: Image):
        psf_filename = image[NORM_PSFEX_KEY]
        return psf_filename

    def _apply_to_images(
        self,
        batch: ImageBatch,
    ) -> ImageBatch:
        for image in batch:
            image_cutout, unc_image_cutout = self.generate_cutouts(image)
            psf_filename = self.get_psf_filename(image)

            psf_photometer = PSFPhotometry(psf_filename=psf_filename)

            flux, fluxunc, _, _, _ = psf_photometer.perform_photometry(
                image_cutout, unc_image_cutout
            )

            image[PSF_FLUX_KEY] = flux
            image[PSF_FLUXUNC_KEY] = fluxunc
            image[MAG_PSF_KEY] = -2.5 * np.log10(flux) + image[ZP_KEY]
            image[MAGERR_PSF_KEY] = 1.086 * fluxunc / flux

        return batch
