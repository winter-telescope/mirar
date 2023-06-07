"""
Module with processors to perform point-spread-function photometry
"""
import logging
import shutil

import numpy as np

from mirar.data import Image, ImageBatch, SourceBatch
from mirar.paths import (
    MAG_PSF_KEY,
    MAGERR_PSF_KEY,
    NORM_PSFEX_KEY,
    PSF_FLUX_KEY,
    PSF_FLUXUNC_KEY,
    ZP_KEY,
)
from mirar.processors.astromatic.psfex import PSFex
from mirar.processors.base_processor import PrerequisiteError
from mirar.processors.photometry.base_photometry import (
    BaseCandidatePhotometry,
    BaseImagePhotometry,
    PSFPhotometry,
)

logger = logging.getLogger(__name__)


def check_psf_phot_prerequisites(processor):
    """
    Function to check prerequisites for running PSF photometry
    Args:
        processor: PSF photometry processor

    """
    mask = [isinstance(x, PSFex) for x in processor.preceding_steps]
    if np.sum(mask) < 1:
        err = (
            f"{processor.__module__} requires {PSFex} as a prerequisite. "
            f"However, the following steps were found: {processor.preceding_steps}."
        )
        logger.error(err)
        raise PrerequisiteError(err)


class CandidatePSFPhotometry(BaseCandidatePhotometry):
    """
    Processor to run PSF photometry on all candidates in candidate table
    """

    base_key = "PSFPHOTDF"

    def __init__(self, zp_colname=ZP_KEY):
        super().__init__()
        self.zp_colname = zp_colname

    def get_psf_filename(self, row):
        """
        Function to get the name of psf file
        Args:
            row: row of a pandas Dataframe

        Returns:

        """
        psf_filename = row[self.psf_file_key]
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

            candidate_table[PSF_FLUX_KEY] = fluxes
            candidate_table[PSF_FLUXUNC_KEY] = fluxuncs
            candidate_table["chipsf"] = minchi2s
            candidate_table["xshift"] = xshifts
            candidate_table["yshift"] = yshifts
            candidate_table[MAG_PSF_KEY] = np.array(
                candidate_table[self.zp_colname], dtype=float
            ) - 2.5 * np.log10(candidate_table[PSF_FLUX_KEY])
            candidate_table[MAGERR_PSF_KEY] = (
                1.086 * candidate_table[PSF_FLUXUNC_KEY] / candidate_table[PSF_FLUX_KEY]
            )

            source_table.set_data(candidate_table)
        if self.photometry_out_temp_dir is not None:
            shutil.rmtree(self.photometry_out_temp_dir)
        return batch

    def check_prerequisites(
        self,
    ):
        check_psf_phot_prerequisites(self)


class ImagePSFPhotometry(BaseImagePhotometry):
    """
    Processor to run PSF photometry at the RA/Dec specified in the header
    """

    base_key = "PSFPHOTIM"

    def __init__(
        self,
        *args,
        zp_colname: str = ZP_KEY,
        **kwargs,
    ):
        super().__init__(*args, **kwargs)

        self.zp_colname = zp_colname

    def get_psf_filename(self, image: Image):
        """
        Function to get PSF file name of an image
        Args:
            image: Image

        Returns:

        """
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
            image[MAG_PSF_KEY] = -2.5 * np.log10(flux) + float(image[self.zp_colname])
            image[MAGERR_PSF_KEY] = 1.086 * fluxunc / flux

        if self.photometry_out_temp_dir is not None:
            shutil.rmtree(self.photometry_out_temp_dir)
        return batch

    def check_prerequisites(
        self,
    ):
        check_psf_phot_prerequisites(self)
