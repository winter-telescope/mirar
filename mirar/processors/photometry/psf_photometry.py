"""
Module with processors to perform point-spread-function photometry
"""

import logging
from pathlib import Path

import numpy as np

from mirar.data import SourceBatch
from mirar.paths import (
    MAG_PSF_KEY,
    MAGERR_PSF_KEY,
    PSF_FLUX_KEY,
    PSF_FLUXUNC_KEY,
    get_output_dir,
)
from mirar.processors.base_processor import PrerequisiteError
from mirar.processors.photometry.base_photometry import BasePhotometryProcessor
from mirar.processors.photometry.utils import (
    get_mags_from_fluxes,
    make_psf_shifted_array,
    psf_photometry,
)

logger = logging.getLogger(__name__)


class PSFPhotometry(BasePhotometryProcessor):
    """
    Processor to run PSF photometry on a source table
    """

    base_key = "PSFPHOT"

    def perform_photometry(
        self,
        image_cutout: np.array,
        unc_image_cutout: np.array,
        psf_filename: str | Path,
    ) -> tuple[float, float, float, float, float]:
        """
        Function to perform PSF photometry on a cutout
        :param image_cutout: cutout of image
        :param unc_image_cutout: cutout of uncertainty image
        :param psf_filename: filename of psf file
        :return: flux, fluxunc, minchi2, xshift, yshift
        """
        if not isinstance(psf_filename, Path):
            psf_filename = Path(psf_filename)
        psfmodels = make_psf_shifted_array(
            psf_filename=psf_filename.as_posix(),
            cutout_size_psf_phot=int(image_cutout.shape[0] / 2),
        )

        flux, fluxunc, minchi2, xshift, yshift = psf_photometry(
            image_cutout,
            unc_image_cutout,
            psfmodels,
            psf_filename=psf_filename.as_posix(),
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

            if self.psf_file_key not in metadata:
                raise PrerequisiteError(
                    f"PSF file key {self.psf_file_key} not in source table."
                    f"Have you run the PSFEx processor, or set the correct key with"
                    f" the psf file name?"
                )
            psf_filename = source_table[self.psf_file_key]
            temp_imagename, temp_unc_imagename = self.save_temp_image_uncimage(metadata)

            for ind, row in candidate_table.iterrows():
                image_cutout, unc_image_cutout = self.generate_cutouts(
                    imagename=temp_imagename,
                    unc_imagename=temp_unc_imagename,
                    data_item=row,
                )
                (
                    flux,
                    fluxunc,
                    minchi2,
                    xshift,
                    yshift,
                ) = self.perform_photometry(
                    image_cutout, unc_image_cutout, psf_filename=psf_filename
                )

                if self.save_cutouts:
                    image_cutout_path = get_output_dir(
                        self.temp_output_sub_dir, self.night_sub_dir
                    ).joinpath(f"image_cutout_{ind}.dat")
                    logger.debug(f"Writing cutout to {image_cutout_path}")
                    np.savetxt(X=image_cutout, fname=image_cutout_path)
                    unc_image_cutout_path = get_output_dir(
                        self.temp_output_sub_dir, self.night_sub_dir
                    ).joinpath(f"unc_image_cutout_{ind}.dat")
                    logger.debug(f"Writing cutout to {unc_image_cutout_path}")
                    np.savetxt(X=unc_image_cutout, fname=unc_image_cutout_path)

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
                flux_list=np.array(fluxes, dtype=float),
                fluxunc_list=np.array(fluxuncs, dtype=float),
                zeropoint=float(source_table[self.zp_key]),
                zeropoint_unc=float(source_table[self.zp_std_key]),
            )

            candidate_table[MAG_PSF_KEY] = magnitudes
            candidate_table[MAGERR_PSF_KEY] = magnitudes_unc

            temp_imagename.unlink()
            temp_unc_imagename.unlink()

            source_table.set_data(candidate_table)

        return batch
