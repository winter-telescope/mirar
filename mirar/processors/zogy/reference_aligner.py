from pathlib import Path
from typing import Callable

import astropy
from scipy.ndimage import shift

from mirar.data import Image, ImageBatch
from mirar.paths import BASE_NAME_KEY, REF_IMG_KEY, SEXTRACTOR_HEADER_KEY
from mirar.processors.astromatic import PSFex, Sextractor
from mirar.processors.reference import logger
from mirar.processors.zogy.zogy import ZOGYPrepare, default_catalog_purifier


class AlignReference(ZOGYPrepare):
    """
    Processor to calculate an offset betwen the science image and its reference image
    (path taken from header) and shift the reference image using scipy.ndimage.shift.
    The shifted reference image is then saved to a directory and path is added to a
    header. A new PSF model and source catalog is generated for the shifted reference.
    """

    base_key = "ref_aligner"

    def __init__(
        self,
        sextractor: Callable[..., Sextractor],
        psfex: Callable[..., PSFex],
        phot_sextractor: Callable[..., Sextractor] = None,
        temp_output_subtract_dir: str = "subtract",
        catalog_purifier: Callable[
            [astropy.table.Table, astropy.table.Table],
            list[astropy.table.Table, astropy.table.Table],
        ] = default_catalog_purifier,
        order: int = 1,
    ):
        super().__init__(catalog_purifier=catalog_purifier)
        self.sextractor = sextractor
        self.psfex = psfex
        self.phot_sextractor = phot_sextractor
        self.temp_output_subtract_dir = temp_output_subtract_dir
        self.order = order

    def _apply_to_images(
        self,
        batch: ImageBatch,
    ) -> ImageBatch:
        for image in batch:
            ref_path = Path(image[REF_IMG_KEY])
            logger.debug(f"Aligning {image[BASE_NAME_KEY]} with {ref_path}")
            ref_image = self.open_fits(ref_path)
            ref_catalog_path = ref_image[SEXTRACTOR_HEADER_KEY]

            sci_catalog_path = image[SEXTRACTOR_HEADER_KEY]

            _, _, median_x_offset, median_y_offset, _ = self.get_ast_fluxscale(
                ref_catalog_name=ref_catalog_path, sci_catalog_name=sci_catalog_path
            )

            logger.debug(
                f"Calculated median x offset: {median_x_offset}"
                f" and median y offset: {median_y_offset} between "
                f"{ref_path} and {image[BASE_NAME_KEY]}"
            )
            # Shift the reference image
            ref_shifted_data = shift(
                ref_image.get_data(),
                (median_y_offset, median_x_offset),
                order=self.order,
                mode="reflect",
                cval=0.0,
                prefilter=True,
            )

            ref_shifted_image = Image(
                data=ref_shifted_data, header=ref_image.get_header()
            )
            ref_shifted_image["CRPIX1"] += median_x_offset
            ref_shifted_image["CRPIX2"] += median_y_offset

            ref_shifted_path = ref_path.parent / (ref_path.stem + "_shifted.fits")
            ref_shifted_image[BASE_NAME_KEY] = ref_shifted_path.name
            logger.debug(f"Saving shifted reference image to {ref_shifted_path}")
            self.save_fits(ref_shifted_image, ref_shifted_path)

            ref_shifted_image = self.open_fits(ref_shifted_path)

            # Run sextractor and psfex and sextractor again on the shifted reference
            ref_sextractor = self.sextractor(
                output_sub_dir=self.temp_output_subtract_dir,
            )
            ref_sextractor.set_night(night_sub_dir=self.night_sub_dir)

            ref_shifted_img = ref_sextractor.apply(ImageBatch(ref_shifted_image))[0]

            ref_psfex = self.psfex(
                output_sub_dir=self.temp_output_subtract_dir, norm_fits=True
            )

            ref_shifted_img = ref_psfex.apply(ImageBatch(ref_shifted_img))[0]

            logger.debug(f"Running photometry on " f"{ref_shifted_img.get_name()}")

            # Run Sextractor again using PSFex model
            ref_psf_phot_sextractor = self.phot_sextractor(
                output_sub_dir=self.temp_output_subtract_dir,
            )
            ref_psf_phot_sextractor.set_night(night_sub_dir=self.night_sub_dir)

            final_ref_shifted_image = ref_psf_phot_sextractor.apply(
                ImageBatch(ref_shifted_img)
            )[0]

            # Save the final resampled, sextracted and psfexed reference image
            self.save_fits(final_ref_shifted_image, ref_shifted_path)

            # Update the header of the science image with the path to the shifted
            # reference image

            image[REF_IMG_KEY] = ref_shifted_path.as_posix()

        return batch
