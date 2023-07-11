import logging
import os.path
from collections.abc import Callable
from pathlib import Path

import numpy as np
from astropy.wcs import WCS

from mirar.data import Image, ImageBatch
from mirar.paths import (
    BASE_NAME_KEY,
    NORM_PSFEX_KEY,
    RAW_IMG_KEY,
    REF_IMG_KEY,
    REF_PSF_KEY,
    get_output_dir,
)
from mirar.processors.astromatic.psfex.psfex import PSFex
from mirar.processors.astromatic.sextractor.sextractor import Sextractor
from mirar.processors.astromatic.swarp.swarp import Swarp
from mirar.processors.base_processor import BaseImageProcessor
from mirar.references.base_reference_generator import BaseReferenceGenerator

logger = logging.getLogger(__name__)


class ProcessReference(BaseImageProcessor):
    base_key = "REFPREP"

    def __init__(
        self,
        ref_image_generator: Callable[..., BaseReferenceGenerator],
        swarp_resampler: Callable[..., Swarp],
        sextractor: Callable[..., Sextractor],
        ref_psfex: Callable[..., PSFex],
        temp_output_subtract_dir: str = "subtract",
    ):
        super().__init__()
        self.ref_image_generator = ref_image_generator
        self.swarp_resampler = swarp_resampler
        self.sextractor = sextractor
        self.psfex = ref_psfex
        self.temp_output_subtract_dir = temp_output_subtract_dir

    def get_sub_output_dir(self):
        return get_output_dir(self.temp_output_subtract_dir, self.night_sub_dir)

    @staticmethod
    def get_image_header_params(image: Image):
        header = image.get_header()

        wcs = WCS(header)

        image_x_cent = header["NAXIS1"] / 2
        image_y_cent = header["NAXIS2"] / 2
        x_imgsize = header["NAXIS1"]
        y_imgsize = header["NAXIS2"]

        [ra_cent, dec_cent] = wcs.all_pix2world(image_x_cent, image_y_cent, 1)

        pixscale = np.abs(header["CD1_1"]) * 3600
        gain = header["GAIN"]

        return (
            image_x_cent,
            image_y_cent,
            ra_cent,
            dec_cent,
            pixscale,
            x_imgsize,
            y_imgsize,
            gain,
        )

    def _apply_to_images(
        self,
        batch: ImageBatch,
    ) -> ImageBatch:
        try:
            os.makedirs(self.get_sub_output_dir())
        except OSError:
            pass

        new_batch = ImageBatch()

        for image in batch:
            ref_writer = self.ref_image_generator(image)

            ref_image_path = ref_writer.write_reference(
                image, output_dir=self.get_sub_output_dir()
            )

            ref_image = self.open_fits(ref_image_path)

            if BASE_NAME_KEY not in ref_image.keys():
                base_name = Path(ref_image_path).name
                logger.debug(f"Setting basename to {base_name})")
                ref_image[BASE_NAME_KEY] = base_name

            if RAW_IMG_KEY not in ref_image.keys():
                ref_image[RAW_IMG_KEY] = str(ref_image_path)

            ref_gain = ref_image["GAIN"]

            sci_x_cent, sci_y_cent = image["NAXIS1"] / 2, image["NAXIS2"] / 2

            header = image.get_header()
            wcs = WCS(header)

            [sci_ra_cent, sci_dec_cent] = wcs.all_pix2world(sci_x_cent, sci_y_cent, 1)
            sci_pixscale = np.abs(image["CD1_1"]) * 3600

            sci_x_imgsize = image["NAXIS1"]
            sci_y_imgsize = image["NAXIS2"]

            sci_gain = image["GAIN"]

            propogate_headerlist = ["TMC_ZP", "TMC_ZPSD"]

            # Resample ref image onto science image
            ref_resampler = self.swarp_resampler(
                pixscale=sci_pixscale,
                x_imgpixsize=sci_x_imgsize,
                y_imgpixsize=sci_y_imgsize,
                center_ra=sci_ra_cent,
                center_dec=sci_dec_cent,
                propogate_headerlist=propogate_headerlist,
                temp_output_sub_dir=self.temp_output_subtract_dir,
                include_scamp=False,
                combine=False,
                gain=ref_gain,
                subtract_bkg=True,
            )

            ref_resampler.set_night(night_sub_dir=self.night_sub_dir)
            resampled_ref_img = ref_resampler.apply(ImageBatch(ref_image))[0]

            resampled_ref_path = os.path.join(
                self.get_sub_output_dir(), resampled_ref_img.get_name()
            )
            self.save_fits(resampled_ref_img, resampled_ref_path)

            (
                _,
                _,
                ref_resamp_ra_cent,
                ref_resamp_dec_cent,
                ref_resamp_pixscale,
                ref_resamp_x_imgsize,
                ref_resamp_y_imgsize,
                ref_resamp_gain,
            ) = self.get_image_header_params(resampled_ref_img)

            # This is a fall back if the ref image resampling by Swarp fails
            sci_resampler = self.swarp_resampler(
                pixscale=ref_resamp_pixscale,
                x_imgpixsize=ref_resamp_x_imgsize,
                y_imgpixsize=ref_resamp_y_imgsize,
                center_ra=ref_resamp_ra_cent,
                center_dec=ref_resamp_dec_cent,
                propogate_headerlist=propogate_headerlist,
                temp_output_sub_dir=self.temp_output_subtract_dir,
                include_scamp=False,
                combine=False,
                gain=sci_gain,
                subtract_bkg=True,
            )
            sci_resampler.set_night(night_sub_dir=self.night_sub_dir)

            resampled_sci_image = sci_resampler.apply(ImageBatch([image]))[0]

            # Detect source in reference image, and save as catalog

            ref_sextractor = self.sextractor(
                output_sub_dir=self.temp_output_subtract_dir, gain=ref_resamp_gain
            )
            ref_sextractor.set_night(night_sub_dir=self.night_sub_dir)

            resampled_ref_sextractor_img = ref_sextractor.apply(
                ImageBatch(resampled_ref_img)
            )[0]

            rrsi_path = os.path.join(
                self.get_sub_output_dir(), resampled_ref_sextractor_img.get_name()
            )

            self.save_fits(image=resampled_ref_sextractor_img, path=rrsi_path)
            logger.debug(f"Saved reference image to {rrsi_path}")

            ref_psfex = self.psfex(
                output_sub_dir=self.temp_output_subtract_dir, norm_fits=True
            )

            resampled_ref_sextractor_img = ref_psfex.apply(
                ImageBatch(resampled_ref_sextractor_img)
            )[0]

            # Copy over header keys from ref to sci
            resampled_sci_image[REF_PSF_KEY] = resampled_ref_sextractor_img[
                NORM_PSFEX_KEY
            ]
            resampled_sci_image[REF_IMG_KEY] = resampled_ref_sextractor_img[
                BASE_NAME_KEY
            ]

            new_batch.append(resampled_sci_image)
        return new_batch


class GetReferenceImage(BaseImageProcessor):
    base_key = "refimg_returner"

    def __init__(
        self,
        ref_image_generator: Callable[..., BaseReferenceGenerator],
        output_sub_dir: str = "ref",
    ):
        super().__init__()
        self.ref_image_generator = ref_image_generator
        self.output_sub_dir = output_sub_dir

    def __str__(self):
        output_sub_dir = get_output_dir(
            dir_root=self.output_sub_dir, sub_dir=self.night_sub_dir
        )
        return (
            f"Processor to get reference images and " f"save them to {output_sub_dir}"
        )

    def _apply_to_images(
        self,
        batch: ImageBatch,
    ) -> ImageBatch:
        ref_batch = ImageBatch()
        for image in batch:
            logger.info(
                f"Output directory for reference building is " f"{self.night_sub_dir}"
            )
            ref_generator = self.ref_image_generator(image)

            output_sub_dir = get_output_dir(
                dir_root=self.output_sub_dir, sub_dir=self.night_sub_dir
            )
            if not output_sub_dir.exists():
                output_sub_dir.mkdir(parents=True, exist_ok=True)

            logger.info(f"Output directory for reference building is {output_sub_dir}")
            ref_image_path = ref_generator.write_reference(
                image,
                output_dir=output_sub_dir.as_posix(),
            )

            ref_image = self.open_fits(ref_image_path)

            ref_batch.append(ref_image)

        return ref_batch
