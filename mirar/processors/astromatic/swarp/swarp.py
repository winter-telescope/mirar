"""
Module relating to `swarp <https://www.astromatic.net/software/swarp`_
"""
import logging
from pathlib import Path
from typing import Literal

import numpy as np
from astropy.wcs import WCS

from mirar.data import Image, ImageBatch
from mirar.errors import ProcessorError
from mirar.paths import (
    BASE_NAME_KEY,
    LATEST_SAVE_KEY,
    LATEST_WEIGHT_SAVE_KEY,
    RAW_IMG_KEY,
    STACKED_COMPONENT_IMAGES_KEY,
    SWARP_FLUX_SCALING_KEY,
    all_astrometric_keywords,
    copy_temp_file,
    get_output_dir,
    get_temp_path,
)
from mirar.processors.astromatic.scamp.scamp import scamp_header_key
from mirar.processors.astromatic.swarp.swarp_wrapper import run_swarp
from mirar.processors.base_processor import BaseImageProcessor

logger = logging.getLogger(__name__)


class SwarpError(ProcessorError):
    """Error relating to swarp"""


class Swarp(BaseImageProcessor):
    """
    Processor to apply Swarp
    """

    base_key = "swarp"

    def __init__(
        self,
        swarp_config_path: str | Path,
        temp_output_sub_dir: str = "swarp",
        pixscale: float | None = None,
        x_imgpixsize: float | None = None,
        y_imgpixsize: float | None = None,
        propogate_headerlist: list[str] | None = None,
        center_type: Literal["MOST", "ALL", "MANUAL"] = "MANUAL",
        center_ra: float | None = None,
        center_dec: float | None = None,
        gain: float | None = None,
        include_scamp: bool = True,
        cache: bool = False,
        subtract_bkg: bool = False,
        flux_scaling_factor: float | None = None,
        calculate_dims_in_swarp: bool = False,
    ):
        """

        Args:
            swarp_config_path: str
                path to config path
            temp_output_sub_dir: str
                output sub-directory
            pixscale: float
                Pixel scale in degrees. If None, set as median of the pixel scales
                 of input images.
            x_imgpixsize: float
                X-dimension in pixels. If None, set as max x-size of input images. If
                you want a stacked image covering all input images, set
                calculate_dims_in_swarp to True instead.
            y_imgpixsize: float
                Y-dimension in pixels. If None, set as max y-size of input images.
            propogate_headerlist: list
                List of header keywords to propagate. Recommended to leave None, the
                processor will take care of it.
            center_ra: float
                Desired central RA of output image. If None, set as the median of the
                input images.
            center_dec:
                Desired central Dec of output image. If None, set as the median of the
                input images.
            gain: float
                Gain
            include_scamp: bool
                Whether to include scamp results or not?
            cache: bool
                Save temporary files?
            subtract_bkg:
                Subtract background?
            flux_scaling_factor:
                Do you want to scale the images by some factor? There are two ways to
                provide this value. First, you can specify the scaling factor for the
                image in the header using the FLXSCALE keyword. This is the preferred
                way for a batch with multiple images. This calculation should happen
                somewhere upstream. The second option, which could be used for cases
                for single images is by providing a float value to this parameter.
                The processor will then create a FLXSCALE card in the header with this
                value. By design, Swarp is forced to look for FLXSCALE in the header
                for scaling, so this keyword should not be used for other purposes.
                The FLXSCALE keyword is not propagated to the resampled/stacked image to
                avoid unwanted re-scalings for future resamplings of the image.
                If you want this logged, you should add it manually to the headers where
                you calculate it using some keyword other than FLXSCALE. If no value is
                provided, it defaults to 1 (no scaling). If multiple images are provided
                without FLXSCALE in their headers, the same scaling will be applied to
                all of them.
        """
        super().__init__()
        self.swarp_config_path = Path(swarp_config_path)
        assert (
            self.swarp_config_path.exists()
        ), f"Swarp config path {self.swarp_config_path} does not exist."
        self.temp_output_sub_dir = temp_output_sub_dir
        self.pixscale = pixscale
        self.propogate_headerlist = propogate_headerlist
        self.center_ra = center_ra
        self.center_dec = center_dec
        self.x_imgpixsize = x_imgpixsize
        self.y_imgpixsize = y_imgpixsize
        self.include_scamp = include_scamp
        self.cache = cache
        self.gain = gain
        self.subtract_bkg = subtract_bkg
        self.flux_scaling_factor = flux_scaling_factor
        self.calculate_dims_in_swarp = calculate_dims_in_swarp
        self.center_type = center_type

    def __str__(self) -> str:
        return "Processor to apply swarp to images, stacking them together."

    def get_swarp_output_dir(self) -> Path:
        """
        Get custom output directory for swarp

        :return: Swarp directory
        """
        return get_output_dir(self.temp_output_sub_dir, self.night_sub_dir)

    def _apply_to_images(
        self,
        batch: ImageBatch,
    ) -> ImageBatch:
        swarp_output_dir = self.get_swarp_output_dir()
        swarp_output_dir.mkdir(parents=True, exist_ok=True)

        swarp_image_list_path = swarp_output_dir.joinpath(
            Path(batch[0][BASE_NAME_KEY]).name + "_swarp_img_list.txt",
        )
        swarp_weight_list_path = swarp_output_dir.joinpath(
            Path(batch[0][BASE_NAME_KEY]).name + "_swarp_weight_list.txt",
        )
        logger.debug(f"Writing file list to {swarp_image_list_path}")

        temp_files = [swarp_image_list_path, swarp_weight_list_path]

        output_image_path = swarp_output_dir.joinpath(
            Path(batch[0][BASE_NAME_KEY]).name + "_stack.fits",
        )

        logger.debug(f"Saving to {output_image_path}")

        all_pixscales = []
        all_imgpixsizes = []
        all_ras = []
        all_decs = []

        pixscale_to_use = self.pixscale
        x_imgpixsize_to_use = self.x_imgpixsize
        y_imgpixsize_to_use = self.y_imgpixsize
        center_ra_to_use = self.center_ra
        center_dec_to_use = self.center_dec

        missing_scale_bool = (
            (pixscale_to_use is None)
            | (x_imgpixsize_to_use is None)
            | (y_imgpixsize_to_use is None)
            | (center_ra_to_use is None)
            | (center_ra_to_use is None)
        )

        with open(swarp_image_list_path, "w", encoding="utf8") as img_list, open(
            swarp_weight_list_path, "w", encoding="utf8"
        ) as weight_list:
            for image in batch:
                if missing_scale_bool:
                    ra, dec, pixscale, imgpixsize = self.get_image_scale(image)
                    all_pixscales.append(pixscale)
                    all_imgpixsizes.append(imgpixsize)
                    all_ras.append(ra)
                    all_decs.append(dec)

                logger.debug(f"{all_ras}, {all_decs}")

                if self.include_scamp:
                    temp_head_path = copy_temp_file(
                        output_dir=swarp_output_dir,
                        file_path=Path(image[scamp_header_key]),
                    )

                if np.logical_and(
                    SWARP_FLUX_SCALING_KEY in image.header.keys(),
                    self.flux_scaling_factor is not None,
                ):
                    err = (
                        f"{SWARP_FLUX_SCALING_KEY} is present in header, and"
                        f"a value for flux_scaling_factor has also been provided. "
                        f"Please use only one."
                    )
                    raise SwarpError(err)

                if SWARP_FLUX_SCALING_KEY not in image.header.keys():
                    if self.flux_scaling_factor is None:
                        image[SWARP_FLUX_SCALING_KEY] = 1
                    else:
                        image[SWARP_FLUX_SCALING_KEY] = self.flux_scaling_factor

                temp_img_path = get_temp_path(swarp_output_dir, image[BASE_NAME_KEY])

                self.save_fits(image, temp_img_path)

                temp_mask_path = self.save_weight_image(image, temp_img_path)

                img_list.write(f"{temp_img_path}\n")
                weight_list.write(f"{temp_mask_path}\n")

                temp_files += [temp_img_path, temp_mask_path]

                if self.include_scamp:
                    temp_files += [temp_head_path]

        if pixscale_to_use is None:
            pixscale_to_use = np.max(all_pixscales)
        if x_imgpixsize_to_use is None:
            x_imgpixsize_to_use = np.max(all_imgpixsizes)
        if y_imgpixsize_to_use is None:
            y_imgpixsize_to_use = np.max(all_imgpixsizes)
        if center_ra_to_use is None:
            center_ra_to_use = np.median(all_ras)
        if center_dec_to_use is None:
            center_dec_to_use = np.median(all_decs)

        if self.calculate_dims_in_swarp:
            x_imgpixsize_to_use = None
            y_imgpixsize_to_use = None

        if self.center_type != "MANUAL":
            center_dec_to_use = None
            center_ra_to_use = None

        output_image_weight_path = output_image_path.with_suffix(".weight.fits")

        run_swarp(
            stack_list_path=swarp_image_list_path,
            swarp_config_path=self.swarp_config_path,
            out_path=output_image_path,
            weight_list_path=swarp_weight_list_path,
            weight_out_path=output_image_weight_path,
            pixscale=pixscale_to_use,
            x_imgpixsize=x_imgpixsize_to_use,
            y_imgpixsize=y_imgpixsize_to_use,
            propogate_headerlist=self.propogate_headerlist,
            center_type=self.center_type,
            center_ra=center_ra_to_use,
            center_dec=center_dec_to_use,
            gain=self.gain,
            subtract_bkg=self.subtract_bkg,
            flux_scaling_keyword=SWARP_FLUX_SCALING_KEY,
            cache=self.cache,
        )

        new_image = self.load_swarp_output(
            output_image_path, output_image_weight_path, batch
        )

        if not self.cache:
            for temp_file in temp_files:
                temp_file.unlink()
                logger.debug(f"Deleted temporary file {temp_file}")

        return ImageBatch([new_image])

    def load_swarp_output(
        self, output_image_path: Path, output_image_weight_path: Path, batch: ImageBatch
    ) -> Image:
        """
        Function to load the output of SWarp into an Image object,
        and transfer header information from the input ImageBatch

        :param output_image_path: Swarp output image path
        :param output_image_weight_path: Swarp output weight image path
        :param batch: batch of component images
        :return: Updated Image object
        """
        new_image = self.open_fits(output_image_path)
        # Add missing keywords that are common in all input images to the
        # header of resampled image, and save again
        # Omit any astrometric keywords

        for key in batch[0].keys():
            if np.any([key not in x.keys() for x in batch]):
                continue
            if np.logical_and(
                np.sum([x[key] == batch[0][key] for x in batch]) == len(batch),
                key.strip() not in all_astrometric_keywords,
            ):
                if key not in new_image.keys():
                    try:
                        new_image[key] = batch[0][key]
                    except ValueError:
                        continue
        new_image["COADDS"] = np.sum([x["COADDS"] for x in batch])

        new_image[RAW_IMG_KEY] = ",".join([x[RAW_IMG_KEY] for x in batch])

        try:
            component_images_list = [x[LATEST_SAVE_KEY] for x in batch]
            new_image[STACKED_COMPONENT_IMAGES_KEY] = ",".join(component_images_list)
        except KeyError:
            pass

        new_image[BASE_NAME_KEY] = output_image_path.name
        new_image[LATEST_WEIGHT_SAVE_KEY] = output_image_weight_path.as_posix()
        return new_image

    @staticmethod
    def get_image_scale(image: Image) -> tuple[float, float, float, float]:
        """
        Function to get the image scale from the image header

        :param image: Image to get the scale from
        :return: RA, Dec, pixel scale, image pixel size
        """

        wcs = WCS(image.get_header())

        cd11 = image["CD1_1"]
        cd21 = image["CD2_1"]

        nxpix = image["NAXIS1"]
        nypix = image["NAXIS2"]

        image_x_cen = nxpix / 2
        image_y_cen = nypix / 2

        [ra, dec] = wcs.all_pix2world(image_x_cen, image_y_cen, 1)

        xscale = np.sqrt(cd11**2 + cd21**2)

        pixscale = xscale * 3600
        imgpixsize = max(nxpix, nypix)

        return ra, dec, pixscale, imgpixsize
