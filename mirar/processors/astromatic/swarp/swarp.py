"""
Module relating to `swarp <https://www.astromatic.net/software/swarp`_
"""
import logging
import os
from pathlib import Path
from typing import Optional

import numpy as np
from astropy.wcs import WCS

from mirar.data import ImageBatch
from mirar.errors import ProcessorError
from mirar.paths import (
    BASE_NAME_KEY,
    LATEST_WEIGHT_SAVE_KEY,
    RAW_IMG_KEY,
    SWARP_FLUX_SCALING_KEY,
    all_astrometric_keywords,
    copy_temp_file,
    get_output_dir,
    get_temp_path,
)
from mirar.processors.astromatic.scamp.scamp import scamp_header_key
from mirar.processors.base_processor import BaseImageProcessor
from mirar.utils import execute

logger = logging.getLogger(__name__)


class SwarpError(ProcessorError):
    """Error relating to swarp"""


def run_swarp(
    stack_list_path: str | Path,
    swarp_config_path: str | Path,
    out_path: str | Path,
    weight_list_path: Optional[str | Path] = None,
    weight_out_path: Optional[str | Path] = None,
    pixscale: Optional[float] = None,
    x_imgpixsize: Optional[float] = None,
    y_imgpixsize: Optional[float] = None,
    propogate_headerlist: Optional[list] = None,
    center_ra: Optional[float] = None,
    center_dec: Optional[float] = None,
    combine: bool = True,
    gain: Optional[float] = None,
    subtract_bkg: bool = False,
    flux_scaling_keyword: str = None,
    cache: bool = False,
):
    """
    Wrapper to resample and stack images with swarp

    Parameters
    ----------
    stack_list_path : string
        Name of file containing the names of files to be stacked
        One file name per line
    weight_list_path : string
        Name of file containing the names of weight files to be stacked
        One file name per line
    swarp_config_path: str
        Path of Swarp config file
    out_path : string
        Path of stacked output file
    weight_out_path: str
        Path of output weight image
    pixscale: float
        Pixelscale in degrees
    x_imgpixsize: float
        X-dimension in pixels
    y_imgpixsize: float
        Y-dimension in pixels
    propogate_headerlist: list
        Headerlist to propagate from header
    center_ra: float
        Central RA
    center_dec: float
        Central Dec
    combine: bool
        Combine and coadd all images? For reasons internal to Swarp, it is strongly
        advised to always set this to True (even if you are running Swarp on only
        one image).
    gain: float
        Gain
    subtract_bkg: bool
        Background subtraction
    flux_scaling_keyword: str
        What flux scaling keyword do you want to use? If None, the default value in
        the config will be used
    """

    swarp_command = (
        f"swarp -c {swarp_config_path} "
        f"@{stack_list_path} "
        f"-IMAGEOUT_NAME {out_path} "
        f"-RESAMPLE Y -RESAMPLE_DIR {os.path.dirname(out_path)} "
    )

    if subtract_bkg:
        swarp_command += "-SUBTRACT_BACK Y "
    else:
        swarp_command += "-SUBTRACT_BACK N "
    if combine:
        swarp_command += "-COMBINE Y -COMBINE_TYPE MEDIAN "
    else:
        swarp_command += "-COMBINE N "

    if weight_list_path is not None:
        swarp_command += f" -WEIGHT_TYPE MAP_WEIGHT -WEIGHT_IMAGE @{weight_list_path} "

    if weight_out_path is not None:
        swarp_command += f" -WEIGHTOUT_NAME {weight_out_path}"

    if pixscale is not None:
        swarp_command += f" -PIXELSCALE_TYPE MANUAL -PIXEL_SCALE {pixscale}"

    if propogate_headerlist is not None:
        swarp_command += " -COPY_KEYWORDS "
        for keyword in propogate_headerlist:
            swarp_command += f"{keyword},"

        # remove final comma
        swarp_command = swarp_command[:-1]

    if np.logical_and(center_ra is not None, center_dec is not None):
        swarp_command += f" -CENTER_TYPE MANUAL -CENTER {center_ra},{center_dec}"

    if x_imgpixsize is not None:
        swarp_command += f" -IMAGE_SIZE {x_imgpixsize}"
        if y_imgpixsize is not None:
            swarp_command += f",{y_imgpixsize}"

    if gain is not None:
        swarp_command += f" -GAIN {gain}"

    if flux_scaling_keyword is not None:
        swarp_command += f" -FSCALE_KEYWORD {flux_scaling_keyword}"

    if not cache:
        swarp_command += " -DELETE_TMPFILES Y"
    else:
        swarp_command += " -DELETE_TMPFILES N"

    execute(swarp_command)


class Swarp(BaseImageProcessor):
    """
    Processor to apply Swarp
    """

    base_key = "swarp"

    def __init__(
        self,
        swarp_config_path: str,
        temp_output_sub_dir: str = "swarp",
        pixscale: Optional[float] = None,
        x_imgpixsize: Optional[float] = None,
        y_imgpixsize: Optional[float] = None,
        propogate_headerlist: Optional[list] = None,
        center_ra: Optional[float] = None,
        center_dec: Optional[float] = None,
        gain: Optional[float] = None,
        include_scamp: bool = True,
        combine: bool = True,
        cache: bool = False,
        subtract_bkg: bool = False,
        flux_scaling_factor: float = None,
        calculate_dims_in_swarp: bool = False,
    ):
        """

        Args:
            swarp_config_path: str
                path to config path
            temp_output_sub_dir: str
                output sub-directory
            pixscale: float
                Pixel scale in degrees
            x_imgpixsize: float
                X-dimension in pixels
            y_imgpixsize: float
                Y-dimension in pixels
            propogate_headerlist: list
                List of header keywords to propagate. Recommended to leave None, the
                processor will take care of it.
            center_ra: float
                Desired central RA of output image
            center_dec:
                Desired central Dec of output image
            gain: float
                Gain
            include_scamp: bool
                Whether to include scamp results or not?
            combine: bool
                Combine and coadd all images? For reasons internal to Swarp, it is
                strongly advised to always set this to True (even if you are running
                Swarp on only one image). The processor will raise an
                error if you try running it on a Batch with multiple images by
                setting combine to False. If you want to resample multiple images,
                just DeBatch them and run Swarp separately.
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
        self.swarp_config = swarp_config_path
        self.temp_output_sub_dir = temp_output_sub_dir
        self.pixscale = pixscale
        self.propogate_headerlist = propogate_headerlist
        self.center_ra = center_ra
        self.center_dec = center_dec
        self.x_imgpixsize = x_imgpixsize
        self.y_imgpixsize = y_imgpixsize
        self.include_scamp = include_scamp
        self.combine = combine
        self.cache = cache
        self.gain = gain
        self.subtract_bkg = subtract_bkg
        self.flux_scaling_factor = flux_scaling_factor
        self.calculate_dims_in_swarp = calculate_dims_in_swarp

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

        # If swarp is run with combine -N option,
        # it just sets the output name as inpname+.resamp.fits

        if self.combine:
            output_image_path = swarp_output_dir.joinpath(
                Path(batch[0][BASE_NAME_KEY]).name + "_stack.fits",
            )
        else:
            if len(batch) > 1:
                err = (
                    f"Attempting to run Swarp with batch-size of length "
                    f"{len(batch)} but with combine=False. Either set combine=True"
                    f", or DeBatch the images"
                )
                logger.error(err)
                raise SwarpError(err)

            logger.warning(
                "You are choosing to run Swarp without combining the image. "
                "This causes swarp to output an intermediate image, "
                "with possibly incorrect FLXSCALE values. Please consider "
                "running with combine=True, it almost always gives the same "
                "result as running it without, but will have consistent "
                "headers."
            )
            output_image_path = swarp_output_dir.joinpath(
                Path(batch[0][BASE_NAME_KEY]).with_suffix(".resamp.fits").name,
            )
        logger.debug(f"Saving to {output_image_path}")

        all_pixscales = []
        all_imgpixsizes = []
        all_ras = []
        all_decs = []

        with open(swarp_image_list_path, "w", encoding="utf8") as img_list, open(
            swarp_weight_list_path, "w", encoding="utf8"
        ) as weight_list:
            for image in batch:
                pixscale_to_use = self.pixscale
                x_imgpixsize_to_use = self.x_imgpixsize
                y_imgpixsize_to_use = self.y_imgpixsize
                center_ra_to_use = self.center_ra
                center_dec_to_use = self.center_dec

                if (
                    (self.pixscale is None)
                    | (self.x_imgpixsize is None)
                    | (self.y_imgpixsize is None)
                    | (self.center_ra is None)
                    | (self.center_dec is None)
                ):
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

        output_image_weight_path = output_image_path.with_suffix(".weight.fits")

        run_swarp(
            stack_list_path=swarp_image_list_path,
            swarp_config_path=self.swarp_config,
            out_path=output_image_path,
            weight_list_path=swarp_weight_list_path,
            weight_out_path=output_image_weight_path,
            pixscale=pixscale_to_use,
            x_imgpixsize=x_imgpixsize_to_use,
            y_imgpixsize=y_imgpixsize_to_use,
            propogate_headerlist=self.propogate_headerlist,
            center_ra=center_ra_to_use,
            center_dec=center_dec_to_use,
            combine=self.combine,
            gain=self.gain,
            subtract_bkg=self.subtract_bkg,
            flux_scaling_keyword=SWARP_FLUX_SCALING_KEY,
            cache=self.cache,
        )

        # Check if output image exists if combine is no.
        # This is the intermediate image that swarp makes
        # Hopefully this is obsolete now and noone uses this
        if not self.combine:
            temp_output_image_path = get_temp_path(
                swarp_output_dir,
                output_image_path.name,
            )

            temp_output_image_weight_path = temp_output_image_path.with_suffix(
                ".weight.fits"
            )

            if temp_output_image_path.exists():
                temp_output_image_path.rename(output_image_path)
                temp_output_image_weight_path.rename(output_image_weight_path)
            else:
                err = (
                    f"Swarp seems to have misbehaved, "
                    f"and not made the correct output file {temp_output_image_path}"
                )
                logger.error(err)
                raise SwarpError(err)

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
                    new_image[key] = batch[0][key]
        new_image["COADDS"] = np.sum([x["COADDS"] for x in batch])

        new_image[RAW_IMG_KEY] = ",".join([x[RAW_IMG_KEY] for x in batch])
        new_image[BASE_NAME_KEY] = output_image_path.name
        new_image[LATEST_WEIGHT_SAVE_KEY] = output_image_weight_path.as_posix()
        self.save_fits(new_image, output_image_path)
        logger.info(f"Saved resampled image to {output_image_path.name}")

        if not self.cache:
            for temp_file in temp_files:
                temp_file.unlink()
                logger.debug(f"Deleted temporary file {temp_file}")

        return ImageBatch([new_image])
