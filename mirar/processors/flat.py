"""
Module containing processors for flat calibration
"""

import logging
import os.path
import sys
from collections.abc import Callable
from copy import copy

import numpy as np
from astropy.convolution import Tophat2DKernel, convolve_fft

from mirar.data import Image, ImageBatch
from mirar.errors import ImageNotFoundError
from mirar.paths import (
    BASE_NAME_KEY,
    COADD_KEY,
    EXPTIME_KEY,
    FLAT_FRAME_KEY,
    LATEST_SAVE_KEY,
    OBSCLASS_KEY,
)
from mirar.processors.base_processor import ProcessorPremadeCache, ProcessorWithCache
from mirar.processors.utils.image_selector import select_from_images

logger = logging.getLogger(__name__)


def get_convolution(data: np.ndarray, kernel_width: int) -> np.ndarray:
    """
    Convolve data with a tophat kernel

    :param data: Image data
    :param kernel_width: Width of the kernel (pixels)
    :return: Smoothed image
    """
    pad_top = np.array([data[0] for _ in range(kernel_width)])
    pad_bottom = np.array([data[-1] for _ in range(kernel_width)])
    extended = np.vstack([pad_top, data, pad_bottom])
    pad_left = np.array([extended.T[0] for _ in range(kernel_width)])
    pad_right = np.array([extended.T[-1] for _ in range(kernel_width)])
    extended = np.hstack([pad_left.T, extended, pad_right.T])

    tophat_kernel = Tophat2DKernel(kernel_width)
    smooth_illumination = convolve_fft(
        extended, tophat_kernel, nan_treatment="interpolate"
    )[kernel_width:-kernel_width, kernel_width:-kernel_width]
    return smooth_illumination


import warnings

from astropy.stats import sigma_clipped_stats


def get_outlier_pixel_mask(img: np.ndarray, thresh: float = 3.0) -> np.ndarray:
    """
    Get oulier pixels that are above or below a threshold
    :param img: np.ndarray
    :param thresh: float
    """
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        _, median, std = sigma_clipped_stats(img, sigma=3.0)
    return (img < (median - thresh * std)) | (img > (median + thresh * std))


def construct_smooth_gradient_for_image(data: np.ndarray) -> np.ndarray:
    """
    Construct a smooth gradient for the image
    :param data: np.ndarray
    :return: np.ndarray
    """
    smooth_img = get_convolution(data, 100)
    return smooth_img


def smooth_and_normalize_image(data: np.ndarray) -> np.ndarray:
    """
    Smooth and normalize the image
    :param data: np.ndarray
    :return: np.ndarray
    """
    smooth_img = get_convolution(data, 100)
    return data / smooth_img


def get_smoothened_outlier_pixel_mask_from_list(
    img_data_list: list[np.ndarray], threshold: float = 3.0
) -> np.ndarray:
    """
    Take a list of images and return a mask of outlier pixels after removing a smooth
    gradient from them
    :param img_data_list: list[np.ndarray]
    :param threshold: float
    :return: np.ndarray
    """
    masks = []
    for data in img_data_list:
        smooth_img = smooth_and_normalize_image(data)
        mask = get_outlier_pixel_mask(smooth_img, thresh=threshold)
        masks.append(mask)
    return np.logical_and.reduce(masks)


class MissingFlatError(ImageNotFoundError):
    """
    Error for when a dark image is missing
    """


def default_select_flat(
    images: ImageBatch,
) -> ImageBatch:
    """
    Select images tagged as flat

    :param images: set of images
    :return: subset of flat images
    """
    return select_from_images(images, key=OBSCLASS_KEY, target_values="flat")


class FlatCalibrator(ProcessorWithCache):
    """
    Processor to apply flat calibration
    """

    base_key = "flat"

    def __init__(
        self,
        *args,
        x_min: int = 0,
        x_max: int = sys.maxsize,
        y_min: int = 0,
        y_max: int = sys.maxsize,
        flat_nan_threshold: float = 0.0,
        select_flat_images: Callable[[ImageBatch], ImageBatch] = default_select_flat,
        flat_mask_key: str = None,
        flat_mode: str = "median",
        **kwargs,
    ):
        super().__init__(*args, **kwargs)
        self.x_min = x_min
        self.x_max = x_max
        self.y_min = y_min
        self.y_max = y_max
        self.flat_nan_threshold = flat_nan_threshold
        self.select_cache_images = select_flat_images
        self.flat_mask_key = flat_mask_key
        self.flat_mode = flat_mode
        if self.flat_mode not in ["median", "pixel", "structure"]:
            raise ValueError(f"Flat mode {self.flat_mode} not supported")

    def description(self) -> str:
        return "Creates a flat image, divides other images by this image."

    def _apply_to_images(
        self,
        batch: ImageBatch,
    ) -> ImageBatch:
        master_flat = self.get_cache_file(batch)
        master_flat_data = master_flat.get_data()

        mask = master_flat_data <= self.flat_nan_threshold

        if np.sum(mask) > 0:
            master_flat_data[mask] = np.nan

        for image in batch:
            data = image.get_data()
            data = data / master_flat_data
            image.set_data(data)
            image[FLAT_FRAME_KEY] = master_flat[LATEST_SAVE_KEY]

        return batch

    def make_image(
        self,
        images: ImageBatch,
    ) -> Image:
        images = self.select_cache_images(images)

        logger.debug(f"Found {len(images)} suitable flats in batch")

        n_frames = len(images)
        if n_frames == 0:
            err = f"Found {n_frames} suitable flats in batch"
            logger.error(err)
            raise MissingFlatError(err)

        nx, ny = images[0].get_data().shape

        flats = np.zeros((nx, ny, n_frames))

        flat_exptimes = []
        for i, img in enumerate(images):
            data = img.get_data().copy()

            if self.flat_mask_key is not None:
                if self.flat_mask_key not in img.header.keys():
                    err = (
                        f"Image {img} does not have a mask with key "
                        f"{self.flat_mask_key}"
                    )
                    logger.error(err)
                    raise KeyError(err)

                mask_file = img[self.flat_mask_key]
                logger.debug(f"Masking flat {img[BASE_NAME_KEY]} with mask {mask_file}")
                if not os.path.exists(mask_file):
                    err = f"Mask file {mask_file} does not exist"
                    logger.error(err)
                    raise FileNotFoundError(err)

                mask_img = self.open_fits(mask_file)
                pixels_to_keep = mask_img.get_data().astype(bool)
                mask = ~pixels_to_keep
                logger.debug(
                    f"Masking {np.sum(mask)} pixels in flat {img[BASE_NAME_KEY]}"
                )
                data[mask] = np.nan

            flat_exptimes.append(img[EXPTIME_KEY])

            median = np.nanmedian(
                data[self.x_min : self.x_max, self.y_min : self.y_max]
            )

            flats[:, :, i] = data / median

        logger.debug(f"Median combining {n_frames} flats")

        master_flat = np.nanmedian(flats, axis=2)

        if self.flat_mode != "median":

            if self.flat_mode == "pixel":

                mask = get_smoothened_outlier_pixel_mask_from_list(
                    [x.get_data() for x in images]
                )

                # flatdata_norm_smooth = get_convolution(master_flat, 100)
                #
                # pixel_variation = master_flat / flatdata_norm_smooth
                #
                # # Clip outliers (they'll get worked out in stacking)
                # std = np.nanstd(pixel_variation)
                # sig = abs(pixel_variation - np.nanmedian(pixel_variation)) / std
                #
                # mask = sig > 1.0

                # mask = get_outlier_pixel_mask(master_flat, thresh=1.5)

                frac = np.sum(mask) / len(mask.flatten())

                logger.info(  # FIXME: Change to debug
                    f"Masking {100.*frac:.1f}% of pixels in image"
                )
                # pixel_variation[mask] = np.nan
                # master_flat = pixel_variation / np.nanmedian(pixel_variation)
                master_flat = np.ones_like(master_flat)
                master_flat[mask] = np.nan

            elif self.flat_mode == "structure":
                flatdata_norm_smooth = get_convolution(master_flat, 100)
                flatdata_norm_smooth[np.isnan(master_flat)] = np.nan

                pixel_variation = master_flat / flatdata_norm_smooth

                # Clip outliers (they'll get worked out in stacking)
                std = np.nanstd(pixel_variation)
                sig = abs(pixel_variation - np.nanmedian(pixel_variation)) / std

                mask = sig > 1.0

                logger.debug(
                    f"Masking {np.sum(mask)} pixels "
                    f"out of {len(mask.flatten()) }in flat"
                )

                master_flat = np.ones_like(master_flat)
                master_flat[mask] = np.nan

                # pixel_variation[mask] = np.nan
                # master_flat = (
                #     pixel_variation
                #     * flatdata_norm_smooth
                #     / np.nanmedian(pixel_variation)
                # )

            else:
                raise ValueError(f"Flat mode {self.flat_mode} not supported")

        master_flat_image = Image(master_flat, header=copy(images[0].get_header()))
        master_flat_image[COADD_KEY] = n_frames

        master_flat_image["INDIVEXP"] = ",".join(
            [str(x) for x in np.unique(flat_exptimes)]
        )
        return master_flat_image


class SkyFlatCalibrator(FlatCalibrator):
    """
    Processor to do flat calibration using sky flats
    """

    def __init__(self, *args, flat_mask_key=None, **kwargs):
        super().__init__(
            *args,
            select_flat_images=self.select_sky_flat,
            flat_mask_key=flat_mask_key,
            **kwargs,
        )

    @staticmethod
    def select_sky_flat(
        images: ImageBatch,
    ) -> ImageBatch:
        """
        Select science images to use as sky flats

        :param images: set of images
        :return: subset of 'sky' images
        """
        return select_from_images(images, key=OBSCLASS_KEY, target_values="science")

    def description(self) -> str:
        return (
            "Processor to create a sky flat image, divides other images by this image."
        )


class MasterFlatCalibrator(ProcessorPremadeCache, FlatCalibrator):
    """Processor to do flat calibration with a master flate"""
