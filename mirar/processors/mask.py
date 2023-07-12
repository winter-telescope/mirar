"""
Module containing processors which mask pixels
"""
import logging
from pathlib import Path
from typing import Callable

import numpy as np
from astropy.coordinates import SkyCoord
from astropy.io import fits
from astropy.wcs import WCS

from mirar.data import Image, ImageBatch
from mirar.paths import BASE_NAME_KEY, FITS_MASK_KEY, get_output_dir
from mirar.processors.base_processor import BaseImageProcessor

logger = logging.getLogger(__name__)

MASK_VALUE = np.nan


class BaseMask(BaseImageProcessor):
    """
    Base class for masking processors
    """

    def __init__(
        self,
        write_masked_pixels_to_file: bool = False,
        output_dir: str | Path = "mask",
        only_write_mask: bool = False,
    ):
        super().__init__()
        self.write_masked_pixels_to_file = write_masked_pixels_to_file
        self.output_dir = output_dir
        self.only_write_mask = only_write_mask

    def get_mask(self, image) -> np.ndarray[bool]:
        """
        Function to get the mask for a given image
        """
        raise NotImplementedError

    def _apply_to_images(
        self,
        batch: ImageBatch,
    ) -> ImageBatch:
        for image in batch:
            data = image.get_data()
            logger.debug(f"Masking {image[BASE_NAME_KEY]}")
            pixels_to_mask = self.get_mask(image)

            if not self.only_write_mask:
                data[pixels_to_mask] = MASK_VALUE
                image.set_data(data)

            logger.debug(
                f"Masked {np.sum(pixels_to_mask)}/{pixels_to_mask.size} pixels "
                f"in {image[BASE_NAME_KEY]}"
            )

            if self.write_masked_pixels_to_file:
                mask_directory = get_output_dir(self.output_dir, self.night_sub_dir)
                if not mask_directory.exists():
                    mask_directory.mkdir(parents=True)
                mask_file_path = mask_directory.joinpath(
                    image[BASE_NAME_KEY]
                ).with_suffix(".mask.fits")

                mask_data = np.ones_like(data)
                mask_data[pixels_to_mask] = 0.0

                mask_image = Image(data=mask_data, header=image.get_header())

                self.save_fits(mask_image, mask_file_path)
                image[FITS_MASK_KEY] = mask_file_path.as_posix()

        return batch


class MaskPixelsFromPath(BaseMask):
    """
    Processor to apply a mask to images using another `mask image'.
    Following the general mirar convention, every zero pixel in the
    mask image will be masked in the science image.
    """

    base_key = "maskfrompath"

    def __init__(
        self,
        mask_path: str | Path = None,
        mask_path_key: str = None,
        write_masked_pixels_to_file: bool = False,
        output_dir: str | Path = "mask",
        only_write_mask: bool = False,
    ):
        super().__init__(
            write_masked_pixels_to_file=write_masked_pixels_to_file,
            output_dir=output_dir,
            only_write_mask=only_write_mask,
        )
        self.mask_path = mask_path
        self.mask_path_key = mask_path_key
        if mask_path is None and mask_path_key is None:
            raise ValueError("Must specify either mask_path or mask_path_key")
        if mask_path is not None and mask_path_key is not None:
            raise ValueError("Must specify either mask_path or mask_path_key, not both")

    def __str__(self) -> str:
        return f"Processor to mask bad pixels using a pre-defined map: {self.mask_path}"

    def get_mask(self, image) -> np.ndarray:
        """
        loads mask if needed, and returns it

        :return: mask
        """
        # if self.mask is None: # why is this needed?
        if self.mask_path is not None:
            mask_img = self.open_fits(self.mask_path)
        elif self.mask_path_key is not None:
            logger.debug(f"Loading mask from {image[self.mask_path_key]}")
            mask_img = self.open_fits(image[self.mask_path_key])
        else:
            raise ValueError("Must specify either mask_path or mask_path_key")
        pixels_to_keep = mask_img.get_data().astype(bool)

        return ~pixels_to_keep


class MaskPixelsFromPathInverted(MaskPixelsFromPath):
    """
    Processor to apply a mask to images using another `mask image'.
    In contrast to the general mirar convention, every non-zero pixel in the
    mask image will be masked in the science image.
    """

    base_key = "maskfrompathinverted"

    def get_mask(self, image) -> np.ndarray:
        """
        Mask pixels which are non-zero in the mask file.
        This is the inverse of MaskPixelsFromPath,
        which masks pixels which are zero in the mask file.

        :param image: image to mask
        :return: Boolean mask
        """
        mask = super().get_mask(image)
        return ~mask


class MaskPixelsFromFunction(BaseMask):
    """
    Processor to apply a mask to images using a function
    """

    base_key = "maskfromfunction"

    def __init__(
        self,
        mask_function: Callable[[Image], np.ndarray],
        write_masked_pixels_to_file: bool = False,
        output_dir: str | Path = "mask",
        only_write_mask: bool = False,
    ):
        super().__init__(
            write_masked_pixels_to_file=write_masked_pixels_to_file,
            output_dir=output_dir,
            only_write_mask=only_write_mask,
        )
        self.mask_function = mask_function

    def get_mask(self, image) -> np.ndarray:
        """
        Function to get the mask for a given image
        """
        return self.mask_function(image)


class MaskAboveThreshold(BaseMask):
    """
    Processor to mask pixels above a threshold
    """

    base_key = "maskthresh"

    def __init__(
        self,
        threshold: float = None,
        threshold_key: str = None,
        write_masked_pixels_to_file: bool = False,
        output_dir: str | Path = "mask",
        only_write_mask: bool = False,
    ):
        """
        :param threshold: threshold to mask above
        :param threshold_key: key to use to get threshold from image header
        """
        super().__init__(
            write_masked_pixels_to_file=write_masked_pixels_to_file,
            output_dir=output_dir,
            only_write_mask=only_write_mask,
        )
        self.threshold = threshold
        self.threshold_key = threshold_key
        self.write_masked_pixels_to_file = write_masked_pixels_to_file
        if threshold is None and threshold_key is None:
            raise ValueError("Must specify either threshold or threshold_key")
        if threshold is not None and threshold_key is not None:
            raise ValueError("Must specify either threshold or threshold_key, not both")

    def __str__(self) -> str:
        return f"Processor to mask pixels above a threshold: {self.threshold}"

    def get_mask(self, image) -> np.ndarray:
        """
        Returns a mask for pixels above a threshold

        :return: mask
        """
        if self.threshold is None:
            self.threshold = image.get_header()[self.threshold_key]
        pixels_to_mask = image.get_data() > self.threshold
        return pixels_to_mask


class MaskPixelsFromWCS(BaseMask):
    """
    Processor to mask pixels from a file where WCS coordinates of masked pixels are
    given
    """

    base_key = "maskwcs"

    def __init__(
        self,
        mask_pixels_ra: float | list[float] = None,
        mask_pixels_dec: float | list[float] = None,
        mask_file_key: str = FITS_MASK_KEY,
        write_masked_pixels_to_file: bool = False,
        output_dir: str | Path = "mask",
        only_write_mask: bool = False,
    ):
        super().__init__(
            write_masked_pixels_to_file=write_masked_pixels_to_file,
            output_dir=output_dir,
            only_write_mask=only_write_mask,
        )
        self.mask_pixels_ra = mask_pixels_ra
        self.mask_pixels_dec = mask_pixels_dec
        self.mask_file_key = mask_file_key

        if self.mask_pixels_ra is not None:
            self.mask_file_key = None

    def __str__(self) -> str:
        return "Processor to mask pixels using a  list of RA/Dec."

    def get_mask(self, image) -> np.ndarray:
        """
        loads mask if needed, and returns it

        :return: mask
        """
        wcs = WCS(image.get_header())
        if self.mask_file_key is not None:
            mask_file_path = image.get_header()[self.mask_file_key]
            with fits.open(mask_file_path) as mask_image:
                pixels_to_keep = mask_image[0].data  # pylint: disable=no-member
                mask_wcs = WCS(mask_image[0].header)  # pylint: disable=no-member

            masked_pixel_x, masked_pixel_y = np.where(pixels_to_keep == 0.0)

            mask_pixel_coords = mask_wcs.pixel_to_world(masked_pixel_y, masked_pixel_x)
            mask_pixels_ra = mask_pixel_coords.ra.deg
            mask_pixels_dec = mask_pixel_coords.dec.deg

        else:
            mask_pixels_ra = self.mask_pixels_ra
            mask_pixels_dec = self.mask_pixels_dec

        logger.debug(f"Masking {mask_pixels_ra} ras and {mask_pixels_dec} decs")

        mask_pixel_coords = SkyCoord(mask_pixels_ra, mask_pixels_dec, unit="deg")
        mask_pixels_x, mask_pixels_y = wcs.world_to_pixel(mask_pixel_coords)
        mask_pixels_x = mask_pixels_x.astype(int)
        mask_pixels_y = mask_pixels_y.astype(int)

        pixels_to_keep = np.ones(image.get_data().shape, dtype=bool)

        mask_in_image = np.logical_and(
            mask_pixels_x > 0.0, mask_pixels_x < pixels_to_keep.shape[1]
        ) & np.logical_and(mask_pixels_y > 0.0, mask_pixels_y < pixels_to_keep.shape[0])
        mask_pixels_x = mask_pixels_x[mask_in_image]
        mask_pixels_y = mask_pixels_y[mask_in_image]

        pixels_to_keep[mask_pixels_y, mask_pixels_x] = False

        return ~pixels_to_keep


class WriteMaskedCoordsToFile(BaseMask):
    """
    Processor to write masked coordinates to a file
    """

    base_key = "writemaskedcoords"

    def __init__(self, output_dir: str | Path = "mask", only_write_mask: bool = False):
        super().__init__(
            write_masked_pixels_to_file=True,
            output_dir=output_dir,
            only_write_mask=only_write_mask,
        )

    def get_mask(self, image) -> np.ndarray:
        pixels_to_mask = np.zeros(image.get_data().shape, dtype=bool)

        # For some reason, MASK_VALUE == np.nan returns False. Issue/Feature of numpy?
        # This is a workaround
        if np.isnan(MASK_VALUE):
            pixels_to_mask[np.isnan(image.get_data())] = True
        else:
            pixels_to_mask[image.get_data() == MASK_VALUE] = True
        return pixels_to_mask


class MaskDatasecPixels(BaseMask):
    """
    Processor to mask the data section of an image
    """

    base_key = "maskdatasec"

    def get_mask(self, image: Image) -> np.ndarray:
        """
        Function to mask the data section of an image
        """
        header = image.header
        data = image.get_data()
        datasec = header["DATASEC"].replace("[", "").replace("]", "").split(",")
        datasec_xmin = int(datasec[0].split(":")[0])
        datasec_xmax = int(datasec[0].split(":")[1])
        datasec_ymin = int(datasec[1].split(":")[0])
        datasec_ymax = int(datasec[1].split(":")[1])

        mask = np.zeros(data.shape)
        mask[:, :datasec_xmin] = 1.0
        mask[:, datasec_xmax:] = 1.0
        mask[:datasec_ymin, :] = 1.0
        mask[datasec_ymax:, :] = 1.0
        return mask.astype(bool)
