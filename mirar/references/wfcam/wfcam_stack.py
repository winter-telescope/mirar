"""
Module for creating stacked reference images from WFCAM (UKIRT/VISTA) surveys.
This module can be used to query THe WFAU service to download NIR images.
Only single extension images overlapping with the specified coordinates are downloaded.
These images are then stacked together.

This module works in the following way -
1. Images are queried using the type of query specified by the user. Query needs to be
inherited from BaseWFCAMQuery.
2. The queried images are filtered to remove images with wildly different zeropoints
and poor seeings. The remaining images are stacked together using SWarp. The median
zeropoint is assigned to the stacked image.
3. The stacked image is optionally saved to the user-specified path and inserted into
a user-specified stack-table in the database by the parent class.
"""

import logging
from collections.abc import Callable
from typing import Type

import numpy as np

from mirar.data import Image, ImageBatch
from mirar.database.base_model import BaseDB
from mirar.paths import BASE_NAME_KEY, ZP_KEY, ZP_STD_KEY
from mirar.processors.astromatic.sextractor.sextractor import Sextractor
from mirar.processors.astromatic.swarp import Swarp
from mirar.processors.base_processor import ImageHandler
from mirar.processors.photcal.photcalibrator import PhotCalibrator
from mirar.references.base_reference_generator import (
    BaseStackReferenceGenerator,
    ReferenceGenerationError,
)
from mirar.references.wfcam.wfcam_query import BaseWFCAMQuery

logger = logging.getLogger(__name__)


def default_filter_wfau_images(image_batch: ImageBatch) -> ImageBatch:
    """
    Function to filter WFAU images based on zeropoints and seeing.
    The images are filtered to remove images with wildly different zeropoints (more than
    0.4 mag from the median zeropoint)
    and seeing <= 0 or > 3.5 arcsec.
    Args:
        :param image_batch: ImageBatch

    Returns:
        :return: Filtered image batch
    """
    image_array = np.array([x for x in image_batch])

    mag_zps = np.array(
        [
            x["MAGZPT"]
            + 2.5 * np.log10(x["EXPTIME"])
            - x["EXTINCT"] * ((x["AMSTART"] + x["AMEND"]) / 2)
            for x in image_batch
        ]
    )
    # magerr_zps = np.array([x["MAGZRR"] for x in ukirt_images])
    median_mag_zp = np.median(mag_zps)
    seeings = np.array([x["SEEING"] for x in image_batch])
    zpmask = np.abs(mag_zps - median_mag_zp) < 0.4
    seeingmask = (seeings < 3.5 / 0.4) & (seeings > 0)

    image_array = image_array[zpmask & seeingmask]

    return ImageBatch(image_array.tolist())


class WFCAMStackedRef(BaseStackReferenceGenerator, ImageHandler):
    """
    Class to query UKIRT images from the WFAU archive and stack them together
    """

    abbreviation = "wfcam_stack_online_ref"

    def __init__(
        self,
        filter_name: str,
        image_resampler_generator: Callable[..., Swarp],
        wfcam_query: BaseWFCAMQuery,
        write_stacked_image: bool = True,
        write_stack_sub_dir: str = "references/ref_stacks",
        write_stack_to_db: bool = False,
        stacks_db_table: Type[BaseDB] = None,
        component_image_sub_dir: str = None,
        references_base_subdir_name: str = "references",
        stack_image_annotator: Callable[[Image], Image] = None,
        photcal_stack: bool = False,
        sextractor_generator: Callable[..., Sextractor] = None,
        phot_calibrator_generator: Callable[..., PhotCalibrator] = None,
        filter_images: Callable[[ImageBatch], ImageBatch] = default_filter_wfau_images,
    ):
        """
        Args:
            :param wfcam_query: Query to be used to query WFCAM images
            :component_image_sub_dir: Subdirectory to save the component images
            :param filter_images: Function to filter the queried images
        """
        BaseStackReferenceGenerator.__init__(
            self,
            filter_name=filter_name,
            write_stack_to_db=write_stack_to_db,
            stacks_db_table=stacks_db_table,
            references_base_subdir_name=references_base_subdir_name,
            image_resampler_generator=image_resampler_generator,
            stack_image_annotator=stack_image_annotator,
            photcal_stack=photcal_stack,
            sextractor_generator=sextractor_generator,
            phot_calibrator_generator=phot_calibrator_generator,
            write_stacked_image=write_stacked_image,
            write_stack_sub_dir=write_stack_sub_dir,
        )
        self.wfcam_query = wfcam_query

        self.component_image_dir = component_image_sub_dir
        self.filter_images = filter_images

    def get_component_images(self, image: Image) -> ImageBatch:
        """
        Function to get the component images for a given image
        Args:
            :param image: Image for which the component images are to be found
            :return: ImageBatch of component images
        """
        wfau_images = self.wfcam_query.run_query(image)

        if len(wfau_images) > 0:
            wfau_images = self.filter_images(wfau_images)

        if len(wfau_images) == 0:
            raise ReferenceGenerationError(
                f"No good WFAU images found for {image[BASE_NAME_KEY]}"
            )

        # change BASENAME to gel well with parallel processing
        for ind, ref_img in enumerate(wfau_images):
            new_basename = (
                f"{ref_img[BASE_NAME_KEY].strip('.fits')}" f"_{image[BASE_NAME_KEY]}"
            )
            ref_img[BASE_NAME_KEY] = new_basename

        # Get the scaling factors
        mag_zps = np.array(
            [
                x["MAGZPT"]
                + 2.5 * np.log10(x["EXPTIME"])
                - x["EXTINCT"] * ((x["AMSTART"] + x["AMEND"]) / 2)
                for x in wfau_images
            ]
        )

        median_mag_zp = np.median(mag_zps)
        scaling_factors = 10 ** (0.4 * (median_mag_zp - mag_zps))

        for ind, image_to_resamp in enumerate(wfau_images):
            image_to_resamp["FLXSCALE"] = scaling_factors[ind]
            image_to_resamp[ZP_KEY] = median_mag_zp
            image_to_resamp[ZP_STD_KEY] = np.std(mag_zps)

        wfau_image_batch = ImageBatch(list(wfau_images))

        return wfau_image_batch
