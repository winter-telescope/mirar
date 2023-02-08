"""
Module with base class for reference image generation
"""
import logging
import os
from pathlib import Path

import numpy as np
from astropy.io import fits

from winterdrp.data import Image
from winterdrp.io import save_hdu_as_fits
from winterdrp.paths import BASE_NAME_KEY, COADD_KEY, PROC_HISTORY_KEY

logger = logging.getLogger(__name__)


class BaseReferenceGenerator:
    """
    Base Ref Image Generator
    """

    @property
    def abbreviation(self):
        """
        Abbreviation for image naming
        """
        raise NotImplementedError()

    def __init__(self, filter_name: str):
        self.filter_name = filter_name

    def get_reference(self, image: Image) -> fits.PrimaryHDU:
        """
        Get loaded ref image for image

        :param image: image
        :return: ref image
        """
        raise NotImplementedError()

    def write_reference(self, image: Image, output_dir: str) -> Path:
        """
        Write reference image to file

        :param image: Image
        :param output_dir: directory to write to
        :return: path of reference image
        """

        base_name = os.path.basename(image[BASE_NAME_KEY])
        logger.debug(f"Base name is {base_name}")

        ref_hdu = self.get_reference(image)

        output_path = Path(
            self.get_output_path(output_dir, base_name).replace(".fits", "")
            + "_ref.fits"
        )

        # This is because Swarp requires the COADDS keyword. I am setting it to
        # zero manually
        if COADD_KEY not in ref_hdu.header.keys():
            logger.debug("Setting COADDS to 1")
            ref_hdu.header[COADD_KEY] = 1
        if PROC_HISTORY_KEY not in ref_hdu.header.keys():
            logger.debug("Setting CALSTEPS to blank")
            ref_hdu.header[PROC_HISTORY_KEY] = ""

        # Remove if needed
        output_path.unlink(missing_ok=True)

        logger.info(f"Saving reference image to {output_path}")
        ref_hdu.header[BASE_NAME_KEY] = os.path.basename(output_path)
        ref_hdu.data[ref_hdu.data == 0] = np.nan  # pylint: disable=no-member
        save_hdu_as_fits(ref_hdu, output_path)

        return output_path

    @staticmethod
    def get_output_path(output_dir: str, base_name: str) -> str:
        """
        Get output path

        :param output_dir: output dir
        :param base_name: Name
        :return: output path
        """
        return os.path.join(output_dir, base_name)
