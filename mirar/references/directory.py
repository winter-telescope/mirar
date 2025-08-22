"""
Module for getting reference images in a directory
"""

import logging
from pathlib import Path

from astropy.io import fits

from mirar.data import Image
from mirar.io import open_fits
from mirar.paths import LATEST_WEIGHT_SAVE_KEY, get_output_dir
from mirar.references.base_reference_generator import BaseReferenceGenerator

logger = logging.getLogger(__name__)


class RefFromDirectory(BaseReferenceGenerator):
    """
    Get locally saved ref with
    """

    abbreviation = "refdirectory"

    def __init__(
        self,
        file_base_name: str,
        **kwargs,
    ):
        super().__init__(**kwargs)
        self.file_base_name = file_base_name

    def get_ref_dir(self) -> Path:
        """
        Get directory where reference images are stored

        :return: Path to directory
        """
        input_image_dir = get_output_dir(self.write_image_dir) / "named"
        input_image_dir.mkdir(parents=True, exist_ok=True)
        return input_image_dir

    def _get_reference(self, image: Image) -> (fits.PrimaryHDU, fits.PrimaryHDU):
        """
        Get reference image from directory

        :param image: Image object
        :return: Tuple of (reference HDU, reference weight HDU or None)
        """

        input_image_dir = self.get_ref_dir()

        ref_path = input_image_dir / f"{self.file_base_name}.fits"

        if not ref_path.exists():
            raise FileNotFoundError(f"Reference image {ref_path} not found")

        ref_weight_hdu = None
        data, header = open_fits(ref_path)
        ref_hdu = fits.PrimaryHDU(data, header).copy()
        if LATEST_WEIGHT_SAVE_KEY in header:  # pylint: disable=no-member
            weight_path = Path(
                header[LATEST_WEIGHT_SAVE_KEY]  # pylint: disable=no-member
            )
            if weight_path.exists():
                with fits.open(weight_path) as wght_hdul:
                    ref_weight_hdu = wght_hdul[0].copy()

        return ref_hdu, ref_weight_hdu
