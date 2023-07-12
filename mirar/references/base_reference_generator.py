"""
Module with base class for reference image generation
"""
import logging
import os
from pathlib import Path
from typing import Type

import numpy as np
from astropy.io import fits
from astropy.time import Time

from mirar.data import Image, ImageBatch
from mirar.errors import ProcessorError
from mirar.io import (
    MissingCoreFieldError,
    check_image_has_core_fields,
    save_hdu_as_fits,
)
from mirar.paths import (
    BASE_NAME_KEY,
    COADD_KEY,
    LATEST_SAVE_KEY,
    LATEST_WEIGHT_SAVE_KEY,
    OBSCLASS_KEY,
    PROC_FAIL_KEY,
    PROC_HISTORY_KEY,
    RAW_IMG_KEY,
    TARGET_KEY,
)
from mirar.processors.sqldatabase.base_model import BaseDB
from mirar.processors.sqldatabase.database_exporter import DatabaseImageExporter

logger = logging.getLogger(__name__)


class ReferenceGenerationError(ProcessorError):
    """Error for reference generation"""


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

    def __init__(
        self,
        filter_name: str,
        write_to_db: bool = False,
        db_table: Type[BaseDB] = None,
        duplicate_protocol: str = "replace",
        q3c_bool: bool = True,
    ):
        self.filter_name = filter_name
        self.write_to_db = write_to_db
        self.write_db_table = db_table
        self.duplicate_protocol = duplicate_protocol
        self.q3c_bool = q3c_bool

        if np.logical_and(self.write_to_db, self.write_db_table is None):
            err = "You have set write_to_db=True but not provided a write_db_table."
            raise ReferenceError(err)

    def get_reference(self, image: Image) -> (fits.PrimaryHDU, fits.PrimaryHDU):
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

        ref_hdu, ref_weight_hdu = self.get_reference(image)

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

        if "FIELDID" in image.header.keys():
            ref_hdu.header["FIELDID"] = image.header["FIELDID"]
        if "SUBDETID" in image.header.keys():
            ref_hdu.header["SUBDETID"] = image.header["SUBDETID"]

        # Remove if needed
        output_path.unlink(missing_ok=True)

        logger.debug(f"Saving reference image to {output_path}")
        ref_hdu.header[BASE_NAME_KEY] = os.path.basename(output_path)
        ref_hdu.header[LATEST_SAVE_KEY] = output_path.as_posix()
        ref_hdu.header[RAW_IMG_KEY] = image[BASE_NAME_KEY]
        ref_hdu.header[OBSCLASS_KEY] = "REF"
        ref_hdu.header[TARGET_KEY] = image[TARGET_KEY]
        ref_hdu.header[PROC_FAIL_KEY] = False

        if ("MJD-OBS" in ref_hdu.header.keys()) & (
            "DATE-OBS" not in ref_hdu.header.keys()
        ):
            ref_hdu.header["DATE-OBS"] = Time(
                ref_hdu.header["MJD-OBS"], format="mjd"
            ).isot

        ref_hdu.data[ref_hdu.data == 0] = np.nan  # pylint: disable=no-member

        if ref_weight_hdu is not None:
            output_weight_path = Path(
                self.get_output_path(output_dir, base_name).replace(".fits", "")
                + "_ref_weight.fits"
            )
            output_weight_path.unlink(missing_ok=True)
            ref_weight_hdu.header[BASE_NAME_KEY] = os.path.basename(output_weight_path)
            ref_weight_hdu.header[LATEST_SAVE_KEY] = output_weight_path.as_posix()

            save_hdu_as_fits(ref_weight_hdu, output_weight_path)
            ref_hdu.header[LATEST_WEIGHT_SAVE_KEY] = output_weight_path.as_posix()

        # Reference images should have all the core fields
        try:
            check_image_has_core_fields(Image(header=ref_hdu.header, data=ref_hdu.data))
        except MissingCoreFieldError as err:
            raise ReferenceGenerationError from err

        save_hdu_as_fits(ref_hdu, output_path)

        if self.write_to_db:
            dbexporter = DatabaseImageExporter(
                db_table=self.write_db_table,
                duplicate_protocol=self.duplicate_protocol,
                q3c_bool=self.q3c_bool,
            )
            ref_image = Image(header=ref_hdu.header, data=ref_hdu.data)
            ref_image_batch = ImageBatch([ref_image])
            _ = dbexporter.apply(ref_image_batch)

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
