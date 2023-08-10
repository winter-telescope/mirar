"""
Module with base class for reference image generation
"""
import logging
import os
from pathlib import Path
from typing import Callable, Type

import numpy as np
from astropy.io import fits
from astropy.time import Time

from mirar.data import Image, ImageBatch
from mirar.data.utils import get_corners_ra_dec_from_header, get_image_center_wcs_coords
from mirar.database.base_model import BaseDB
from mirar.errors import ProcessorError
from mirar.io import (
    MissingCoreFieldError,
    check_image_has_core_fields,
    open_fits,
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
    get_output_dir,
)
from mirar.processors.astromatic.sextractor.sextractor import Sextractor
from mirar.processors.astromatic.swarp import Swarp
from mirar.processors.database.database_inserter import DatabaseImageInserter
from mirar.processors.photcal import PhotCalibrator

logger = logging.getLogger(__name__)


class ReferenceGenerationError(ProcessorError):
    """Error for reference generation"""


class BaseReferenceGenerator:
    """
    Base Reference Image Generator.
    Subclasses should implement get_reference, which
    returns a reference image HDU, and a reference weight image HDU.
    Classes using this generator should call get_reference_image to get the reference
    image as an Image object.
    Users can optionally choose to save the reference images to a file and/or database.
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
        write_image: bool = True,
        write_image_sub_dir: str = "references",
        write_to_db: bool = False,
        db_table: Type[BaseDB] = None,
        duplicate_protocol: str = "replace",
        q3c_bool: bool = True,
    ):
        """
        filter_name: filter name
        write_image: write reference image to file?
        write_image_sub_dir: directory to write reference image to, required if
        write_image=True, defaults to 'references'
        write_to_db: write reference image to database?
        db_table: database table to write to, required if write_to_db=True
        duplicate_protocol: protocol for handling duplicate entries in database
        q3c_bool: use q3c for database queries?
        """
        self.filter_name = filter_name
        self.write_to_db = write_to_db
        self.write_db_table = db_table
        self.duplicate_protocol = duplicate_protocol
        self.q3c_bool = q3c_bool
        self.write_image = write_image
        self.write_image_dir = write_image_sub_dir

        if np.logical_and(self.write_to_db, self.write_db_table is None):
            err = "You have set write_to_db=True but not provided a write_db_table."
            raise ReferenceError(err)

    def _get_reference(self, image: Image) -> (fits.PrimaryHDU, fits.PrimaryHDU):
        """
        Get loaded ref image for image

        :param image: image
        :return: reference image HDU, reference weight image HDU
        """
        raise NotImplementedError()

    def get_reference_image(self, image: Image) -> Image:
        """
        Get reference image corresponding to an image

        :param image: Image
        :return: reference image
        """

        base_name = os.path.basename(image[BASE_NAME_KEY])
        logger.debug(f"Base name is {base_name}")

        ref_hdu, ref_weight_hdu = self._get_reference(image)

        output_dir = get_output_dir(self.write_image_dir)
        output_dir.mkdir(exist_ok=True, parents=True)
        output_path = Path(output_dir).joinpath(base_name.replace(".fits", "_ref.fits"))

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

        logger.debug(f"Saving reference image to {output_path}")
        ref_hdu.header[BASE_NAME_KEY] = os.path.basename(output_path)
        ref_hdu.header[RAW_IMG_KEY] = os.path.basename(output_path)
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

        ref_image = Image(header=ref_hdu.header, data=ref_hdu.data)

        # Reference images should have all the core fields
        try:
            check_image_has_core_fields(ref_image)
        except MissingCoreFieldError as err:
            raise ReferenceGenerationError from err

        if ref_weight_hdu is not None:
            output_weight_path = Path(
                str(self.get_output_path(output_dir, base_name)).replace(".fits", "")
                + "_ref.weight.fits"
            )
            output_weight_path.unlink(missing_ok=True)
            ref_weight_hdu.header[BASE_NAME_KEY] = os.path.basename(output_weight_path)
            ref_weight_hdu.header[OBSCLASS_KEY] = "WEIGHT"
            ref_weight_hdu.header[TARGET_KEY] = image[TARGET_KEY]
            ref_weight_hdu.header[PROC_FAIL_KEY] = False
            if self.write_image:
                ref_weight_hdu.header[LATEST_SAVE_KEY] = output_weight_path.as_posix()
                save_hdu_as_fits(ref_weight_hdu, output_weight_path)
                ref_hdu.header[LATEST_WEIGHT_SAVE_KEY] = output_weight_path.as_posix()

        if self.write_image:
            ref_hdu.header[LATEST_SAVE_KEY] = output_path.as_posix()
            save_hdu_as_fits(ref_hdu, output_path)

        if self.write_to_db:
            dbexporter = DatabaseImageInserter(
                db_table=self.write_db_table,
                duplicate_protocol=self.duplicate_protocol,
            )
            ref_image_batch = ImageBatch([ref_image])
            _ = dbexporter.apply(ref_image_batch)

        return ref_image

    @staticmethod
    def get_output_path(output_dir: str | Path, base_name: str) -> Path:
        """
        Get output path

        :param output_dir: output dir
        :param base_name: Name
        :return: output path
        """
        return Path(output_dir).joinpath(base_name)


class BaseStackReferenceGenerator(BaseReferenceGenerator):
    """
    Base processor where multiple images need to be stacked together
    """

    def __init__(
        self,
        filter_name: str,
        image_resampler_generator: Callable[..., Swarp],
        write_stacked_image: bool = True,
        write_stack_sub_dir: str = "references/stack",
        write_stack_to_db: bool = False,
        stacks_db_table: Type[BaseDB] = None,
        duplicate_protocol: str = "replace",
        q3c_bool: bool = False,
        stack_image_annotator: Callable[[Image, Image], Image] = None,
        references_base_subdir_name: str = "references",
        photcal_stack: bool = False,
        sextractor_generator: Callable[..., Sextractor] = None,
        phot_calibrator_generator: Callable[..., PhotCalibrator] = None,
    ):
        """
        Initialise
        filter_name: Filter name
        image_resampler_generator: FUnction returning a Swarp object for resampling
        write_stacked_image: Write stacked image to file?
        write_stack_sub_dir: Subdirectory to write stacked image to. Defaults to
        "references/stack"
        write_stack_to_db: Write stacked image to db
        stacks_db_table: Table to write stack to, required if write_stack_to_db is True
        duplicate_protocol: Duplicate protocol
        q3c_bool: Use q3c
        stack_image_annotator: Function to optionally annotate the stack image with
        pipeline-specific information after it has been stacked
        references_base_subdir_name: Base name for the references subdirectory, with
        components etc.
        photcal_stack: Do you want to photometrically calibrate the stack? If not, make
        sure zeropoints are being calculated in the get_components method.
        sextractor_generator: Function returning a Sextractor object for source
        extraction, required if photcal_stack is True
        phot_calibrator_generator: Function returning a PhotCalibrator object for
        photometric calibration, required if photcal_stack is True
        """
        super().__init__(
            filter_name=filter_name,
            write_to_db=write_stack_to_db,
            db_table=stacks_db_table,
            duplicate_protocol=duplicate_protocol,
            q3c_bool=q3c_bool,
            write_image=write_stacked_image,
            write_image_sub_dir=write_stack_sub_dir,
        )
        self.image_resampler_generator = image_resampler_generator
        self.stack_image_annotator = stack_image_annotator
        self.references_base_subdir_name = references_base_subdir_name
        self.photcal_stack = photcal_stack
        self.sextractor_generator = sextractor_generator
        self.phot_calibrator_generator = phot_calibrator_generator
        if self.photcal_stack:
            if (self.sextractor_generator is None) | (
                self.phot_calibrator_generator is None
            ):
                raise ValueError(
                    "Must provide sextractor and phot_calibrator "
                    "generators if photcal_stack is True"
                )

    def get_component_images(self, image: Image) -> ImageBatch:
        """
        Get component reference images that will be stacked together.
        :param image: image
        :return: ImageBatch of component reference images to be used for stacking
        """
        raise NotImplementedError

    def _get_reference(self, image: Image) -> (fits.PrimaryHDU, fits.PrimaryHDU):
        """
        Get loaded ref image for image

        :param image: image
        :return: ref image HDU, ref weight HDU
        """
        component_image_batch = self.get_component_images(image)
        resampler = self.image_resampler_generator()

        resampler.set_night(night_sub_dir=self.references_base_subdir_name)
        resampled_batch = resampler.apply(component_image_batch)
        stacked_image = resampled_batch[0]

        reference_weight_path = resampled_batch[0].header[LATEST_WEIGHT_SAVE_KEY]
        reference_weight_data, reference_weight_header = open_fits(
            reference_weight_path
        )

        if self.photcal_stack:
            ref_sextractor = self.sextractor_generator(resampled_batch[0])
            ref_sextractor.set_night(night_sub_dir=self.references_base_subdir_name)
            resampled_batch = ref_sextractor.apply(resampled_batch)

            phot_calibrator = self.phot_calibrator_generator(resampled_batch[0])
            phot_calibrator.set_night(night_sub_dir=self.references_base_subdir_name)
            phot_calibrator.set_preceding_steps([ref_sextractor])
            photcaled_batch = phot_calibrator.apply(resampled_batch)

            stacked_image = photcaled_batch[0]

        ra_cent, dec_cent = get_image_center_wcs_coords(image=stacked_image, origin=1)

        (
            (ra0_0, dec0_0),
            (ra0_1, dec0_1),
            (ra1_0, dec1_0),
            (ra1_1, dec1_1),
        ) = get_corners_ra_dec_from_header(stacked_image.header)

        stacked_image.header["RA_CENT"] = ra_cent
        stacked_image.header["DEC_CENT"] = dec_cent
        stacked_image.header["RA0_0"] = ra0_0
        stacked_image.header["RA0_1"] = ra0_1
        stacked_image.header["RA1_0"] = ra1_0
        stacked_image.header["RA1_1"] = ra1_1
        stacked_image.header["DEC0_0"] = dec0_0
        stacked_image.header["DEC0_1"] = dec0_1
        stacked_image.header["DEC1_0"] = dec1_0
        stacked_image.header["DEC1_1"] = dec1_1

        if self.stack_image_annotator is not None:
            stacked_image = self.stack_image_annotator(stacked_image, image)

        stack_hdu = fits.PrimaryHDU(
            data=stacked_image.get_data(), header=stacked_image.header
        )
        stack_weight_hdu = fits.PrimaryHDU(
            data=reference_weight_data, header=reference_weight_header
        )
        return stack_hdu, stack_weight_hdu
