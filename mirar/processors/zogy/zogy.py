"""
Script to run ZOGY.

This is performed by the :class:`mirar.processors.zogy.zogy.ZOGY`
processor. This processor requires several header keys to be present.

In most cases, a processor chain will first require
:class:`mirar.processors.zogy.zogy.ZOGYPrepare`
in order to set these relevant header paths.
"""

import logging
from collections.abc import Callable
from copy import copy
from pathlib import Path

import astropy.table
import astropy.units as u
import numpy as np
from astropy.coordinates import SkyCoord
from astropy.io import fits
from astropy.stats import sigma_clipped_stats
from astropy.table import Table

from mirar.data import Image, ImageBatch
from mirar.data.utils import write_regions_file
from mirar.errors import ProcessorError
from mirar.paths import (
    BASE_NAME_KEY,
    LATEST_WEIGHT_SAVE_KEY,
    MAGLIM_KEY,
    NORM_PSFEX_KEY,
    OBSCLASS_KEY,
    RAW_IMG_KEY,
    REF_IMG_KEY,
    RMS_COUNTS_KEY,
    SCI_IMG_KEY,
    SCOR_IMG_KEY,
    SCOR_MEAN_KEY,
    SCOR_MEDIAN_KEY,
    SCOR_STD_KEY,
    SEXTRACTOR_HEADER_KEY,
    UNC_IMG_KEY,
    core_fields,
    get_output_dir,
)
from mirar.processors.base_processor import BaseImageProcessor, PrerequisiteError
from mirar.processors.zogy.pyzogy import pyzogy
from mirar.utils.ldac_tools import get_table_from_ldac

logger = logging.getLogger(__name__)


class ZOGYError(ProcessorError):
    """Error derived from running ZOGY"""


def default_catalog_purifier(sci_catalog: Table, ref_catalog: Table):
    """
    TODO: This should be in wirc?

    :param sci_catalog:
    :param ref_catalog:
    :return:
    """
    good_sci_sources = (
        (sci_catalog["FLAGS"] == 0)
        & (sci_catalog["SNR_WIN"] > 5)
        & (sci_catalog["FWHM_WORLD"] < 4.0 / 3600)
        & (sci_catalog["FWHM_WORLD"] > 0.5 / 3600)
        & (sci_catalog["SNR_WIN"] < 1000)
    )

    good_ref_sources = (
        (ref_catalog["FLAGS"] == 0)
        & (ref_catalog["SNR_WIN"] > 5)
        & (ref_catalog["FWHM_WORLD"] < 5.0 / 3600)
        & (ref_catalog["FWHM_WORLD"] > 0.5 / 3600)
        & (ref_catalog["SNR_WIN"] < 1000)
    )
    return good_sci_sources, good_ref_sources


class ZOGYPrepare(BaseImageProcessor):
    """
    Class to create a cross-matched source catalogue for ZOGY,
    save the catalogue and then add relevant info to the header.

    After processing by this class, an
     :class:`~mirar.data.image_data.Image`
     can then by processed by an
    :class:`mirar.processors.zogy.zogy.ZOGY` object.
    """

    base_key = "ZOGYPREP"

    def __init__(
        self,
        output_sub_dir: str = "sub",
        sci_zp_header_key: str = "ZP",
        ref_zp_header_key: str | None = None,
        catalog_purifier: Callable[
            [astropy.table.Table, astropy.table.Table],
            [astropy.table.Table, astropy.table.Table],
        ] = default_catalog_purifier,
        write_region_bool: bool = True,
        crossmatch_radius_arcsec: float = 1.0,
    ):
        super().__init__()
        self.output_sub_dir = output_sub_dir
        if ref_zp_header_key is None:
            ref_zp_header_key = sci_zp_header_key
        self.sci_zp_header_key = sci_zp_header_key
        self.ref_zp_header_key = ref_zp_header_key
        self.catalog_purifier = catalog_purifier
        self.write_region_bool = write_region_bool
        self.crossmatch_radius_arcsec = crossmatch_radius_arcsec

    def get_sub_output_dir(self) -> Path:
        """
        Get output directory for this processor

        :return: Path to output directory
        """
        return Path(get_output_dir(self.output_sub_dir, self.night_sub_dir))

    def get_path(self, name: str | Path) -> Path:
        """
        Get output path for a file of name

        :param name: Name of file
        :return: Full output path
        """
        return self.get_sub_output_dir().joinpath(name)

    def get_ast_fluxscale(
        self, ref_catalog_name: str | Path, sci_catalog_name: str | Path
    ) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
        """
        Cross match science and reference data catalogs to get flux scaling factor and
        astometric uncertainties

        :param ref_catalog_name: path of ref catalog
        :param sci_catalog_name: path of sci catalog
        :return: returns astrometric uncertainties in x/y and flux scale
        """
        logger.debug(f"Reference catalog is at {ref_catalog_name}")
        ref_catalog = get_table_from_ldac(ref_catalog_name)
        sci_catalog = get_table_from_ldac(sci_catalog_name)

        logging.info(
            f"Number of total sources SCI: {len(ref_catalog)}, REF: {len(sci_catalog)}"
        )
        good_sci_sources, good_ref_sources = self.catalog_purifier(
            sci_catalog, ref_catalog
        )
        logger.debug(
            f"Number of good sources SCI: {np.sum(good_sci_sources)} "
            f"REF: {np.sum(good_ref_sources)}"
        )

        # Filter catalog for quality
        ref_catalog = ref_catalog[good_ref_sources]
        sci_catalog = sci_catalog[good_sci_sources]

        # Optionally write region file
        if self.write_region_bool:
            ref_reg_name = ref_catalog_name + ".goodsources.reg"
            sci_reg_name = sci_catalog_name + ".goodsources.reg"

            write_regions_file(
                regions_path=ref_reg_name,
                x_coords=ref_catalog["X_IMAGE"],
                y_coords=ref_catalog["Y_IMAGE"],
            )

            write_regions_file(
                regions_path=sci_reg_name,
                x_coords=sci_catalog["X_IMAGE"],
                y_coords=sci_catalog["Y_IMAGE"],
            )

        sci_coords = SkyCoord(
            ra=sci_catalog["ALPHAWIN_J2000"],
            dec=sci_catalog["DELTAWIN_J2000"],
            frame="icrs",
        )
        ref_coords = SkyCoord(
            ra=ref_catalog["ALPHAWIN_J2000"],
            dec=ref_catalog["DELTAWIN_J2000"],
            frame="icrs",
        )

        # Cross match the catalogs
        idx_ref, idx_sci, d2d, _ = sci_coords.search_around_sky(
            ref_coords, self.crossmatch_radius_arcsec * u.arcsec
        )

        if len(d2d) == 0:
            err = (
                "No stars matched between science and reference data catalogs. "
                "Likely there is a huge mismatch in the sensitivites of the two."
            )
            logger.error(err)
            raise ZOGYError(err)

        xpos_sci = sci_catalog["XWIN_IMAGE"]
        ypos_sci = sci_catalog["YWIN_IMAGE"]
        xpos_ref = ref_catalog["XWIN_IMAGE"]
        ypos_ref = ref_catalog["YWIN_IMAGE"]
        sci_flux_auto = sci_catalog["FLUX_AUTO"]
        ref_flux_auto = ref_catalog["FLUX_AUTO"]

        logger.debug(f"Number of cross-matched sources is {len(d2d)}")

        ast_unc_x = np.std(xpos_sci[idx_sci] - xpos_ref[idx_ref])
        ast_unc_y = np.std(ypos_sci[idx_sci] - ypos_ref[idx_ref])

        logger.debug(f"Astrometric uncertainties are X: {ast_unc_x} Y: {ast_unc_y}")

        _, flux_scale_median, flux_scale_std = sigma_clipped_stats(
            sci_flux_auto[idx_sci] / ref_flux_auto[idx_ref]
        )
        flux_scale = flux_scale_median
        logger.debug(
            f"Flux scaling for reference is {flux_scale:.5f} +/- {flux_scale_std:.5f}"
        )
        return ast_unc_x, ast_unc_y, flux_scale

    @staticmethod
    def get_rms_image(image: Image, rms: float) -> Image:
        """Get an RMS image from a regular image

        :param image: An :class:`~mirar.data.image_data.Image`
        :param rms: rms of the image
        :return: An RMS :class:`~mirar.data.image_data.Image`
        """
        gain = image["GAIN"]

        poisson_noise = np.copy(image.get_data()) / gain
        poisson_noise[poisson_noise < 0.0] = 0.0
        rms_image = Image(
            data=np.sqrt(poisson_noise + rms**2), header=image.get_header()
        )
        return rms_image

    def _apply_to_images(self, batch: ImageBatch) -> ImageBatch:
        for image_ind, image in enumerate(batch):
            ref_img_path = image[REF_IMG_KEY]
            sci_img_path = image[BASE_NAME_KEY]

            ref_img = self.open_fits(self.get_path(ref_img_path))

            sci_x_imgsize = int(image["NAXIS1"])
            sci_y_imgsize = int(image["NAXIS2"])
            ref_x_imgsize = int(ref_img["NAXIS1"])
            ref_y_imgsize = int(ref_img["NAXIS2"])

            assert sci_x_imgsize == ref_x_imgsize, (
                "Science and reference images " "have different x dimensions"
            )
            assert sci_y_imgsize == ref_y_imgsize, (
                "Science and reference images " "have different y dimensions"
            )

            trimming_required = (sci_x_imgsize % 2 == 1) or (sci_y_imgsize % 2 == 1)
            if sci_x_imgsize % 2 == 1:
                sci_x_imgsize -= 1
            if sci_y_imgsize % 2 == 1:
                sci_y_imgsize -= 1

            if trimming_required:
                logger.debug(
                    "Trimming science image, original size was "
                    f"{ref_x_imgsize} x {ref_y_imgsize}"
                    f" new size is {sci_x_imgsize} x {sci_y_imgsize}"
                )
                image.set_data(image.get_data()[0:sci_y_imgsize, 0:sci_x_imgsize])
                image["NAXIS1"] = sci_x_imgsize
                image["NAXIS2"] = sci_y_imgsize

                logger.debug("Trimming reference image")
                ref_img.set_data(ref_img.get_data()[0:sci_y_imgsize, 0:sci_x_imgsize])
                ref_img["NAXIS1"] = sci_x_imgsize
                ref_img["NAXIS2"] = sci_y_imgsize
                logger.debug(f"Saving trimmed reference image to {ref_img_path}")
                self.save_fits(ref_img, ref_img_path)

                logger.debug(f"Trimming science weight image")
                weight_path = image[LATEST_WEIGHT_SAVE_KEY]
                with fits.open(
                    self.get_path(weight_path), "update", memmap=False
                ) as weight_img:
                    weight_data = weight_img[0].data  # pylint: disable=no-member
                    weight_img[0].data = weight_data[:sci_y_imgsize, :sci_x_imgsize]
                    weight_img[0].header[
                        "NAXIS1"
                    ] = sci_x_imgsize  # pylint: disable=no-member
                    weight_img[0].header[
                        "NAXIS2"
                    ] = sci_y_imgsize  # pylint: disable=no-member

                logger.debug("Trimming reference weight image")
                ref_weight_path = ref_img[LATEST_WEIGHT_SAVE_KEY]
                with fits.open(
                    self.get_path(ref_weight_path), "update", memmap=False
                ) as weight_img:
                    weight_data = weight_img[0].data  # pylint: disable=no-member
                    weight_img[0].data = weight_data[:sci_y_imgsize, :sci_x_imgsize]
                    weight_img[0].header["NAXIS1"] = sci_x_imgsize
                    weight_img[0].header["NAXIS2"] = sci_y_imgsize

            ref_catalog_path = ref_img[SEXTRACTOR_HEADER_KEY]
            ref_weight_path = ref_img[LATEST_WEIGHT_SAVE_KEY]

            sci_catalog_path = image[SEXTRACTOR_HEADER_KEY]
            sci_weight_path = image[LATEST_WEIGHT_SAVE_KEY]

            ref_weight_data = self.open_fits(self.get_path(ref_weight_path))
            sci_weight_data = self.open_fits(self.get_path(sci_weight_path))

            image_mask = (sci_weight_data.get_data() == 0.0) | (
                ref_weight_data.get_data() == 0.0
            )

            image_data = image.get_data()
            image_data[image_mask] = np.nan
            image.set_data(image_data)
            ref_data = ref_img.get_data()
            ref_data[image_mask] = np.nan
            ref_img.set_data(ref_data)

            # Create S correlation ('scorr') weight image
            scorr_weight_data = sci_weight_data.get_data() * ref_weight_data.get_data()
            scorr_header = sci_weight_data.get_header()
            scorr_weight_img = Image(data=scorr_weight_data, header=scorr_header)
            scorr_weight_path = sci_img_path.replace(".fits", ".scorr.weight.fits")
            scorr_weight_img[OBSCLASS_KEY] = "SCORR"

            for key in core_fields:
                if key not in scorr_weight_img:
                    scorr_weight_img[key] = image[key]

            self.save_fits(
                scorr_weight_img,
                path=self.get_path(scorr_weight_path),
            )

            ast_unc_x, ast_unc_y, flux_scale = self.get_ast_fluxscale(
                ref_catalog_path, sci_catalog_path
            )

            # Scale reference image
            ref_data *= flux_scale
            ref_img.set_data(ref_data)

            ref_unscaled_zp = ref_img[self.ref_zp_header_key]
            ref_img[self.ref_zp_header_key] = float(
                ref_img[self.ref_zp_header_key]
            ) + 2.5 * np.log10(flux_scale)

            logger.debug(
                f"Zeropoints are reference : {ref_unscaled_zp}, "
                f"scaled reference : {ref_img[self.ref_zp_header_key]} and "
                f"science : {image[self.sci_zp_header_key]}"
            )

            sci_rms = 0.5 * (
                np.percentile(image_data[~image_mask], 84.13)
                - np.percentile(image_data[~image_mask], 15.86)
            )
            ref_rms = 0.5 * (
                np.percentile(ref_data[~image_mask], 84.13)
                - np.percentile(ref_data[~image_mask], 15.86)
            )
            logger.debug(
                f"Science RMS is {sci_rms:.2f}. Reference RMS is {ref_rms:.2f}"
            )

            # Calculate uncertainty images
            sci_rms_image = self.get_rms_image(image, sci_rms)
            ref_rms_image = self.get_rms_image(ref_img, ref_rms)

            sci_rms_path = self.get_path(sci_img_path + ".unc.fits")
            ref_rms_path = self.get_path(ref_img_path + ".unc.fits")

            self.save_fits(sci_rms_image, path=sci_rms_path)
            image[UNC_IMG_KEY] = sci_rms_path.as_posix()

            self.save_fits(ref_rms_image, path=ref_rms_path)
            ref_img[UNC_IMG_KEY] = ref_rms_path.as_posix()
            ref_img[RMS_COUNTS_KEY] = ref_rms
            # Save scaled reference image
            ref_scaled_path = ref_img_path.replace(".fits", ".scaled.fits")
            self.save_fits(ref_img, path=self.get_path(ref_scaled_path))

            # Header keywords only required by ZOGY
            image[RMS_COUNTS_KEY] = sci_rms
            image["ASTUNCX"] = ast_unc_x
            image["ASTUNCY"] = ast_unc_y
            image["REFFS"] = flux_scale
            image[REF_IMG_KEY] = str(ref_scaled_path)
            image["SCORMASK"] = str(scorr_weight_path)

            batch[image_ind] = image
        return batch


class ZOGY(ZOGYPrepare):
    """
    :class:`mirar.processors.base_processor.BaseProcessor` class to run
    the ZOGY algorithm using the
    :func:mirar.processors.zogy.pyzogy.pyzogy` function.
    """

    base_key = "ZOGY"
    max_n_cpu = 1

    def __init__(
        self,
        *args,
        output_sub_dir: str = "sub",
        sci_zp_header_key: str = "ZP",
        **kwargs,
    ):
        super().__init__(*args, **kwargs)
        self.output_sub_dir = output_sub_dir
        self.sci_zp_header_key = sci_zp_header_key

    def _apply_to_images(
        self,
        batch: ImageBatch,
    ) -> ImageBatch:
        diff_batch = ImageBatch()
        for image in batch:
            ref_image_path = self.get_path(image[REF_IMG_KEY])
            ref_image = self.open_fits(ref_image_path)

            sci_rms = image[RMS_COUNTS_KEY]
            ref_rms = ref_image[RMS_COUNTS_KEY]
            sci_psf_path = self.get_path(image[NORM_PSFEX_KEY])
            ref_psf_path = self.get_path(ref_image[NORM_PSFEX_KEY])
            sci_rms_path = image[UNC_IMG_KEY]
            ref_rms_path = ref_image[UNC_IMG_KEY]
            ast_unc_x = image["ASTUNCX"]
            ast_unc_y = image["ASTUNCY"]

            # temp_files = [sci_image_path, ref_image_path, sci_rms_path, ref_rms_path]

            logger.debug(f"Ast unc x is {ast_unc_x:.2f} and y is {ast_unc_y:.2f}")
            logger.debug(f"Running zogy on image {image[BASE_NAME_KEY]}")

            # Load the PSFs into memory
            with fits.open(sci_psf_path, memmap=False) as img_psf_f:
                new_psf = img_psf_f[0].data  # pylint: disable=no-member

            with fits.open(ref_psf_path, memmap=False) as ref_psf_f:
                ref_psf = ref_psf_f[0].data  # pylint: disable=no-member

            # Load the sigma images into memory
            with fits.open(sci_rms_path, memmap=False) as img_sigma_f:
                new_sigma = img_sigma_f[0].data  # pylint: disable=no-member

            with fits.open(ref_rms_path, memmap=False) as ref_sigma_f:
                ref_sigma = ref_sigma_f[0].data  # pylint: disable=no-member

            diff_data, diff_psf_data, scorr_data = pyzogy(
                new_data=image.get_data(),
                ref_data=ref_image.get_data(),
                new_psf=new_psf,
                ref_psf=ref_psf,
                new_sigma=new_sigma,
                ref_sigma=ref_sigma,
                new_avg_unc=sci_rms,
                ref_avg_unc=ref_rms,
                dx=ast_unc_x,
                dy=ast_unc_y,
            )

            sci_image_path = self.get_path(image[BASE_NAME_KEY])
            diff_image_path = Path(sci_image_path).with_suffix(".diff.fits")
            diff_psf_path = diff_image_path.with_suffix(".psf")

            scorr_image_path = Path(sci_image_path).with_suffix(".scorr.fits")

            scorr_mean, scorr_median, scorr_std = sigma_clipped_stats(
                scorr_data, mask_value=np.nan
            )

            logger.debug(
                f"Scorr mean, median, STD is {scorr_mean}, {scorr_median}, {scorr_std}"
            )

            sci_rms_image = self.open_fits(self.get_path(sci_rms_path))
            ref_rms_image = self.open_fits(self.get_path(ref_rms_path))

            diff_rms_data = np.sqrt(
                sci_rms_image.get_data() ** 2 + ref_rms_image.get_data() ** 2
            )
            _, diff_rms_median, _ = sigma_clipped_stats(
                diff_rms_data[~np.isnan(diff_rms_data)], mask_value=np.nan
            )
            diff_rms_path = diff_image_path.with_suffix(".unc.fits")

            diff = Image(data=diff_data, header=copy(image.get_header()))

            diff[NORM_PSFEX_KEY] = diff_psf_path.as_posix()
            diff[SCOR_IMG_KEY] = scorr_image_path.as_posix()
            diff[UNC_IMG_KEY] = diff_rms_path.as_posix()
            noise = np.sqrt(
                np.nansum(np.square(diff_psf_data) * np.square(diff_rms_median))
            ) / np.nansum(np.square(diff_psf_data))

            diff[MAGLIM_KEY] = -2.5 * np.log10(noise * 5) + float(
                diff[self.sci_zp_header_key]
            )
            key_map = {
                SCOR_MEAN_KEY: scorr_mean,
                SCOR_MEDIAN_KEY: scorr_median,
                SCOR_STD_KEY: scorr_std,
            }
            for key, value in key_map.items():
                if np.isnan(value):
                    value = None
                diff[key] = value

            self.save_fits(image=diff, path=self.get_path(diff_image_path))

            psf_header = fits.Header({"SIMPLE": True})
            psf_header[BASE_NAME_KEY] = diff_psf_path.name
            psf_header[RAW_IMG_KEY] = diff_psf_path.as_posix()

            for key in core_fields:
                if key not in psf_header:
                    psf_header[key] = diff[key]

            self.save_fits(
                image=Image(diff_psf_data, psf_header),
                path=self.get_path(diff_psf_path),
            )

            scorr = Image(scorr_data, header=image.header.copy())
            scorr[LATEST_WEIGHT_SAVE_KEY] = diff["SCORMASK"]

            self.save_fits(image=scorr, path=self.get_path(scorr_image_path))

            self.save_fits(
                image=Image(data=diff_rms_data, header=copy(image.get_header())),
                path=self.get_path(diff_rms_path),
            )

            diff[BASE_NAME_KEY] = diff_image_path.name
            diff[NORM_PSFEX_KEY] = diff_psf_path.as_posix()
            diff[UNC_IMG_KEY] = diff_rms_path.as_posix()
            diff[SCI_IMG_KEY] = image[BASE_NAME_KEY]
            diff_batch.append(diff)
        return diff_batch

    def check_prerequisites(
        self,
    ):
        check = np.sum([isinstance(x, ZOGYPrepare) for x in self.preceding_steps])
        if check != 1:
            raise PrerequisiteError("ZOGYPrepare must be run before ZOGY")
