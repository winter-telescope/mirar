"""
Script to run ZOGY.

This is performed by the :class:`mirar.processors.zogy.zogy.ZOGY`
processor. This processor requires several header keys to be present.

In most cases, a processor chain will first require
:class:`mirar.processors.zogy.zogy.ZOGYPrepare`
in order to set these relevant header paths.
"""
import logging
import os
from collections.abc import Callable
from pathlib import Path

import astropy.table
import astropy.units as u
import numpy as np
from astropy.coordinates import SkyCoord
from astropy.io import fits
from astropy.stats import sigma_clipped_stats
from astropy.table import Table

from mirar.data import Image, ImageBatch
from mirar.errors import ProcessorError
from mirar.paths import (
    BASE_NAME_KEY,
    LATEST_WEIGHT_SAVE_KEY,
    NORM_PSFEX_KEY,
    OBSCLASS_KEY,
    RAW_IMG_KEY,
    REF_IMG_KEY,
    REF_PSF_KEY,
    UNC_IMG_KEY,
    core_fields,
    get_output_dir,
)
from mirar.processors.base_processor import BaseImageProcessor
from mirar.processors.candidates.utils.regions_writer import write_regions_file
from mirar.processors.zogy.pyzogy import pyzogy
from mirar.utils.ldac_tools import get_table_from_ldac

logger = logging.getLogger(__name__)


class ZOGYError(ProcessorError):
    """Error derived from running ZOGY"""


def default_wirc_catalog_purifier(sci_catalog: Table, ref_catalog: Table):
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


def default_summer_catalog_purifier(sci_catalog: Table, ref_catalog: Table):
    # Need to do this because the summer data is typically much
    # shallower than the PS1 data, and only the brightest
    # sources in PS1 xmatch to it.
    good_sci_sources = (
        (sci_catalog["FLAGS"] == 0)
        & (sci_catalog["SNR_WIN"] > 5)
        & (sci_catalog["FWHM_WORLD"] < 4.0 / 3600)
        & (sci_catalog["FWHM_WORLD"] > 0.5 / 3600)
        & (sci_catalog["SNR_WIN"] < 1000)
    )

    good_ref_sources = (
        (ref_catalog["SNR_WIN"] > 5)
        & (ref_catalog["FWHM_WORLD"] < 5.0 / 3600)
        & (ref_catalog["FWHM_WORLD"] > 0.5 / 3600)
    )

    return good_sci_sources, good_ref_sources


def default_sedmv2_catalog_purifier(sci_catalog, ref_catalog):
    good_sci_sources = (
        (sci_catalog["FLAGS"] == 0)
        & (sci_catalog["SNR_WIN"] > 5)
        & (sci_catalog["FWHM_WORLD"] < 4.0 / 3600)
        & (sci_catalog["FWHM_WORLD"] > 0.5 / 3600)
        & (sci_catalog["SNR_WIN"] < 1000)
    )

    good_ref_sources = (
        (ref_catalog["SNR_WIN"] > 5)
        & (ref_catalog["FWHM_WORLD"] < 5.0 / 3600)
        & (ref_catalog["FWHM_WORLD"] > 0.5 / 3600)
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
        ] = default_wirc_catalog_purifier,
        write_region_bool: bool = False,
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
        logger.info(f"Reference catalog is at {ref_catalog_name}")
        ref_catalog = get_table_from_ldac(ref_catalog_name)
        sci_catalog = get_table_from_ldac(sci_catalog_name)

        logging.info(
            f"Number of total sources SCI: {len(ref_catalog)}, REF: {len(sci_catalog)}"
        )
        good_sci_sources, good_ref_sources = self.catalog_purifier(
            sci_catalog, ref_catalog
        )
        logger.info(
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

        logger.info(f"Number of cross-matched sources is {len(d2d)}")

        ast_unc_x = np.std(xpos_sci[idx_sci] - xpos_ref[idx_ref])
        ast_unc_y = np.std(ypos_sci[idx_sci] - ypos_ref[idx_ref])

        logger.info(f"Astrometric uncertainties are X: {ast_unc_x} Y: {ast_unc_y}")

        _, flux_scale_median, flux_scale_std = sigma_clipped_stats(
            sci_flux_auto[idx_sci] / ref_flux_auto[idx_ref]
        )
        flux_scale = flux_scale_median
        logger.info(
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
        poisson_noise[poisson_noise < 0] = 0
        rms_image = Image(
            data=np.sqrt(poisson_noise + rms**2), header=image.get_header()
        )
        return rms_image

    def _apply_to_images(self, batch: ImageBatch) -> ImageBatch:
        for image in batch:
            ref_img_path = image[REF_IMG_KEY]
            sci_img_path = image[BASE_NAME_KEY]

            ref_img = self.open_fits(self.get_path(ref_img_path))

            ref_catalog_path = ref_img["SRCCAT"]  # convert to key
            ref_weight_path = ref_img[LATEST_WEIGHT_SAVE_KEY]
            # change in Swarp too, also fix (get_mask/write_mask)?

            sci_catalog_path = image["SRCCAT"]  # convert to key
            sci_weight_path = image[LATEST_WEIGHT_SAVE_KEY]

            ref_weight_data = self.open_fits(self.get_path(ref_weight_path))
            sci_weight_data = self.open_fits(self.get_path(sci_weight_path))

            image_mask = (sci_weight_data.get_data() == 0.0) | (
                ref_weight_data.get_data() == 0.0
            )

            image_data = image.get_data()
            image_data[image_mask] = 0.0
            image.set_data(image_data)
            ref_data = ref_img.get_data()
            ref_data[image_mask] = 0.0
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

            ref_scaled_path = ref_img_path + ".scaled"
            self.save_fits(ref_img, path=self.get_path(ref_scaled_path))

            logger.debug(
                f"Zeropoints are reference : {ref_unscaled_zp}, "
                f"scaled reference : {ref_img[self.ref_zp_header_key]} and "
                f"science : {image[self.ref_zp_header_key]}"
            )

            # Scale is 1 by construction for science image
            sci_scaled_path = self.get_path(sci_img_path + ".scaled")
            self.save_fits(image, path=sci_scaled_path)

            sci_rms = 0.5 * (
                np.percentile(image_data[image_data != 0.0], 84.13)
                - np.percentile(image_data[image_data != 0.0], 15.86)
            )
            ref_rms = 0.5 * (
                np.percentile(ref_data[ref_data != 0.0], 84.13)
                - np.percentile(ref_data[ref_data != 0.0], 15.86)
            )
            logger.debug(
                f"Science RMS is {sci_rms:.2f}. Reference RMS is {ref_rms:.2f}"
            )

            # Calculate uncertainty images
            sci_rms_image = self.get_rms_image(image, sci_rms)
            ref_rms_image = self.get_rms_image(ref_img, ref_rms)

            sci_rms_path = sci_img_path + ".unc"
            ref_rms_path = ref_img_path + ".unc"

            self.save_fits(
                sci_rms_image, path=os.path.join(self.get_path(sci_rms_path))
            )

            self.save_fits(
                ref_rms_image, path=os.path.join(self.get_path(ref_rms_path))
            )

            image["SCIRMS"] = sci_rms
            image["REFRMS"] = ref_rms
            image["ASTUNCX"] = ast_unc_x
            image["ASTUNCY"] = ast_unc_y
            image["REFFS"] = flux_scale
            image["SCUNCPTH"] = str(sci_rms_path)
            image["RFUNCPTH"] = str(ref_rms_path)
            image["SCISCL"] = str(sci_scaled_path)
            image["REFSCL"] = str(ref_scaled_path)
            image["SCORMASK"] = str(scorr_weight_path)

        return batch


class ZOGY(ZOGYPrepare):
    """
    :class:`mirar.processors.base_processor.BaseProcessor` class to run
    the ZOGY algorithm using the
    :func:mirar.processors.zogy.pyzogy.pyzogy` function.
    """

    base_key = "ZOGY"

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
            sci_image_path = self.get_path(image["SCISCL"])
            ref_image_path = self.get_path(image["REFSCL"])
            sci_rms = image["SCIRMS"]
            ref_rms = image["REFRMS"]
            sci_psf_path = self.get_path(image[NORM_PSFEX_KEY])
            ref_psf_path = self.get_path(image[REF_PSF_KEY])
            sci_rms_path = self.get_path(image["SCUNCPTH"])
            ref_rms_path = self.get_path(image["RFUNCPTH"])
            ast_unc_x = image["ASTUNCX"]
            ast_unc_y = image["ASTUNCY"]

            # temp_files = [sci_image_path, ref_image_path, sci_rms_path, ref_rms_path]

            diff_data, diff_psf_data, scorr_data = pyzogy(
                new_image_path=sci_image_path,
                ref_image_path=ref_image_path,
                new_psf_path=sci_psf_path,
                ref_psf_path=ref_psf_path,
                new_sigma_path=sci_rms_path,
                ref_sigma_path=ref_rms_path,
                new_avg_unc=sci_rms,
                ref_avg_unc=ref_rms,
                dx=ast_unc_x,
                dy=ast_unc_y,
            )

            diff = Image(data=diff_data, header=image.get_header())

            diff_image_path = Path(sci_image_path).with_suffix(".diff.fits")
            diff_psf_path = diff_image_path.with_suffix(".psf")

            scorr_image_path = Path(sci_image_path).with_suffix(".scorr.fits")

            scorr_mean, scorr_median, scorr_std = sigma_clipped_stats(scorr_data)

            logger.info(
                f"Scorr mean, median, STD is {scorr_mean}, {scorr_median}, {scorr_std}"
            )

            sci_rms_image = self.open_fits(self.get_path(sci_rms_path))
            ref_rms_image = self.open_fits(self.get_path(ref_rms_path))

            diff_rms_data = np.sqrt(
                sci_rms_image.get_data() ** 2 + ref_rms_image.get_data() ** 2
            )
            _, diff_rms_median, _ = sigma_clipped_stats(diff_rms_data)
            diff_rms_path = diff_image_path.with_suffix(".unc")

            image["DIFFIMG"] = diff_image_path.as_posix()
            image["DIFFPSF"] = diff_psf_path.as_posix()
            image["DIFFSCR"] = scorr_image_path.as_posix()
            image["DIFFUNC"] = diff_rms_path.as_posix()
            noise = np.sqrt(
                np.nansum(np.square(diff_psf_data) * np.square(diff_rms_median))
            ) / np.nansum(np.square(diff_psf_data))
            image["DIFFMLIM"] = -2.5 * np.log10(noise * 5) + float(
                image[self.sci_zp_header_key]
            )
            image["SCORMEAN"] = scorr_mean
            image["SCORMED"] = scorr_median
            image["SCORSTD"] = scorr_std

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
            scorr[LATEST_WEIGHT_SAVE_KEY] = image["SCORMASK"]

            self.save_fits(image=scorr, path=self.get_path(scorr_image_path))

            self.save_fits(
                image=Image(data=diff_rms_data, header=image.get_header()),
                path=self.get_path(diff_rms_path),
            )

            diff[BASE_NAME_KEY] = diff_image_path.name
            diff[NORM_PSFEX_KEY] = diff_psf_path.as_posix()
            diff[UNC_IMG_KEY] = diff_rms_path.as_posix()

            diff_batch.append(diff)
        return diff_batch
