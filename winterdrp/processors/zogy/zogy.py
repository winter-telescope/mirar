import logging

import astropy.table
from collections.abc import Callable
from winterdrp.processors.base_processor import BaseImageProcessor
from winterdrp.paths import get_output_dir, latest_mask_save_key
from astropy.io import fits
import numpy as np
from winterdrp.utils.ldac_tools import get_table_from_ldac
from astropy.coordinates import SkyCoord
from astropy.stats import sigma_clipped_stats
import astropy.units as u
from winterdrp.io import open_fits
from winterdrp.processors.zogy.pyzogy import pyzogy
from winterdrp.paths import norm_psfex_header_key, base_name_key, ref_img_key, ref_psf_key
import os
from winterdrp.processors.candidates.utils.regions_writer import write_regions_file
from winterdrp.errors import ProcessorError
from winterdrp.data import ImageBatch, Image
from pathlib import Path

logger = logging.getLogger(__name__)


class ZOGYError(ProcessorError):
    pass


def default_wirc_catalog_purifier(sci_catalog, ref_catalog):
    good_sci_sources = (sci_catalog['FLAGS'] == 0) & (sci_catalog['SNR_WIN'] > 5) & (
            sci_catalog['FWHM_WORLD'] < 4. / 3600) & (sci_catalog['FWHM_WORLD'] > 0.5 / 3600) & (
                               sci_catalog['SNR_WIN'] < 1000)

    good_ref_sources = (ref_catalog['FLAGS'] == 0) & (ref_catalog['SNR_WIN'] > 5) & (
             ref_catalog['FWHM_WORLD'] < 5. / 3600) & (ref_catalog['FWHM_WORLD'] > 0.5 / 3600) & (
                                ref_catalog['SNR_WIN'] < 1000)
    return good_sci_sources, good_ref_sources


def default_summer_catalog_purifier(sci_catalog, ref_catalog):
    # Need to do this because the summer data is typically much shallower than the PS1 data, and only the brightest
    # sources in PS1 xmatch to it.
    good_sci_sources = (sci_catalog['FLAGS'] == 0) & (sci_catalog['SNR_WIN'] > 5) & (
            sci_catalog['FWHM_WORLD'] < 4. / 3600) & (sci_catalog['FWHM_WORLD'] > 0.5 / 3600) & (
                               sci_catalog['SNR_WIN'] < 1000)

    good_ref_sources = (ref_catalog['SNR_WIN'] > 5) & (
            ref_catalog['FWHM_WORLD'] < 5. / 3600) & (ref_catalog['FWHM_WORLD'] > 0.5 / 3600)

    return good_sci_sources, good_ref_sources


class ZOGYPrepare(BaseImageProcessor):
    base_key = "ZOGYPREP"

    def __init__(self,
                 output_sub_dir: str = "sub",
                 sci_zp_header_key: str = "ZP",
                 catalog_purifier: Callable[[astropy.table.Table, astropy.table.Table],[astropy.table.Table, astropy.table.Table]] = default_wirc_catalog_purifier,
                 *args,
                 **kwargs):
        super(ZOGYPrepare, self).__init__(*args, **kwargs)
        self.output_sub_dir = output_sub_dir
        self.sci_zp_header_key = sci_zp_header_key
        self.catalog_purifier = catalog_purifier

    def get_sub_output_dir(self) -> Path:
        return Path(get_output_dir(self.output_sub_dir, self.night_sub_dir))

    def get_path(self, name: str) -> Path:
        return self.get_sub_output_dir().joinpath(name)

    def get_ast_fluxscale(
            self,
            ref_catalog_name: str,
            sci_catalog_name: str
    ) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
        # Cross match science and reference data catalogs to get flux scaling factor and astometric uncertainties
        logger.info(f'Reference catalog is at {ref_catalog_name}')
        ref_catalog = get_table_from_ldac(ref_catalog_name)
        sci_catalog = get_table_from_ldac(sci_catalog_name)

        good_sci_sources, good_ref_sources = self.catalog_purifier(sci_catalog, ref_catalog)
        logger.info(f'Number of good sources SCI: {np.sum(good_sci_sources)} REF: {np.sum(good_ref_sources)}')
        ref_catalog = ref_catalog[good_ref_sources]
        sci_catalog = sci_catalog[good_sci_sources]

        ref_reg_name = ref_catalog_name + '.goodsources.reg'
        sci_reg_name = sci_catalog_name + '.goodsources.reg'
        write_regions_file(regions_path=ref_reg_name,
                           x_coords=ref_catalog['X_IMAGE'],
                           y_coords=ref_catalog['Y_IMAGE'])

        write_regions_file(regions_path=sci_reg_name,
                           x_coords=sci_catalog['X_IMAGE'],
                           y_coords=sci_catalog['Y_IMAGE'])

        sci_coords = SkyCoord(ra=sci_catalog['ALPHAWIN_J2000'], dec=sci_catalog['DELTAWIN_J2000'], frame='icrs')
        ref_coords = SkyCoord(ra=ref_catalog['ALPHAWIN_J2000'], dec=ref_catalog['DELTAWIN_J2000'], frame='icrs')

        # Cross match the catalogs
        idx_ref, idx_sci, d2d, d3d = sci_coords.search_around_sky(ref_coords, 1.0 * u.arcsec)

        if len(d2d) == 0:
            err = 'No stars matched between science and reference data catalogs. Likely there is a huge mismatch in ' \
                  'the sensitivites of the two'
            logger.error(err)
            raise ZOGYError(err)

        xpos_sci = sci_catalog['XWIN_IMAGE']
        ypos_sci = sci_catalog['YWIN_IMAGE']
        xpos_ref = ref_catalog['XWIN_IMAGE']
        ypos_ref = ref_catalog['YWIN_IMAGE']
        sci_flux_auto = sci_catalog['FLUX_AUTO']
        ref_flux_auto = ref_catalog['FLUX_AUTO']

        logger.info('Number of cross-matched sources is %d' % (len(d2d)))

        ast_unc_x = np.std(xpos_sci[idx_sci] - xpos_ref[idx_ref])
        ast_unc_y = np.std(ypos_sci[idx_sci] - ypos_ref[idx_ref])

        logger.info(f'Astrometric uncertainties are X: {ast_unc_x} Y: {ast_unc_y}')

        flux_scale_mean, flux_scale_median, flux_scale_std = sigma_clipped_stats(
            sci_flux_auto[idx_sci] / ref_flux_auto[idx_ref])
        flux_scale = flux_scale_median
        logger.info('Flux scaling for reference is %.5f +/- %.5f' % (flux_scale, flux_scale_std))
        return ast_unc_x, ast_unc_y, flux_scale

    @staticmethod
    def get_rms_image(image: Image, rms: float) -> Image:
        gain = image["GAIN"]
        poisson_noise = np.copy(image.get_data()) / gain
        poisson_noise[poisson_noise < 0] = 0
        rms_image = Image(data=np.sqrt(poisson_noise + rms ** 2), header=image.get_header())
        return rms_image

    def _apply_to_images(
            self,
            batch: ImageBatch
    ) -> ImageBatch:
        for ind, image in enumerate(batch):
            ref_img_path = image[ref_img_key]
            sci_img_path = image[base_name_key]

            ref_img = self.open_fits(self.get_path(ref_img_path))

            ref_catalog_path = ref_img["SRCCAT"]
            ref_mask_path = ref_img["MASKPATH"]

            sci_catalog_path = image["SRCCAT"]
            sci_mask_path = image["MASKPATH"]

            ref_weight_img = self.open_fits(self.get_path(ref_mask_path))
            sci_weight_img = self.open_fits(self.get_path(sci_mask_path))

            image_mask = (sci_weight_img == 0) | (ref_weight_img == 0)
            image.data[image_mask] = 0
            ref_img.data[image_mask] = 0

            scorr_mask_data = sci_weight_img.data * ref_weight_img.data
            scorr_header = sci_weight_img.get_header()

            scorr_weight_img = Image(data=scorr_mask_data, header=scorr_header)

            scorr_mask_path = sci_img_path.replace('.fits', '.scorr.weight.fits')
            self.save_fits(scorr_weight_img,
                           path=self.get_path(scorr_mask_path))

            ast_unc_x, ast_unc_y, flux_scale = self.get_ast_fluxscale(ref_catalog_path, sci_catalog_path)

            ref_img.data *= flux_scale
            ref_unscaled_zp = ref_img['ZP']
            ref_img['ZP'] = float(ref_img['ZP']) + 2.5 * np.log10(flux_scale)

            ref_scaled_path = ref_img_path + '.scaled'
            self.save_fits(ref_img, path=self.get_path(ref_scaled_path))

            logger.info(f"Zeropoints are reference : {ref_unscaled_zp}, "
                        f"scaled reference : {ref_img['ZP']} and "
                        f"science : {image[self.sci_zp_header_key]}")
            # assert(False)
            sci_scaled_path = self.get_path(sci_img_path + '.scaled')

            self.save_fits(
                image,
                path=sci_scaled_path
            )

            sci_rms = 0.5 * (
                    np.percentile(image.data[image.data != 0], 84.13)
                    - np.percentile(image.data[image.data != 0], 15.86)
            )
            ref_rms = 0.5 * (
                    np.percentile(ref_img.data[ref_img.data != 0], 84.13)
                    - np.percentile(ref_img.data[ref_img.data != 0], 15.86)
            )
            logger.info(f'Science RMS is {sci_rms:.2f}. Reference RMS is {ref_rms:.2f}')

            sci_weight_path = self.get_path(image[latest_mask_save_key])
            ref_weight_path = self.get_path(image[latest_mask_save_key])
            sci_weight_img, _ = open_fits(sci_weight_path)
            ref_weight_img, _ = open_fits(ref_weight_path)

            image_mask = (sci_weight_img == 0) | (ref_weight_img == 0)
            image.data[image_mask] = 0
            ref_img.data[image_mask] = 0

            # Calculate uncertainty images
            sci_rms_image = self.get_rms_image(image, sci_rms)
            ref_rms_image = self.get_rms_image(ref_img, ref_rms)

            sci_rms_path = sci_img_path + '.unc'
            ref_rms_path = ref_img_path + '.unc'

            self.save_fits(sci_rms_image,
                           path=os.path.join(self.get_path(sci_rms_path))
                           )

            self.save_fits(ref_rms_image,
                           path=os.path.join(self.get_path(ref_rms_path))
                           )

            image['SCIRMS'] = sci_rms
            image['REFRMS'] = ref_rms
            image['ASTUNCX'] = ast_unc_x
            image['ASTUNCY'] = ast_unc_y
            image['REFFS'] = flux_scale
            image['SCUNCPTH'] = str(sci_rms_path)
            image["RFUNCPTH"] = str(ref_rms_path)
            image["SCISCL"] = str(sci_scaled_path)
            image["REFSCL"] = str(ref_scaled_path)
            image['SCORMASK'] = str(scorr_mask_path)

        return batch


class ZOGY(ZOGYPrepare):
    base_key = "ZOGY"

    def __init__(self,
                 output_sub_dir: str = "sub",
                 sci_zp_header_key: str = "ZP",
                 *args,
                 **kwargs):
        super(ZOGY, self).__init__(*args, **kwargs)
        self.output_sub_dir = output_sub_dir
        self.sci_zp_header_key = sci_zp_header_key

    def _apply_to_images(
            self,
            batch: ImageBatch,
    ) -> ImageBatch:
        diff_batch = ImageBatch()
        for ind, image in enumerate(batch):
            sci_image_path = self.get_path(image["SCISCL"])
            ref_image_path = self.get_path(image["REFSCL"])
            sci_rms = image["SCIRMS"]
            ref_rms = image["REFRMS"]
            sci_psf_path = self.get_path(image[norm_psfex_header_key])
            ref_psf_path = self.get_path(image[ref_psf_key])
            sci_rms_path = self.get_path(image["SCUNCPTH"])
            ref_rms_path = self.get_path(image["RFUNCPTH"])
            ast_unc_x = image["ASTUNCX"]
            ast_unc_y = image["ASTUNCY"]

            temp_files = [sci_image_path, ref_image_path, sci_rms_path, ref_rms_path]

            print(image["SCISCL"])

            diff_data, diff_psf_data, s_corr_data = pyzogy(
                new_image_path=sci_image_path,
                ref_image_path=ref_image_path,
                new_psf_path=sci_psf_path,
                ref_psf_path=ref_psf_path,
                new_sigma_path=sci_rms_path,
                ref_sigma_path=ref_rms_path,
                new_avg_unc=sci_rms,
                ref_avg_unc=ref_rms,
                dx=ast_unc_x,
                dy=ast_unc_y
            )

            diff = Image(
                data=diff_data,
                header=image.get_header()
            )

            diff_image_path = sci_image_path.replace('.fits', '') + '.diff.fits'
            diff_psf_path = diff_image_path + '.psf'
            scorr_image_path = sci_image_path.replace('.fits', '') + '.scorr.fits'

            scorr_mean, scorr_median, scorr_std = sigma_clipped_stats(s_corr)
            logger.info(f"Scorr mean, median, STD is {scorr_mean}, {scorr_median}, {scorr_std}")
            sci_rms_data, _ = self.open_fits(self.get_path(sci_rms_path))
            ref_rms_data, _ = self.open_fits(self.get_path(ref_rms_path))
            diff_rms_image = np.sqrt(sci_rms_data ** 2 + ref_rms_data ** 2)
            diff_rms_mean, diff_rms_median, diff_rms_std = sigma_clipped_stats(diff_rms_image)
            diff_rms_path = diff_image_path + '.unc'

            image["DIFFIMG"] = diff_image_path
            image["DIFFPSF"] = diff_psf_path
            image["DIFFSCR"] = scorr_image_path
            image["DIFFUNC"] = diff_rms_path
            noise = np.sqrt(np.nansum(np.square(diff_psf) * np.square(diff_rms_median))) / np.nansum(
                np.square(diff_psf))
            image["DIFFMLIM"] = -2.5 * np.log10(noise * 5) + float(image[self.sci_zp_header_key])
            image["SCORMEAN"] = scorr_mean
            image["SCORMED"] = scorr_median
            image["SCORSTD"] = scorr_std

            self.save_fits(image=diff,
                           path=self.get_path(diff_image_path))

            self.save_fits(image=diff_psf,
                           path=self.get_path(diff_psf_path))

            s_corr[latest_mask_save_key] = image['SCORMASK']
            self.save_fits(image=s_corr,
                           path=self.get_path(scorr_image_path))

            self.save_fits(image=diff_rms_image,
                           path=self.get_path(diff_rms_path))
            # logger.info(f"{diff_rms_std}, {diff_rms_median}, {diff_std} DIFFMAGLIM {header['DIFFMLIM']}")
            # # for temp_file in temp_files:

            # for temp_file in temp_files:

            #     os.remove(temp_file)
            #     logger.info(f"Deleted temporary file {temp_file}")
            diff_batch.append(diff)
        return diff_batch
