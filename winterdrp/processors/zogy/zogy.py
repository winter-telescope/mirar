import logging

from winterdrp.processors.base_processor import BaseProcessor
from winterdrp.paths import get_output_dir, latest_mask_save_key
from astropy.io import fits
import numpy as np
from winterdrp.utils.ldac_tools import get_table_from_ldac
from astropy.coordinates import SkyCoord
from astropy.stats import sigma_clipped_stats
import astropy.units as u
from winterdrp.io import open_fits
from winterdrp.processors.zogy.py_zogy import py_zogy
import os

logger = logging.getLogger(__name__)

class ZOGY(BaseProcessor):

    def __init__(self,
                 output_sub_dir: str = "sub",
                 *args,
                 **kwargs):
        super(ZOGY, self).__init__(*args, **kwargs)
        self.output_sub_dir = output_sub_dir

    def _apply_to_images(
            self,
            images: list[np.ndarray],
            headers: list[fits.Header],
    ) -> tuple[list[np.ndarray], list[fits.Header]]:

        for ind, header in enumerate(headers):
            sci_image_path = header["SCISCL"]
            ref_image_path = header["REFSCL"]
            sci_rms = header["SCIRMS"]
            ref_rms = header["REFRMS"]
            sci_psf = header["SCIPSF"]
            ref_psf = header["REFPSF"]
            sci_rms_image = header["SCUNCPTH"]
            ref_rms_image = header["RFUNCPTH"]
            ast_unc_x = header["ASTUNCX"]
            ast_unc_y = header["ASTUNCY"]

            D, P_D, S_corr = py_zogy(sci_image_path, ref_image_path, sci_psf, ref_psf, sci_rms_image,
                                             ref_rms_image, sci_rms, ref_rms, dx=ast_unc_x, dy=ast_unc_y)

            diff_image_path = sci_image_path + '.diff'
            self.save_fits(data=D,
                           header=header,
                           path=os.path.join(self.output_sub_dir, diff_image_path))

            diff_psf_path = sci_image_path + '.psf'
            self.save_fits(data=P_D,
                           header=None,
                           path=os.path.join(self.output_sub_dir, diff_psf_path))

            scorr_image_path = sci_image_path + '.scorr'
            self.save_fits(data=S_corr,
                           header=header,
                           path=os.path.join(self.output_sub_dir, scorr_image_path))

            header["DIFFIMG"] = diff_image_path
            header["DIFFPSF"] = diff_psf_path
            header["DIFFSCR"] = scorr_image_path

        return images, headers


class ZOGYPrepare(BaseProcessor):

    def __init__(self,
                 output_sub_dir: str = "sub",
                 *args,
                 **kwargs):
        super(ZOGYPrepare, self).__init__(*args, **kwargs)
        self.output_sub_dir = output_sub_dir

    def get_sub_output_dir(self):
        return get_output_dir(self.output_sub_dir, self.night_sub_dir)

    @staticmethod
    def get_ast_fluxscale(ref_catalog_name: str,
                          sci_catalog_name: str) -> tuple[float, float, float]:
        # Cross match science and reference image catalogs to get flux scaling factor and astometric uncertainties
        sci_catalog = get_table_from_ldac(ref_catalog_name)
        ref_catalog = get_table_from_ldac(sci_catalog_name)

        good_sci_sources = (sci_catalog['FLAGS'] == 0) & (sci_catalog['SNR_WIN'] > 5) & (
                sci_catalog['FWHM_WORLD'] < 4. / 3600) & (sci_catalog['FWHM_WORLD'] > 0.5 / 3600) & (
                                   sci_catalog['SNR_WIN'] < 1000)
        good_ref_sources = (ref_catalog['FLAGS'] == 0) & (ref_catalog['SNR_WIN'] > 5) & (
                ref_catalog['FWHM_WORLD'] < 5. / 3600) & (ref_catalog['FWHM_WORLD'] > 0.5 / 3600) & (
                                   ref_catalog['SNR_WIN'] < 1000)

        print(f'Number of good sources SCI: {np.sum(good_sci_sources)} REF: {np.sum(good_ref_sources)}')
        ref_catalog = ref_catalog[good_ref_sources]
        sci_catalog = sci_catalog[good_sci_sources]

        sci_coords = SkyCoord(ra=sci_catalog['ALPHAWIN_J2000'], dec=sci_catalog['DELTAWIN_J2000'], frame='icrs')
        ref_coords = SkyCoord(ra=ref_catalog['ALPHAWIN_J2000'], dec=ref_catalog['DELTAWIN_J2000'], frame='icrs')

        # Cross match the catalogs
        idx_ref, idx_sci, d2d, d3d = sci_coords.search_around_sky(ref_coords, 1.0 * u.arcsec)

        xpos_sci = sci_catalog['XWIN_IMAGE']
        ypos_sci = sci_catalog['YWIN_IMAGE']
        xpos_ref = ref_catalog['XWIN_IMAGE']
        ypos_ref = ref_catalog['YWIN_IMAGE']
        sci_flux_auto = sci_catalog['FLUX_AUTO']
        ref_flux_auto = ref_catalog['FLUX_AUTO']

        logger.info('Number of cross-matched sources is %d' % (len(d2d)))

        ast_unc_x = np.std(xpos_sci[idx_sci] - xpos_ref[idx_ref])
        ast_unc_y = np.std(ypos_sci[idx_sci] - ypos_ref[idx_ref])

        # print(np.median(xpos_sci[idx_sci] - xpos_ref[idx_ref]))
        with open('goodmatches.reg', 'w') as f:
            f.write('wcs\n')
            for i in range(len(xpos_ref[idx_ref])):
                f.write('point (%s,%s) #point=cross\n' % (xpos_sci[idx_sci][i], ypos_sci[idx_sci][i]))

        logger.info(f'Astrometric uncertainties are X: {ast_unc_x} Y: {ast_unc_y}')
        #print('Mean of astrometric uncertainties is X:%.2f Y:%.2f' % (
        #    np.mean(xpos_sci[idx_sci] - xpos_ref[idx_ref]), np.mean(xpos_sci[idx_sci] - xpos_ref[idx_ref])))
        flux_scale_mean, flux_scale_median, flux_scale_std = sigma_clipped_stats(
            sci_flux_auto[idx_sci] / ref_flux_auto[idx_ref])
        flux_scale = flux_scale_median
        print('Flux scaling for reference is %.5f +/- %.5f' % (flux_scale, flux_scale_std))
        return ast_unc_x, ast_unc_y, flux_scale

    @staticmethod
    def get_rms_image(image: np.ndarray, header: fits.Header, rms: float) -> np.ndarray:
        gain = header["GAIN"]
        poisson_noise = np.copy(image) / gain
        poisson_noise[poisson_noise < 0] = 0
        rms_image = np.sqrt(poisson_noise + rms ** 2)
        return rms_image

    def _apply_to_images(
            self,
            images: list[np.ndarray],
            headers: list[fits.Header],
    ) -> tuple[list[np.ndarray], list[fits.Header]]:
        for ind, header in enumerate(headers):
            ref_img_path = header["REFIMG"]
            sci_img_path = header["BASENAME"]
            sci_data = images[ind]

            ref_data, ref_header = self.open_fits(os.path.join(self.output_sub_dir, ref_img_path))
            ref_catalog_path = ref_header["SRCCAT"]
            sci_catalog_path = header["SRCCAT"]

            ast_unc_x, ast_unc_y, flux_scale = self.get_ast_fluxscale(ref_catalog_path, sci_catalog_path)

            ref_data = ref_data * flux_scale

            ref_scaled_path = ref_img_path + '.scaled'
            self.save_fits(data=ref_data,
                           header=ref_header,
                           path=os.path.join(self.output_sub_dir,ref_scaled_path)
                           )

            sci_scaled_path = sci_img_path + '.scaled'
            self.save_fits(data=sci_data,
                           header=header,
                           path=os.path.join(self.output_sub_dir,sci_scaled_path)
                           )

            sci_rms = 0.5 * (
                    np.percentile(sci_data[sci_data != 0], 84.13) - np.percentile(sci_data[sci_data != 0], 15.86))
            ref_rms = 0.5 * (
                    np.percentile(ref_data[ref_data != 0], 84.13) - np.percentile(ref_data[ref_data != 0], 15.86))
            print('Science RMS is %.2f. Reference RMS is %.2f' % (sci_rms, ref_rms))

            sci_weight_path = header[latest_mask_save_key]
            ref_weight_path = header[latest_mask_save_key]
            sci_weight_data, _ = open_fits(sci_weight_path)
            ref_weight_data, _ = open_fits(ref_weight_path)

            image_mask = (sci_weight_data == 0) | (ref_weight_data == 0)
            sci_data[image_mask] = 0
            ref_data[image_mask] = 0

            # Calculate uncertainty images
            sci_rms_image = self.get_rms_image(sci_data, header, sci_rms)
            ref_rms_image = self.get_rms_image(ref_data, ref_header, ref_rms)

            sci_rms_path = sci_img_path + '.unc'
            ref_rms_path = ref_img_path + '.unc'

            self.save_fits(data=sci_rms_image,
                           header=header,
                           path=os.path.join(self.output_sub_dir,sci_rms_path)
                           )

            self.save_fits(data=ref_rms_image,
                           header=ref_header,
                           path=os.path.join(self.output_sub_dir,ref_rms_path)
                           )

            header['SCIRMS'] = sci_rms
            header['REFRMS'] = ref_rms
            header['ASTUNCX'] = ast_unc_x
            header['ASTUNCY'] = ast_unc_y
            header['REFFS'] = flux_scale
            header['SCUNCPTH'] = sci_rms_path
            header["RFUNCPTH"] = ref_rms_path
            header["SCISCL"] = sci_scaled_path
            header["REFSCL"] = ref_scaled_path

        return images, headers
