import astropy.io.fits
import numpy as np
import os
import logging
from winterdrp.processors.base_processor import BaseProcessor, ProcessorWithCache
from winterdrp.paths import get_output_dir, copy_temp_file, get_temp_path, get_untemp_path
from collections.abc import Callable
from winterdrp.paths import cal_output_dir
from winterdrp.catalog.base_catalog import BaseCatalog
from winterdrp.processors.astromatic.sextractor.sextractor import Sextractor, sextractor_header_key
from winterdrp.utils.ldac_tools import get_table_from_ldac
from astropy.coordinates import SkyCoord
import astropy.units as u
from astropy.stats import sigma_clip, sigma_clipped_stats

logger = logging.getLogger(__name__)


class PhotCalibrator(BaseProcessor):
    base_key = 'photcalibrator'

    def __init__(self,
                 ref_catalog_generator: Callable[[astropy.io.fits.Header], BaseCatalog],
                 temp_output_sub_dir: str = "phot",
                 redo: bool = True,
                 x_lower_limit: float = 100,
                 x_upper_limit: float = 2800,
                 y_lower_limit: float = 100,
                 y_upper_limit: float = 2800,
                 fwhm_threshold_arcsec: float = 4.0,
                 num_matches_threshold: int = 5,
                 *args,
                 **kwargs):
        super(PhotCalibrator, self).__init__(*args, **kwargs)
        self.redo = redo
        self.ref_catalog_generator = ref_catalog_generator
        self.temp_output_sub_dir = temp_output_sub_dir
        self.x_lower_limit = x_lower_limit
        self.x_upper_limit = x_upper_limit
        self.y_lower_limit = y_lower_limit
        self.y_upper_limit = y_upper_limit
        self.fwhm_threshold_arcsec = fwhm_threshold_arcsec
        self.num_matches_threshold = num_matches_threshold

    def get_phot_output_dir(self):
        return get_output_dir(self.temp_output_sub_dir, self.night_sub_dir)

    def calculate_zeropoint(
            self,
            ref_cat_path,
            img_cat_path) -> list[dict]:
        ref_cat = get_table_from_ldac(ref_cat_path)
        img_cat = get_table_from_ldac(img_cat_path)

        if len(ref_cat) == 0:
            logger.info('No sources found in reference catalog')
            return [{'Error': -1}]

        ref_coords = SkyCoord(ra=ref_cat['RA'],dec=ref_cat['DEC'],unit=(u.deg,u.deg))

        clean_mask = (img_cat['FLAGS'] == 0) & (img_cat['FWHM_WORLD'] < self.fwhm_threshold_arcsec / 3600) \
                     & (img_cat['X_IMAGE'] > self.x_lower_limit) & \
                     (img_cat['X_IMAGE'] < self.x_upper_limit) & \
                     (img_cat['Y_IMAGE'] > self.y_lower_limit) & \
                     (img_cat['Y_IMAGE'] < self.y_upper_limit)

        clean_img_cat = img_cat[clean_mask]
        logger.debug(f'Found {len(clean_img_cat)} clean sources in image.')

        clean_img_coords = SkyCoord(ra=clean_img_cat['ALPHAWIN_J2000'], dec=clean_img_cat['DELTAWIN_J2000'],
                                    unit=(u.deg, u.deg))

        if 0 == len(clean_img_coords):
            logger.info('No clean sources found in image')
            return [{'Error': -2}]

        idx, d2d, d3d = ref_coords.match_to_catalog_sky(clean_img_coords)
        match_mask = d2d < 1.0 * u.arcsec
        matched_ref_cat = ref_cat[match_mask]
        matched_img_cat = clean_img_cat[idx[match_mask]]
        logger.info(f'Cross-matched {len(matched_img_cat)} sources from catalog to the image.')

        if len(matched_img_cat)< self.num_matches_threshold:
            logger.info(f'Not enough cross-matched sources found to calculate a reliable zeropoint.')
            return [{'Error': -3}]

        apertures = np.array([2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0])  # aperture diameters
        zeropoints = []
        for i in range(len(apertures)):
            offsets = np.ma.array(matched_ref_cat['magnitude'] - matched_img_cat['MAG_APER'][:, i])
            cl_offset = sigma_clip(offsets)
            num_stars = np.sum(np.invert(cl_offset.mask))
            # print(np.median(cl_offset))
            zp_mean, zp_med, zp_std = sigma_clipped_stats(offsets)
            if np.isnan(zp_mean):
                zp_mean = -99
            if np.isnan(zp_med):
                zp_med = -99
            if np.isnan(zp_std):
                zp_std = -99
            zero_dict = {'diameter': apertures[i], 'zp_mean': zp_mean, 'zp_median': zp_med, 'zp_std': zp_std,
                        'nstars': num_stars, 'mag_cat': matched_ref_cat['magnitude'][np.invert(cl_offset.mask)],
                        'mag_apers': matched_img_cat['MAG_APER'][:, i][np.invert(cl_offset.mask)]}
            zeropoints.append(zero_dict)

        offsets = np.ma.array(matched_ref_cat['magnitude'] - matched_img_cat['MAG_AUTO'])
        cl_offset = sigma_clip(offsets, sigma=3)
        num_stars = np.sum(np.invert(cl_offset.mask))
        zp_mean, zp_med, zp_std = sigma_clipped_stats(offsets, sigma=3)
        zero_auto_mag_cat = matched_ref_cat['magnitude'][np.invert(cl_offset.mask)]
        zero_auto_mag_img = matched_img_cat['MAG_AUTO'][np.invert(cl_offset.mask)]
        zeropoints.append({'diameter':'AUTO','zp_mean':zp_mean,'zp_median':zp_med,'zp_std':zp_std,
                           'nstars': num_stars, 'mag_cat': zero_auto_mag_cat,
                           'mag_apers': zero_auto_mag_img
                           })

        return zeropoints


    def _apply_to_images(
            self,
            images: list[np.ndarray],
            headers: list[astropy.io.fits.Header],
    ) -> tuple[list[np.ndarray], list[astropy.io.fits.Header]]:

        phot_output_dir = self.get_phot_output_dir()

        try:
            os.makedirs(phot_output_dir)
        except OSError:
            pass

        for header in headers:
            ref_catalog = self.ref_catalog_generator(header)
            ref_cat_path = ref_catalog.write_catalog(header,output_dir=phot_output_dir)
            temp_cat_path = copy_temp_file(
                output_dir=phot_output_dir,
                file_path=header[sextractor_header_key]
            )

            zp_dicts = self.calculate_zeropoint(ref_cat_path,temp_cat_path)

            if 'Error' in zp_dicts[0].keys():
                logger.info(f'Failed to run photometric calibration for ')
                header['PHOTCAL'] = 'FAILED'
                continue

            for zpvals in zp_dicts:
                header['ZP_%s' % (zpvals['diameter'])] = zpvals['zp_mean']
                header['ZP_%s_std' % (zpvals['diameter'])] = zpvals['zp_std']
                header['ZP_%s_nstars' % (zpvals['diameter'])] = zpvals['nstars']

            #header.add_history('Calibrated to SDSS')
            header['PHOTCAL'] = 'SUCCESS'

        return images, headers


    def check_prerequisites(
            self,
    ):
        check = np.sum([isinstance(x, Sextractor) for x in self.preceding_steps])
        if check < 1:
            err = f"{self.__module__} requires {Sextractor} as a prerequisite. " \
                  f"However, the following steps were found: {self.preceding_steps}."
            logger.error(err)
            raise ValueError


