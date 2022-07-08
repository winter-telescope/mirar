import logging

import pandas as pd

from winterdrp.processors.base_processor import BaseCandidateGenerator
import numpy as np
from astropy.io import fits
from winterdrp.processors.astromatic.sextractor.sourceextractor import run_sextractor_dual
from winterdrp.utils.ldac_tools import get_table_from_ldac
from winterdrp.paths import get_output_dir
import io
import gzip
import os

# TODO : Move photometry to its own thing like catalogs, user can choose whichever way they want to do photometry
logger = logging.getLogger(__name__)


class DetectCandidates(BaseCandidateGenerator):
    base_key = "DETCANDS"

    def __init__(
            self,
            cand_det_sextractor_config: str,
            cand_det_sextractor_filter: str,
            cand_det_sextractor_nnw: str,
            cand_det_sextractor_params: str,
            output_sub_dir: str = "candidates",
            *args,
            **kwargs
    ):
        super(DetectCandidates, self).__init__(*args, **kwargs)
        self.output_sub_dir = output_sub_dir
        self.cand_det_sextractor_config = cand_det_sextractor_config
        self.cand_det_sextractor_filter = cand_det_sextractor_filter
        self.cand_det_sextractor_nnw = cand_det_sextractor_nnw
        self.cand_det_sextractor_params = cand_det_sextractor_params

    def get_sub_output_dir(self):
        return get_output_dir(self.output_sub_dir, self.night_sub_dir)

    @staticmethod
    def make_alert_cutouts(
            image_path: str,
            position,
            half_size
    ):
        data = fits.getdata(image_path)
        y_image_size, x_image_size = np.shape(data)
        x, y = position
        # logger.debug(f'{x},{y},{np.shape(data)}')
        if y < half_size:
            cutout = data[0:y + half_size + 1, x - half_size:x + half_size + 1]
            n_pix = half_size - y
            cutout = np.pad(cutout, ((n_pix, 0), (0, 0)), 'constant')

        elif y + half_size + 1 > y_image_size:
            cutout = data[y - half_size: y_image_size, x - half_size: x + half_size + 1]
            n_pix = (half_size + y + 1) - y_image_size
            cutout = np.pad(cutout, ((0, n_pix), (0, 0)), 'constant')

        elif x < half_size:
            cutout = data[y - half_size: y + half_size + 1, 0:x + half_size + 1]
            n_pix = half_size - x
            cutout = np.pad(cutout, ((0, 0), (n_pix, 0)), 'constant')
        elif x + half_size > x_image_size:
            cutout = data[y - half_size:y + half_size + 1, x - half_size:x_image_size]
            n_pix = (half_size + x + 1) - x_image_size
            cutout = np.pad(cutout, ((0, 0), (0, n_pix)), 'constant')
        else:
            cutout = data[y - half_size:y + half_size + 1, x - half_size:x + half_size + 1]
        return cutout

    @staticmethod
    def makebitims(image):
        ######################################################
        # make bit images of the cutouts for the marshal
        #
        # Inputs:
        # image: input image cutout
        #
        # Returns:
        # buf2: a gzipped fits file of the cutout image as
        #  a BytesIO object
        ######################################################

        # open buffer and store image in memory
        buf = io.BytesIO()
        buf2 = io.BytesIO()
        fits.writeto(buf, image)
        with gzip.open(buf2, 'wb') as fz:
            fz.write(buf.getvalue())

        return buf2

    @staticmethod
    def psf_photometry(diff_cutout, diff_unc_cutout, psfmodels):
        numpsfmodels = psfmodels.shape[2]

        chi2s, psf_fluxes, psf_flux_uncs = [], [], []
        for ind in range(numpsfmodels):
            psfmodel = psfmodels[:, :, ind]
            psf_flux = np.sum(psfmodel * diff_cutout) / np.sum(np.square(psfmodel))
            psf_flux_unc = np.sqrt(np.sum(np.square(psfmodel) * np.square(diff_unc_cutout))) / np.sum(
                np.square(psfmodel))
            deg_freedom = np.size(diff_cutout) - 1
            chi2 = np.sum(np.square(diff_cutout - psfmodel * psf_flux) / np.square(diff_unc_cutout)) / deg_freedom

            psf_fluxes.append(psf_flux)
            psf_flux_uncs.append(psf_flux_unc)
            chi2s.append(chi2)

        minchi2_ind = np.argmin(chi2s)
        minchi2 = np.min(chi2s)
        best_fit_psf_flux = psf_fluxes[minchi2_ind]
        best_fit_psf_fluxunc = psf_flux_uncs[minchi2_ind]

        best_fit_psfmodel = psfmodels[:, :, minchi2_ind]
        ys, xs = np.where(best_fit_psfmodel == np.max(best_fit_psfmodel))
        yshift = ys[0] - 6
        xshift = xs[0] - 6

        return best_fit_psf_flux, best_fit_psf_fluxunc, minchi2, xshift, yshift

    @staticmethod
    def make_psf_shifted_array(psf_filename, delta, ngrid=25):
        psf = fits.getdata(psf_filename)
        normpsf = psf / np.sum(psf)
        ngrid = 25
        xs = np.linspace(-2, 2, 5)
        gx, gy = np.meshgrid(xs, xs)
        gx = np.ndarray.flatten(gx)
        gy = np.ndarray.flatten(gy)

        padpsfs = np.zeros((60, 60, ngrid))
        for i in range(ngrid):
            padpsfs[int(10 + gy[i]):int(51 + gy[i]), int(10 + gx[i]):int(51 + gx[i]), i] = normpsf

        normpsfmax = np.max(normpsf)
        x1, x2 = np.where(padpsfs[:, :, 12] == normpsfmax)
        x1 = int(x1)
        x2 = int(x2)
        cutout_size_psf_phot = 20
        psfmodels = padpsfs[x1 - cutout_size_psf_phot:x1 + cutout_size_psf_phot + 1,
                    x2 - cutout_size_psf_phot:x2 + cutout_size_psf_phot + 1]

        return psfmodels

    def generate_candidates_table(self, scorr_catalog_name, sci_resamp_imagename, ref_resamp_imagename, diff_filename,
                                  diff_scorr_filename, diff_psf_filename, diff_unc_filename) -> pd.DataFrame:
        det_srcs = get_table_from_ldac(scorr_catalog_name)

        logger.info(f'Found {len(det_srcs)} candidates in image {diff_filename}.')
        det_srcs['X_IMAGE'] = det_srcs['X_IMAGE'] - 1
        det_srcs['Y_IMAGE'] = det_srcs['Y_IMAGE'] - 1

        scorr_data = fits.getdata(diff_scorr_filename)
        xpeaks, ypeaks = det_srcs['XPEAK_IMAGE'] - 1, det_srcs['YPEAK_IMAGE'] - 1
        scorr_peaks = scorr_data[ypeaks, xpeaks]
        det_srcs['Scorr_peak'] = scorr_peaks

        cutout_size_psf_phot = 20
        cutout_size_display = 40

        psfmodels = self.make_psf_shifted_array(diff_psf_filename, cutout_size_psf_phot)

        det_srcs['psf_flux'] = np.zeros(len(det_srcs))
        det_srcs['psf_fluxunc'] = np.zeros(len(det_srcs))
        det_srcs['minchi2'] = np.zeros(len(det_srcs))
        det_srcs['xshift'] = np.zeros(len(det_srcs))
        det_srcs['yshift'] = np.zeros(len(det_srcs))

        display_sci_ims = []
        display_ref_ims = []
        display_diff_ims = []

        for ind, src in enumerate(det_srcs):
            logger.info(f'Cand # {ind}')
            xpeak, ypeak = int(xpeaks[ind]), int(ypeaks[ind])
            scorr_peak = scorr_peaks[ind]

            diff_cutout = self.make_alert_cutouts(diff_filename, (xpeak, ypeak), cutout_size_psf_phot)
            diff_unc_cutout = self.make_alert_cutouts(diff_unc_filename, (xpeak, ypeak), cutout_size_psf_phot)
            # logger.debug(f'Cutout dimensions {diff_cutout.shape}')
            psf_flux, psf_fluxunc, minchi2, xshift, yshift = self.psf_photometry(diff_cutout, diff_unc_cutout,
                                                                                 psfmodels)

            src['psf_flux'] = psf_flux
            src['psf_fluxunc'] = psf_fluxunc
            src['minchi2'] = minchi2
            src['xshift'] = xshift
            src['yshift'] = yshift

            display_sci_cutout = self.make_alert_cutouts(sci_resamp_imagename, (xpeak, ypeak), cutout_size_display)
            display_ref_cutout = self.make_alert_cutouts(ref_resamp_imagename, (xpeak, ypeak), cutout_size_display)
            display_diff_cutout = self.make_alert_cutouts(diff_filename, (xpeak, ypeak), cutout_size_display)

            display_sci_bit = self.makebitims(display_sci_cutout.astype(np.float32))
            display_ref_bit = self.makebitims(display_ref_cutout.astype(np.float32))
            display_diff_bit = self.makebitims(display_diff_cutout.astype(np.float32))

            display_sci_ims.append(display_sci_bit)
            display_ref_ims.append(display_ref_bit)
            display_diff_ims.append(display_diff_bit)

        det_srcs['SciBitIm'] = display_sci_ims
        det_srcs['RefBitIm'] = display_ref_ims
        det_srcs['DiffBitIm'] = display_diff_ims

        diff_zp = float(fits.getval(diff_filename, 'TMC_ZP'))
        det_srcs['psf_mag'] = diff_zp - 2.5 * np.log10(det_srcs['psf_flux'])
        det_srcs['psf_magerr'] = 1.086 * det_srcs['psf_fluxunc'] / det_srcs['psf_flux']
        det_srcs['ZP'] = diff_zp

        # det_srcs.write(sci_resamp_imagename + '.candidates.dat', format='ascii', overwrite=True)
        det_srcs = det_srcs.to_pandas()
        det_srcs.to_pickle(sci_resamp_imagename + '.candidates.pkl')
        return det_srcs

    def _apply_to_images(
            self,
            images: list[np.ndarray],
            headers: list[fits.Header],
    ) -> pd.DataFrame:

        all_cands_list = []
        for ind, header in enumerate(headers):
            image = images[ind]

            scorr_image_path = os.path.join(self.get_sub_output_dir(), header["DIFFSCR"])
            diff_image_path = os.path.join(self.get_sub_output_dir(), header["DIFFIMG"])
            diff_psf_path = os.path.join(self.get_sub_output_dir(), header["DIFFPSF"])
            diff_unc_path = os.path.join(self.get_sub_output_dir(), header["DIFFUNC"])

            cands_catalog_name = diff_image_path.replace('.fits', '.dets')
            cands_catalog_name = run_sextractor_dual(
                det_image=scorr_image_path,
                measure_image=diff_image_path,
                output_dir=self.get_sub_output_dir(),
                catalog_name=cands_catalog_name,
                config=self.cand_det_sextractor_config,
                parameters_name=self.cand_det_sextractor_params,
                filter_name=self.cand_det_sextractor_filter,
                starnnw_name=self.cand_det_sextractor_nnw,
                gain=1.0
            )

            sci_image_path = os.path.join(self.get_sub_output_dir(), header['BASENAME'])
            ref_image_path = os.path.join(self.get_sub_output_dir(), header['REFIMG'])
            cands_table = self.generate_candidates_table(
                scorr_catalog_name=cands_catalog_name,
                sci_resamp_imagename=sci_image_path,
                ref_resamp_imagename=ref_image_path,
                diff_filename=diff_image_path,
                diff_scorr_filename=scorr_image_path,
                diff_psf_filename=diff_psf_path,
                diff_unc_filename=diff_unc_path
            )

            x_shape, y_shape = image.shape
            cands_table['X_SHAPE'] = x_shape
            cands_table['Y_SHAPE'] = y_shape
            all_cands_list.append(cands_table)

        return pd.concat(all_cands_list)