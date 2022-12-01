import numpy as np
import pandas as pd
from astropy.io import fits

from winterdrp.data import SourceBatch
from winterdrp.processors.base_processor import BaseDataframeProcessor
from winterdrp.processors.photometry.utils import make_cutouts


class PSFPhotometry(BaseDataframeProcessor):
    def __init__(self, cutout_size_psf_phot: float = 20, *args, **kwargs):
        super(BaseDataframeProcessor, self).__init__(*args, **kwargs)
        self.cutout_size_psf_phot = cutout_size_psf_phot

    @staticmethod
    def psf_photometry(diff_cutout, diff_unc_cutout, psfmodels):
        numpsfmodels = psfmodels.shape[2]

        chi2s, psf_fluxes, psf_flux_uncs = [], [], []
        for ind in range(numpsfmodels):
            psfmodel = psfmodels[:, :, ind]
            psf_flux = np.sum(psfmodel * diff_cutout) / np.sum(np.square(psfmodel))
            psf_flux_unc = np.sqrt(
                np.sum(np.square(psfmodel) * np.square(diff_unc_cutout))
            ) / np.sum(np.square(psfmodel))
            deg_freedom = np.size(diff_cutout) - 1
            chi2 = (
                np.sum(
                    np.square(diff_cutout - psfmodel * psf_flux)
                    / np.square(diff_unc_cutout)
                )
                / deg_freedom
            )

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
    def make_psf_shifted_array(psf_filename):
        psf = fits.getdata(psf_filename)
        normpsf = psf / np.sum(psf)
        ngrid = 81
        xs = np.linspace(-4, 4, 9)
        gx, gy = np.meshgrid(xs, xs)
        gx = np.ndarray.flatten(gx)
        gy = np.ndarray.flatten(gy)

        padpsfs = np.zeros((60, 60, ngrid))
        for i in range(ngrid):
            padpsfs[
                int(10 + gy[i]) : int(51 + gy[i]), int(10 + gx[i]) : int(51 + gx[i]), i
            ] = normpsf

        normpsfmax = np.max(normpsf)
        x1, x2 = np.where(padpsfs[:, :, 12] == normpsfmax)
        x1 = int(x1)
        x2 = int(x2)
        cutout_size_psf_phot = 20
        psfmodels = padpsfs[
            x1 - cutout_size_psf_phot : x1 + cutout_size_psf_phot + 1,
            x2 - cutout_size_psf_phot : x2 + cutout_size_psf_phot + 1,
        ]

        return psfmodels

    def _apply_to_candidates(
        self,
        batch: SourceBatch,
    ) -> SourceBatch:

        for source_table in batch:
            candidate_table = source_table.get_data()

            fluxes, fluxuncs, minchi2s, xshifts, yshifts = [], [], [], [], []

            for ind in range(len(candidate_table)):
                row = candidate_table.iloc[ind]
                xpeak, ypeak = row["xpeak"], row["ypeak"]
                diff_filename = row["diffimname"]
                diff_psf_filename = row["diffpsfname"]
                diff_unc_filename = row["diffuncname"]
                diff_cutout = make_cutouts(
                    diff_filename, (xpeak, ypeak), self.cutout_size_psf_phot
                )
                diff_unc_cutout = make_cutouts(
                    diff_unc_filename, (xpeak, ypeak), self.cutout_size_psf_phot
                )
                psfmodels = self.make_psf_shifted_array(
                    diff_psf_filename, self.cutout_size_psf_phot
                )
                flux, fluxunc, minchi2, xshift, yshift = self.psf_photometry(
                    diff_cutout, diff_unc_cutout, psfmodels
                )
                fluxes.append(flux)
                fluxuncs.append(fluxunc)
                minchi2s.append(minchi2)
                xshifts.append(xshift)
                yshifts.append(yshift)

            candidate_table["psf_flux"] = fluxes
            candidate_table["psf_fluxunc"] = fluxuncs
            candidate_table["chipsf"] = minchi2s
            candidate_table["xshift"] = xshifts
            candidate_table["yshift"] = yshifts
            candidate_table["magpsf"] = candidate_table["magzpsci"] - 2.5 * np.log10(
                candidate_table["psf_flux"]
            )
            candidate_table["sigmapsf"] = (
                1.086 * candidate_table["psf_fluxunc"] / candidate_table["psf_flux"]
            )

            source_table.set_data(candidate_table)
        return batch
