from winterdrp.processors.base_processor import BaseDataframeProcessor
import pandas as pd
import numpy as np
from winterdrp.processors.photometry.utils import make_cutouts
import matplotlib.pyplot as plt
from astropy.stats import sigma_clipped_stats
from matplotlib.patches import Circle
from photutils import CircularAperture, CircularAnnulus, aperture_photometry
from winterdrp.data import SourceBatch


class AperturePhotometry(BaseDataframeProcessor):

    def __init__(self,
                 aper_diameters: list[float] = None,
                 cutout_size_aper_phot: float = 30,
                 bkg_in_diameters: list[float] = None,
                 bkg_out_diameters: list[float] = None,
                 col_suffix_list: list[str] = None,
                 *args,
                 **kwargs):
        if aper_diameters is None:
            aper_diameters = [10.]
        if bkg_in_diameters is None:
            bkg_in_diameters = [25.]
        if bkg_out_diameters is None:
            bkg_out_diameters = [40.]
        super(BaseDataframeProcessor, self).__init__(*args, **kwargs)
        self.cutout_size_aper_phot = cutout_size_aper_phot
        self.aper_diameters = aper_diameters
        self.bkg_in_diameters = bkg_in_diameters
        self.bkg_out_diameters = bkg_out_diameters
        self.col_suffix_list = col_suffix_list
        if self.col_suffix_list is None:
            self.col_suffix_list = self.aper_diameters

    @staticmethod
    def aper_photometry(diff_cutout, diff_unc_cutout, aper_diameter, bkg_in_diameter, bkg_out_diameter,
                        plot=False):

        #     w = WCS(header)
        #     x,y = w.all_world2pix(ra,dec,0)
        x, y = int(diff_cutout.shape[0] / 2), int(diff_cutout.shape[1] / 2)
        if plot:
            fig, ax = plt.subplots()
            m, s = np.nanmean(diff_cutout), np.nanstd(diff_cutout)
            im = ax.imshow(diff_cutout, interpolation='nearest', cmap='gray',
                           vmin=m - s, vmax=m + 10 * s, origin='lower')
            # c = Circle(xy=(x_img, y_img),radius=15)

            c = Circle(xy=(x, y), radius=aper_diameter / 2)
            c1 = Circle(xy=(x, y), radius=bkg_in_diameter / 2)
            c2 = Circle(xy=(x, y), radius=bkg_out_diameter / 2)
            c.set_facecolor('none')
            c.set_edgecolor('red')
            c1.set_facecolor('none')
            c1.set_edgecolor('red')
            c2.set_facecolor('none')
            c2.set_edgecolor('red')
            ax.add_artist(c)
            ax.add_artist(c1)
            ax.add_artist(c2)
            ax.set_xlim(x - 30, x + 30)
            ax.set_ylim(y - 30, y + 30)

        aperture = CircularAperture((x, y), r=aper_diameter)
        annulus_aperture = CircularAnnulus((x, y), r_in=bkg_in_diameter / 2, r_out=bkg_out_diameter / 2)

        annulus_masks = annulus_aperture.to_mask(method='center')
        annulus_data = annulus_masks.multiply(diff_cutout)
        mask = annulus_masks.data
        annulus_data_1d = annulus_data[mask > 0]
        bkg_mean, bkg_median, bkg_std = sigma_clipped_stats(annulus_data_1d, sigma=2)
        bkg = np.zeros(diff_cutout.shape) + bkg_median
        bkg_error = np.zeros(diff_cutout.shape) + bkg_std

        aperture_mask = aperture.to_mask(method='center')
        aperture_unc_data = aperture_mask.multiply(diff_unc_cutout)
        # effective_gain = header['GAIN']
        # error = calc_total_error(data, bkg_error, effective_gain)
        # phot_table = aperture_photometry(diff_cutout - bkg, aperture, error=error)
        # counts_err = phot_table['aperture_sum_err'][0]
        error = np.sqrt(np.sum(aperture_unc_data ** 2))
        phot_table = aperture_photometry(diff_cutout - bkg, aperture)
        counts = phot_table['aperture_sum'][0]
        counts_err = error
        return counts, counts_err

    def _apply_to_candidates(
            self,
            batch: SourceBatch,
    ) -> SourceBatch:

        for source_table in batch:
            candidate_table = source_table.get_data()

            for ind, aper_diam in enumerate(self.aper_diameters):
                bkg_in_diameter = self.bkg_in_diameters[ind]
                bkg_out_diameter = self.bkg_out_diameters[ind]
                suffix = self.col_suffix_list[ind]

                fluxes, fluxuncs = [], []

                for cand_ind in range(len(candidate_table)):
                    row = candidate_table.iloc[cand_ind]
                    ximage, yimage = int(row['X_IMAGE']) - 1, int(row['Y_IMAGE']) - 1
                    diff_filename = row['diffimname']
                    diff_unc_filename = row['diffuncname']
                    diff_cutout = make_cutouts(diff_filename, (ximage, yimage), self.cutout_size_aper_phot)
                    diff_unc_cutout = make_cutouts(diff_unc_filename, (ximage, yimage), self.cutout_size_aper_phot)

                    flux, fluxunc = self.aper_photometry(diff_cutout, diff_unc_cutout, aper_diam, bkg_in_diameter,
                                                         bkg_out_diameter)
                    fluxes.append(flux)
                    fluxuncs.append(fluxunc)
                candidate_table[f'fluxap{suffix}'] = fluxes
                candidate_table[f'fluxuncap{suffix}'] = fluxuncs

                candidate_table[f'magap{suffix}'] = (
                        candidate_table['magzpsci']
                        - 2.5 * np.log10(candidate_table[f'fluxap{suffix}'])
                )
                candidate_table[f'sigmagap{suffix}'] = (
                        1.086 * (candidate_table[f'fluxuncap{suffix}'] / candidate_table[f'fluxap{suffix}']))

            source_table.set_data(candidate_table)
        return batch
