import copy
import astropy.io.fits
import numpy as np
import os
import logging
import pandas as pd
from collections.abc import Callable
from winterdrp.processors.base_processor import ProcessorWithCache
from winterdrp.paths import cal_output_dir

logger = logging.getLogger(__name__)


class DarkCalibrator(ProcessorWithCache):
    base_name = "master_dark"
    base_key = "dark"

    def __init__(
            self,
            use_normed_dark: bool = False,
            select_cache_images: Callable[[pd.DataFrame], list] = None,
            *args,
            **kwargs
    ):
        super().__init__(*args, **kwargs)
        self.use_normed_dark = use_normed_dark

        if select_cache_images is None:

            def select_cache_images(x):
                return self.select_from_log(x, "dark")

        self.select_cache_images = select_cache_images

    def get_file_path(
            self,
            header: astropy.io.fits.Header,
            sub_dir: str = ""
    ) -> str:

        cal_dir = cal_output_dir(sub_dir=sub_dir)

        exptime = header['EXPTIME']

        if self.use_normed_dark:
            name = f"{self.base_name}_normed.fits"
        else:
            name = f"{self.base_name}_{exptime:.0f}s.fits"

        return os.path.join(cal_dir, name)

    def _apply_to_images(
            self,
            images: list[np.ndarray],
            headers: list[astropy.io.fits.Header],
    ) -> tuple[list[np.ndarray], list[astropy.io.fits.Header]]:

        for i, data in enumerate(images):
            header = headers[i]
            master_dark, _ = self.load_cache_file(self.get_file_path(header, sub_dir=sub_dir))
            data = data - (master_dark * header["EXPTIME"])
            header["CALSTEPS"] += "dark,"
            images[i] = data
            headers[i] = header

        return images, headers

    def make_cache_files(
            self,
            image_paths: list[str],
    ):
        logger.info(f'Found {len(image_paths)} dark frames')

        img, header = self.open_fits(image_paths[0])

        nx = header['NAXIS1']
        ny = header['NAXIS2']

        logger.info("Making one 'master_dark' for each exposure time.")

        exp_list = []

        for dark in image_paths:
            img, header = self.open_fits(dark)
            exp_list.append(header['EXPTIME'])

        exp_times = sorted(list(set(exp_list)))

        logger.info(f'Found {len(exp_times)} different exposure times: {exp_times}')

        image_paths = np.array(image_paths)

        dark_loop = []

        # Either create one normalised master dark, or loop over each exposure time

        if not self.use_normed_dark:
            for exp in exp_times:
                mask = np.array([x == exp for x in exp_list])

                dark_loop.append((
                    image_paths[mask],
                    exp,
                    "Median stacked dark",
                    )
                )
        else:
            logger.info("Making one additional 'master_dark' combining all exposure times.")
            dark_loop.append((
                image_paths,
                1.0,
                "Median stacked normalised dark",
            ))

        # Loop over each set of darks and create a master_dark

        for (cut_image_list, exp_time, history) in dark_loop:

            n_frames = len(cut_image_list)

            darks = np.zeros((ny, nx, len(cut_image_list)))

            for i, dark in enumerate(cut_image_list):
                img, header = self.open_fits(dark)
                dark_exptime = header['EXPTIME']
                logger.debug(f'Read dark {i + 1}/{n_frames} with exposure time {dark_exptime}')

                # Iteratively apply corrections
                for p in self.preceding_steps:
                    img, header = p.apply(list(img), list(header))

                darks[:, :, i] = img / dark_exptime

            master_dark = np.nanmedian(darks, axis=2)

            primary_header = copy.deepcopy(header)
            primary_header["HISTORY"] = history
            primary_header['EXPTIME'] = exp_time

            master_dark_path = self.get_file_path(
                header=primary_header,
                sub_dir=self.night
            )

            logger.info(f"Saving stacked 'master dark' "
                        f"combining {n_frames} exposures to {master_dark_path}")

            self.save_fits(master_dark, primary_header, master_dark_path)

# base_mdark_name = "master_dark"
#
# def mdark_name(exptime, norm=False):
#     if norm:
#         return f"{base_mdark_name}_normed.fits"
#     else:
#         return f"{base_mdark_name}_{exptime:.0f}s.fits"
#
#
# def make_master_dark(dark_list, cal_dir, open_fits, subtract_bias, make_norm=True):
#
#     if len(dark_list) > 0:
#
#         logger.info(f'Found {len(dark_list)} dark frames')
#
#         with open_fits(dark_list[0]) as img:
#             header = img[0].header
#
#         nx = header['NAXIS1']
#         ny = header['NAXIS2']
#
#         master_bias = load_master_bias(cal_dir, open_fits, header)
#
#         logger.info("Making one 'master_dark' for each exposure time.")
#
#         explist = []
#
#         for dark in dark_list:
#             with open_fits(dark) as img:
#                 header = img[0].header
#             explist.append(header['EXPTIME'])
#
#         exps = sorted(list(set(explist)))
#
#         logger.info(f'Found {len(exps)} different exposure times: {exps}')
#
#         dark_list = np.array(dark_list)
#
#         dark_loop = []
#
#         for exp in exps:
#
#             mask = np.array([x == exp for x in explist])
#
#             dark_loop.append((
#                 dark_list[mask],
#                 exp,
#                 "Median stacked dark",
#                 os.path.join(cal_dir, mdark_name(exp, norm=False))
#             ))
#
#         # Optionally loop over all darks
#
#         if make_norm:
#             logger.info("Making one additional 'master_dark' combining all exposure times.")
#             dark_loop.append((
#                 dark_list,
#                 1.0,
#                 "Median stacked normalised dark",
#                 os.path.join(cal_dir, mdark_name(None, norm=True))
#             ))
#
#         # Loop over each set of darks and create a master_dark
#
#         for (cutdarklist, exptime, history, mdark_path) in dark_loop:
#
#             nframes = len(cutdarklist)
#
#             darks = np.zeros((ny, nx, len(cutdarklist)))
#
#             for i, dark in enumerate(cutdarklist):
#                 with open_fits(dark) as img:
#                     dark_exptime = img[0].header['EXPTIME']
#                     logger.debug(f'Read dark {i + 1}/{nframes} with exposure time {dark_exptime}')
#                     darks[:, :, i] = subtract_bias(img)[0].data * exptime/dark_exptime
#
#             master_dark = np.nanmedian(darks, axis=2)
#
#             with open_fits(dark_list[0]) as img:
#                 primary_header = img[0].header
#
#             proc_hdu = create_fits(master_dark)  # Create a new HDU with the processed image data
#             proc_hdu.header = primary_header       # Copy over the header from the raw file
#
#             proc_hdu.header.add_history(history)
#             proc_hdu.header['EXPTIME'] = exptime
#
#             logger.info(f"Saving stacked 'master dark' combining {nframes} exposures to {mdark_path}")
#
#             proc_hdu.writeto(mdark_path, overwrite=True)
#         return 0
#
#     else:
#         logger.warning("No dark images provided. No master dark created.")
#
#
# def load_master_darks(cal_dir, open_fits, header=None, use_norm=False):
#
#     master_dark_paths = glob(f'{cal_dir}/{base_mdark_name}*.fits')
#
#     mdark_norm_path = mdark_name(None, norm=True)
#
#     if len(master_dark_paths) == 0:
#
#         try:
#             nx = header['NAXIS1']
#             ny = header['NAXIS2']
#
#             master_darks = np.zeros((ny,nx))
#
#             logger.warning("No master dark found. No dark correction will be applied.")
#
#         except (TypeError, KeyError) as e:
#             err = "No master dark files found, and no header info provided to create a dummy image."
#             logger.error(err)
#             raise FileNotFoundError(err)
#
#     elif use_norm:
#
#         with open_fits(mdark_norm_path) as img:
#             master_darks = img[0].data
#
#     else:
#
#         master_darks = dict()
#
#         for mdpath in master_dark_paths:
#             if mdark_norm_path not in mdpath:
#
#                 exp = os.path.basename(mdpath).split("_")[-1][:-6]
#                 with open_fits(mdpath) as img:
#                     master_darks[float(exp)] = img[0].data
#
#     return master_darks
#
#
# def select_master_dark(all_master_darks, header):
#
#     if isinstance(all_master_darks, np.ndarray):
#         master_dark = all_master_darks
#
#     elif isinstance(all_master_darks, dict):
#         exp = header['EXPTIME']
#         try:
#             master_dark = all_master_darks[float(exp)]
#         except KeyError:
#             err = f'Unrecognised key {exp}. Available dark exposure times are {all_master_darks.keys()}'
#             logger.error(err)
#             raise KeyError(err)
#     else:
#         err = f"Unrecognised Type for all_master_darks ({type(all_master_darks)}). " \
#               f"Was expecting 'numpy.ndarray' or 'dict'."
#         logger.error(err)
#         raise TypeError(err)
#
#     return master_dark
#