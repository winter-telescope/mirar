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
    ) -> str:

        cal_dir = cal_output_dir(sub_dir=self.night_sub_dir)

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
            master_dark, _ = self.load_cache_file(self.get_file_path(header))
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
            )

            try:
                os.makedirs(os.path.dirname(master_dark_path))
            except OSError:
                pass

            logger.info(f"Saving stacked 'master dark' "
                        f"combining {n_frames} exposures to {master_dark_path}")

            self.save_fits(master_dark, primary_header, master_dark_path)
