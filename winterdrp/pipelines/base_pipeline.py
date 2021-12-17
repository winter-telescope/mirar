import logging
import os
import numpy as np
from astropy.io import fits
from winterdrp.paths import cal_output_dir, parse_image_list, reduced_img_dir, reduced_img_path, raw_img_dir
from winterdrp.preprocessing import BiasCalibrator, DarkCalibrator, FlatCalibrator
from winterdrp.io import create_fits


logger = logging.getLogger(__name__)


class Pipeline:

    subclasses = {}

    # Set up elements to use
    bias = True
    dark = True
    flat = True
    stack = True
    dither = True
    astrometry = ("GAIA", 9., 13.)
    photometry_cal = dict()

    def __init__(self, *args, **kwargs):
        self.bias_calibrator = [None, BiasCalibrator(self.open_fits, *args, **kwargs)][self.bias]
        self.dark_calibrator = [None, DarkCalibrator(self.open_fits, *args, **kwargs)][self.dark]
        self.flats_calibrator = [None, FlatCalibrator(self.open_fits, *args, **kwargs)][self.flat]

    def open_fits(self, path):
        img = fits.open(path)
        img = self.reformat_raw_data(img)
        return img

    @staticmethod
    def reformat_raw_data(img):
        return img

    def make_calibration_files(self, sub_dir):

        cal_dict = parse_image_list(sub_dir)

        cal_dir = cal_output_dir(sub_dir)

        # Make calibration directory, unless it already exists

        try:
            os.makedirs(cal_dir)
        except OSError:
            pass

        logger.info(f"Making calibration files for directory {raw_img_dir(sub_dir)}")

        if self.bias:
            make_master_bias(
                cal_dict["bias"],
                cal_dir=cal_dir,
                open_fits=self.open_fits
            )

        if self.dark:
            self.dark_calibrator.make_calibration_files(
                image_list=cal_dict["dark"],
                sub_dir=sub_dir,
                subtract_bias=self.subtract_bias
            )

        if self.flat:
            make_master_flats(
                cal_dict["flats"],
                cal_dir=cal_dir,
                open_fits=self.open_fits,
                subtract_bias=self.subtract_bias,
                subtract_dark=self.subtract_dark
            )

    def reduce_image(self, img, sub_dir=""):

        img = self.subtract_bias(img, sub_dir=sub_dir)
        img = self.subtract_dark(img, sub_dir=sub_dir)
        img = self.divide_flat(img, sub_dir=sub_dir)

        return img

    def preprocess_images(self, sub_dir="", raw_image_list=None, reprocess=True):

        if raw_image_list is None:
            raw_image_list = parse_image_list(sub_dir, group_by_object=False)

        # Try making output directory, unless it exists

        output_dir = reduced_img_dir(sub_dir)

        try:
            os.makedirs(output_dir)
        except OSError:
            pass

        nframes = len(raw_image_list)

        proccessed_list = []

        # Loop over science images

        for i, raw_img_path in enumerate(raw_image_list):

            img_name = os.path.basename(raw_img_path)

            logger.debug(f"Processing image {i + 1}/{nframes} ({img_name})")

            output_path = reduced_img_path(img_name, sub_dir=sub_dir)

            if np.logical_and(os.path.exists(output_path), reprocess is False):
                logger.debug(f"Skipping image {img_name}, because it has already "
                             f"been processed and 'reprocess' is False.")
                continue

            with self.open_fits(raw_img_path) as img:
                header = img[0].header

                if header['OBSTYPE'] not in ['science', "object"]:
                    logger.debug(f'Obstype is not science, skipping {raw_img_path}')
                    continue

                data_redux = self.reduce_image(img, sub_dir=sub_dir)[0].data

                proc_hdu = create_fits(data_redux, header=header, history=None)

                proc_hdu.header['BZERO'] = 0

                # Write the reduced frame to disk

                logger.debug(f"Saving processed image to {output_path}")
                proccessed_list.append(output_path)
                proc_hdu.writeto(output_path, overwrite=True)

        return proccessed_list

    def subtract_bias(self, img, sub_dir=""):

        if self.bias:
            img = self.bias_calibrator.apply_calibration(img, sub_dir=sub_dir)

        return img

    def subtract_dark(self, img, sub_dir=""):

        if self.dark:
            img = self.dark_calibrator.apply_calibration(img, sub_dir=sub_dir)

        return img

    def divide_flat(self, img, sub_dir=""):

        if self.flat:
            img = self.flats_calibrator.apply_calibration(img, sub_dir=sub_dir)

        return img

