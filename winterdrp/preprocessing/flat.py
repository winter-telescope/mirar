import numpy as np
import os
from winterdrp.io import create_fits
import logging
from winterdrp.preprocessing.base_calibrator import BaseCalibrator
from winterdrp.paths import cal_output_dir

logger = logging.getLogger(__name__)


class FlatCalibrator(BaseCalibrator):

    base_name = "master_flat"

    def __init__(self, open_fits, x_min=0., x_max=np.inf, y_min=0., y_max=np.inf, flat_nan_threshold=np.nan, *args, **kwargs):
        BaseCalibrator.__init__(self, open_fits=open_fits)
        self.x_min = x_min
        self.x_max = x_max
        self.y_min = y_min
        self.y_max = y_max
        self.flat_nan_threshold = flat_nan_threshold

    def get_file_path(self, header, sub_dir=""):
        cal_dir = cal_output_dir(sub_dir=sub_dir)
        filtername = header['FILTER']
        name = f"{self.base_name}_{filtername}.fits"
        return os.path.join(cal_dir, name)

    def apply_calibration(self, img, sub_dir=""):
        data = img[0].data
        header = img[0].header
        master_flat = self.load_calibrator_file(self.get_file_path(header, sub_dir=sub_dir))
        if np.any(master_flat < self.flat_nan_threshold):
            master_flat[master_flat < self.flat_nan_threshold] = np.nan
        img[0].data = data/master_flat
        return img

    def make_calibration_files(self, image_list, sub_dir="", subtract_bias=None, subtract_dark=None, **kwargs):

        image_list = image_list

        logger.info(f'Found {len(image_list)} flat frames')

        with self.open_fits(image_list[0]) as img:
            header = img[0].header

        nx = header['NAXIS1']
        ny = header['NAXIS2']

        filter_list = []

        for flat in image_list:
            with self.open_fits(flat) as img:
                header = img[0].header
            filter_list.append(header['FILTER'])

        image_list = np.array(image_list)

        for f in list(set(filter_list)):

            mask = np.array([x == f for x in image_list])

            cut_flat_list = image_list[mask]

            n_frames = np.sum(mask)

            logger.info(f'Found {n_frames} frames for filer {f}')

            flats = np.zeros((ny, nx, n_frames))

            for i, flat in enumerate(cut_flat_list):
                logger.debug(f'Reading flat {i + 1}/{n_frames}')

                with self.open_fits(flat) as img:
                    data = subtract_dark(subtract_bias(img))[0].data
                    median = np.nanmedian(data[self.x_min:self.x_max, self.y_min:self.y_max])
                    flats[:, :, i] = data / median

            logger.info(f'Median combining {n_frames} flats')

            master_flat = np.nanmedian(flats, axis=2)

            with self.open_fits(image_list[0]) as img:
                primary_header = img[0].header

            proc_hdu = create_fits(master_flat, header=primary_header, history='Stacked flat-fielded')
            # Create a new HDU with the processed image data

            primary_header['BZERO'] = 0

            master_flat_path = self.get_file_path(header, sub_dir=sub_dir)

            logger.info(f"Saving stacked 'master flat' for filter {f} to {master_flat_path}")

            self.save_fits(proc_hdu, master_flat_path)

#
# base_mflat_name = "master_flat"
#
#
# def mflat_name(filtername):
#     return f"{base_mflat_name}_{filtername}.fits"
#
#
# def make_master_flats(flat_list, cal_dir, open_fits, subtract_bias, subtract_dark, xlolim=500, xuplim=3500, ylolim=500, yuplim=3500):
#
#     if len(flat_list) > 0:
#
#         logger.info(f'Found {len(flat_list)} flat frames')
#
#         with open_fits(flat_list[0]) as img:
#             header = img[0].header
#
#         nx = header['NAXIS1']
#         ny = header['NAXIS2']
#
#         master_bias = load_master_bias(cal_dir, open_fits, header)
#
#         filter_list = []
#
#         for flat in flat_list:
#             with open_fits(flat) as img:
#                 header = img[0].header
#             filter_list.append(header['FILTER'])
#
#         flat_list = np.array(flat_list)
#
#         for f in list(set(filter_list)):
#
#             mask = np.array([x == f for x in flat_list])
#
#             cutflatlist = flat_list[mask]
#
#             nframes = np.sum(mask)
#
#             logger.info(f'Found {nframes} frames for filer {f}')
#
#             flats = np.zeros((ny, nx, nframes))
#
#             for i, flat in enumerate(cutflatlist):
#                 logger.debug(f'Reading flat {i+1}/{nframes}')
#
#                 with open_fits(flat) as img:
#                     # data = img[0].data
#                     data = subtract_dark(subtract_bias(img))[0].data
#                     median = np.nanmedian(data[xlolim:xuplim, ylolim:yuplim])
#                     flats[:, :, i] = data/median
#
#             logger.info(f'Median combining {nframes} flats')
#
#             master_flat = np.nanmedian(flats,axis=2)
#
#             with open_fits(flat_list[0]) as img:
#                 primary_header = img[0].header
#
#             proc_hdu = create_fits(master_flat)  # Create a new HDU with the processed image data
#             proc_hdu.header = primary_header       # Copy over the header from the raw file
#             proc_hdu.header.add_history('Stacked flat-fielded')
#
#             primary_header['BZERO'] = 0
#
#             mflat_path = os.path.join(cal_dir, mflat_name(f))
#
#             logger.info(f"Saving stacked 'master flat' for filter {f} to {mflat_path}")
#
#             proc_hdu.writeto(mflat_path, overwrite=True)
#
#         return 0
#
#     else:
#         logger.warning("No flat images provided. Proceeding without flat-fielding correction.")
#
#
# def select_sky_flats():
#     raise NotImplementedError
#
#
# def load_master_flats(cal_dir, open_fits, header=None):
#
#     master_flat_paths = glob(f'{cal_dir}/{base_mflat_name}*.fits')
#
#     if len(master_flat_paths) == 0:
#
#         try:
#             nx = header['NAXIS1']
#             ny = header['NAXIS2']
#
#             master_flats = np.zeros((ny,nx))
#
#             logger.warning("No master flat found. No flat-fielding correction will be applied.")
#
#         except (TypeError, KeyError) as e:
#             err = "No master flat files found, and no header info provided to create a dummy image."
#             logger.error(err)
#             raise FileNotFoundError(err)
#
#     else:
#
#         master_flats = dict()
#
#         for mfpath in master_flat_paths:
#
#             f = os.path.basename(mfpath).split("_")[-1][:-5]
#
#             with open_fits(mfpath) as img:
#                 master_flats[f] = img[0].data
#
#     return master_flats
#
#
# def select_master_flat(all_master_flats, header=None, flat_nan_threshold=0.0):
#
#     if isinstance(all_master_flats, np.ndarray):
#         master_flat = all_master_flats
#
#     elif isinstance(all_master_flats, dict):
#         f = header['FILTER']
#         try:
#             master_flat = all_master_flats[f]
#         except KeyError:
#             err = f"Unrecognised key {f}. Available filters are {all_master_flats.keys()}"
#             logger.error(err)
#             raise KeyError(err)
#     else:
#         err = f"Unrecognised Type for all_master_flats ({type(all_master_flats)}). " \
#               f"Was expecting 'numpy.ndarray' or 'dict'."
#         logger.error(err)
#         raise TypeError(err)
#
#      # Mask pixels below a threshold
#
#     masked_mflat = np.copy(master_flat)
#     if np.any(master_flat < flat_nan_threshold):
#         masked_mflat[master_flat < flat_nan_threshold] = np.nan
#
#     return masked_mflat