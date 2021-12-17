import logging

logger = logging.getLogger(__name__)
#
#
# def make_calibration_images(object_list, subdir="", flat_method=None):
#
#     cal_dir = cal_output_dir(subdir)
#
#     # Make calibration directory, unless it already exists
#
#     try:
#         os.makedirs(cal_dir)
#     except OSError:
#         pass
#
#     make_master_bias(object_list["bias"], cal_dir=cal_dir)
#     make_master_dark(object_list["dark"], cal_dir=cal_dir)
#     make_master_flats(object_list["dark"], cal_dir=cal_dir)
#
#     # if flat_method is None:
#     #     logger.warning("Flat-fielding method not specified. No master flats will be built.")
#     #
#     # else:
#     #     if flat_method in ["dome"]:
#     #         flats = object_list["flat"]
#     #     elif flat_method in ["sky"]:
#     #         flats = select_sky_flats(subdir)
#     #     else:
#     #         msg = f"Selected 'flat_method' ({flat_method}) not recognised." \
#     #               f"Please specify 'dome' or 'sky, or None to skip this step."
#     #
#     #         logger.error(msg)
#     #         raise ValueError(msg)


class BaseCalibrator:

    base_name = None

    def __init__(self, open_fits, **kwargs):
        self.cache = dict()
        self.open_fits = open_fits

    def apply_calibration(self, img):
        raise NotImplementedError

    @staticmethod
    def get_file_path(header, sub_dir=""):
        raise NotImplementedError

    def load_calibrator_file(self, path):

        if path in self.cache:
            img = self.cache[path]
        else:
            img = self.open_fits(path)
            self.cache[path] = img
        return img

    def make_calibration_files(self, image_list, sub_dir="", **kwargs):
        raise NotImplementedError

    def save_fits(self, img, path):
        self.cache[path] = img.data
        logger.info(f"Saving to {path}")
        img.writeto(path, overwrite=True)
