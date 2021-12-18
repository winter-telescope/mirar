import logging
import numpy as np

logger = logging.getLogger(__name__)


class BaseCalibrator:

    base_name = None

    def __init__(self, open_fits, *args, **kwargs):
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
        return np.copy(img)

    def make_calibration_files(self, image_list, sub_dir="", **kwargs):
        raise NotImplementedError

    def save_fits(self, img, path):
        self.cache[path] = img.data
        logger.info(f"Saving to {path}")
        img.writeto(path, overwrite=True)
