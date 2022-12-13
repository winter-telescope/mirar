import logging
import os
from pathlib import Path

import astropy.io.fits
import numpy as np

from winterdrp.data import Image
from winterdrp.paths import BASE_NAME_KEY, COADD_KEY, PROC_HISTORY_KEY
from winterdrp.utils.fits_tools import save_HDU_as_fits

logger = logging.getLogger(__name__)


class BaseReferenceGenerator:
    @property
    def abbreviation(self):
        raise NotImplementedError()

    def __init__(self, filter_name: str):
        self.filter_name = filter_name

    @staticmethod
    def get_reference(image: Image) -> astropy.io.fits.PrimaryHDU:
        raise NotImplementedError()

    def write_reference(self, image: Image, output_dir: str) -> Path:

        base_name = os.path.basename(image[BASE_NAME_KEY])
        logger.debug(f"Base name is {base_name}")

        refHDU = self.get_reference(image)

        output_path = Path(
            self.get_output_path(output_dir, base_name).replace(".fits", "")
            + "_ref.fits"
        )

        # This is because Swarp requires the COADDS keyword. I am setting it to zero manually
        if "COADDS" not in refHDU.header.keys():
            logger.debug("Setting COADDS to 1")
            refHDU.header[COADD_KEY] = 1
        if "CALSTEPS" not in refHDU.header.keys():
            logger.debug("Setting CALSTEPS to blank")
            refHDU.header[PROC_HISTORY_KEY] = ""

        # Remove if needed
        output_path.unlink(missing_ok=True)

        logger.info(f"Saving reference image to {output_path}")
        refHDU.header[BASE_NAME_KEY] = os.path.basename(output_path)
        refHDU.data[refHDU.data == 0] = np.nan
        save_HDU_as_fits(refHDU, output_path)

        return output_path

    @staticmethod
    def get_output_path(output_dir: str, base_name: str) -> str:
        return os.path.join(output_dir, base_name)
