from astropy.io.fits import HDUList
from winterdrp.pipelines.base_pipeline import Pipeline
import os

wirc_flats_dir = os.path.join(os.path.dirname(os.path.abspath(__file__)))


class WircPipeline(Pipeline):

    astrometry_cal = ("GAIA", 9., 13.)
    photometry_cal = {
        "J": ()
    }

    # Set up elements to use

    image_steps = [
        "dark",
        "flat",
        # "stack",
        # "dither"
    ]

    header_keys = [
        "UTSHUT",
        'OBJECT',
        "FILTER",
        "EXPTIME",
        "COADDS",
    ]

    batch_split_keys = ["OBJECT", "FILTER"]

    @staticmethod
    def reformat_raw_data(
            img: HDUList
    ) -> HDUList:
        header = img[0].header
        header["FILTER"] = header["AFT"].split("__")[0]
        header["OBSCLASS"] = ["calibration", "science"][header["OBSTYPE"] == "object"]
        img[0].header = header
        return img

    # def apply_reduction(self, raw_image_list):
    #     return
