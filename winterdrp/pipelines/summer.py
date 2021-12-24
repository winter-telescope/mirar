from winterdrp.pipelines.base_pipeline import Pipeline


def update_summer(img):

    header = img[0].header

    header["FILTER"] = header["AFT"].split("__")[0]

    img[0].header = header

    return img


class SummerPipeline(Pipeline):

    # Set up elements to use
    bias = False
    dark = True
    flat = False
    stack = True
    dither = True
    astrometry_cal = ("GAIA", 9., 13.)
    photometry_cal = {
        "J": ()
    }

    header_keys = [
        'FILENAME',
        "FIELD",
        "FILTER",
        "EXPTIME",
        "COADDS",
        "UT_SHUT",
    ]

    batch_split_keys = ["FIELD", "FILTER"]

    @staticmethod
    def reformat_raw_data(img):
        header = img[0].header
        # header["FILTER"] = header["AFT"].split("__")[0]
        # print(header)
        header["FILENAME"] = img.split("__")[0]
        img[0].header = header
        return img

    def apply_reduction(self, raw_image_list):
        return
