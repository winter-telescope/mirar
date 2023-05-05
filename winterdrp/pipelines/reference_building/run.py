import logging
import sys
from pathlib import Path

import astropy.units as u
import numpy as np
import pandas as pd
from astropy.io import ascii, fits

from winterdrp.data import Dataset, Image, ImageBatch
from winterdrp.paths import BASE_NAME_KEY, core_fields
from winterdrp.pipelines.reference_building.ir_refbuild_pipeline import (
    IRRefBuildPipeline,
)
from winterdrp.references.ukirt import get_query_coordinates_from_header


def dummy_image_generator(
    cent_ra, cent_dec, fieldid, nx=3, ny=4, full_ra_size_deg=1, full_dec_size_deg=1.2
) -> ImageBatch:
    # winter_northern_fields = winter_fields[winter_fields['Dec'] > -30]
    # ind = 0
    # cent_ra, cent_dec = winter_northern_fields.loc[ind]['RA'], \
    #           winter_northern_fields.loc[ind]['Dec']

    subimg_half_ra_deg = full_ra_size_deg / (2 * nx)
    subimg_half_dec_deg = full_dec_size_deg / (2 * ny)

    split_ras, split_decs, all_images = [], [], []
    for iind, i in enumerate(np.arange(-nx + 1, nx + 1, 2)):
        for jind, j in enumerate(np.arange(-ny + 1, ny + 1, 2)):
            ra = cent_ra + subimg_half_ra_deg * i / np.cos(cent_dec * np.pi / 180)
            if ra < 0:
                ra += 360
            dec = cent_dec + subimg_half_dec_deg * j

            hdu = fits.PrimaryHDU()
            hdu.header["NAXIS"] = 2
            hdu.header["NAXIS1"] = 3
            hdu.header["NAXIS2"] = 3
            hdu.header["CTYPE1"] = "RA--TAN"
            hdu.header["CTYPE2"] = "DEC-TAN"

            hdu.header["CRVAL1"] = ra
            hdu.header["CRVAL2"] = dec
            hdu.header["CRPIX1"] = 1.5
            hdu.header["CRPIX2"] = 1.5
            hdu.header["FILTER"] = "J"
            hdu.header["CD1_1"] = (
                2 * subimg_half_ra_deg / hdu.header["NAXIS1"]
            ) / np.cos(cent_dec * np.pi / 180)
            hdu.header["CD2_2"] = (2 * subimg_half_dec_deg) / hdu.header["NAXIS2"]

            data = np.zeros((3, 3))

            for k in core_fields:
                hdu.header[k] = ""

            hdu.header["FIELDID"] = int(fieldid)
            subdet = int(iind * len(np.arange(-ny + 1, ny + 1, 2)) + jind)
            hdu.header["SUBDET"] = subdet
            hdu.header[BASE_NAME_KEY] = f"field{fieldid}_subdet{subdet}.fits"

            image = Image(header=hdu.header, data=data)

            all_images.append(image)

    image_batch = ImageBatch(all_images)

    return image_batch


if __name__ == "__main__":
    fields_file_dir = Path(__file__).parent.joinpath("files")
    fields_file = fields_file_dir.joinpath("WINTER_fields.txt")
    winter_fields = pd.read_csv(fields_file, delim_whitespace=True)

    winter_northern_fields = winter_fields[
        (winter_fields["Dec"] > -40) & (winter_fields["Dec"] < 60)
    ].reset_index(drop=True)

    ind = 0

    cent_ra, cent_dec = (
        winter_northern_fields.loc[ind]["RA"],
        winter_northern_fields.loc[ind]["Dec"],
    )

    split_image_batch = dummy_image_generator(
        cent_ra=cent_ra,
        cent_dec=cent_dec,
        fieldid=winter_northern_fields.loc[ind]["ID"],
    )

    for image in split_image_batch:
        print(get_query_coordinates_from_header(image.header, 4))

    pipeline = IRRefBuildPipeline(
        night="ir_refbuild", selected_configurations="default"
    )

    log = logging.getLogger("winterdrp")
    handler = logging.StreamHandler(sys.stdout)
    formatter = logging.Formatter(
        "%(name)s [l %(lineno)d] - %(levelname)s - %(message)s"
    )
    handler.setFormatter(formatter)
    log.addHandler(handler)
    log.setLevel(logging.DEBUG)

    res, errorstack = pipeline.reduce_images(
        dataset=Dataset([[ImageBatch(x)] for x in split_image_batch][0]),
        catch_all_errors=True,
    )
