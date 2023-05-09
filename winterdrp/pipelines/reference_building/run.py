"""
Module to run the IR reference building pipeline on WINTER images
"""
import logging
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from astropy.io import fits

from winterdrp.data import Dataset, Image, ImageBatch
from winterdrp.paths import BASE_NAME_KEY, core_fields
from winterdrp.pipelines.reference_building.ir_refbuild_pipeline import (
    IRRefBuildPipeline,
)

logger = logging.getLogger(__name__)
winter_fields_file_dir = Path(__file__).parent.joinpath("files")
winter_fields_file = winter_fields_file_dir.joinpath("WINTER_fields.txt")
winter_subfields_file = winter_fields_file_dir.joinpath("WINTER_subfields.txt")


def dummy_split_image_batch_generator(
    cent_ra, cent_dec, fieldid, nx=3, ny=4, full_ra_size_deg=1, full_dec_size_deg=1.2
) -> ImageBatch:
    """
    Make a dummy image batch by splitting an image with the specified dimensions
    Args:
        cent_ra: Center RA of the parent image
        cent_dec: Center Dec of the parent image
        fieldid: field ID
        nx: Number of subimages in the x direction
        ny: Number of subimages in the y direction
        full_ra_size_deg: RA-Size of parent image in degrees
        full_dec_size_deg: Dec-Size of parent image in degrees

    Returns:

    """
    subimg_half_ra_deg = full_ra_size_deg / (2 * nx)
    subimg_half_dec_deg = full_dec_size_deg / (2 * ny)

    all_images = []
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
            hdu.header["SUBDETID"] = subdet
            hdu.header[BASE_NAME_KEY] = f"field{fieldid}_subdet{subdet}.fits"

            image = Image(header=hdu.header, data=data)

            all_images.append(image)

    image_batch = ImageBatch(all_images)

    return image_batch


def write_subfields_file(
    fields_filename,
    subfields_filename,
    full_ra_size_deg=1,
    full_dec_size_deg=1.2,
    nx=3,
    ny=4,
):
    """
    Write a subfields file for the specified fields file
    Args:
        fields_filename:
        subfields_filename:
        full_ra_size_deg: RA-Size of the field in degrees
        full_dec_size_deg: Dec-Size of the field in degrees
        nx: Number of sub-fields in the x direction
        ny: Number of sub-fields in the y direction

    Returns:

    """
    winter_fields = pd.read_csv(fields_filename, delim_whitespace=True)
    subimg_half_ra_deg = full_ra_size_deg / (2 * nx)
    subimg_half_dec_deg = full_dec_size_deg / (2 * ny)

    with open(subfields_filename, "w") as f:
        f.write(
            "FieldID, SubdetID, RA_cent, Dec_cent, RA0_0, Dec0_0, "
            "RA0_1, Dec0_1, RA1_0, Dec1_0, RA1_1, Dec1_1\n"
        )
    for _, row in winter_fields.iterrows():
        fieldid = row["ID"]
        cent_ra = row["RA"]
        cent_dec = row["Dec"]
        for iind, i in enumerate(np.arange(-nx + 1, nx + 1, 2)):
            for jind, j in enumerate(np.arange(-ny + 1, ny + 1, 2)):
                subdetid = int(iind * len(np.arange(-ny + 1, ny + 1, 2)) + jind)
                ra = cent_ra + subimg_half_ra_deg * i / np.cos(cent_dec * np.pi / 180)
                if ra < 0:
                    ra += 360
                dec = cent_dec + subimg_half_dec_deg * j

                half_delta_ra = subimg_half_ra_deg / np.cos(cent_dec * np.pi / 180)
                half_delta_dec = subimg_half_dec_deg

                ra0_0 = ra - half_delta_ra
                ra1_0 = ra - half_delta_ra
                ra0_1 = ra + half_delta_ra
                ra1_1 = ra + half_delta_ra
                dec0_0 = dec - half_delta_dec
                dec0_1 = dec - half_delta_dec
                dec1_0 = dec + half_delta_dec
                dec1_1 = dec + half_delta_dec

                with open(subfields_filename, "a") as f:
                    f.write(
                        f"{int(fieldid)}, {subdetid}, "
                        f"{ra}, {dec}, "
                        f"{ra0_0}, {dec0_0}, {ra0_1}, {dec0_1}, "
                        f"{ra1_0}, {dec1_0}, {ra1_1}, {dec1_1}\n"
                    )


def run_winter_reference_build_pipeline(
    winter_fields: pd.DataFrame,
    nx: int = 3,
    ny: int = 4,
    full_ra_size_deg=1,
    full_dec_size_deg=1.2,
):
    """
    Run the reference build pipeline on the winter fields
    Args:
        winter_fields:
        nx:
        ny:
        full_ra_size_deg:
        full_dec_size_deg:

    Returns:

    """
    if not winter_subfields_file.exists():
        logger.info(f"Writing subfields file to {winter_subfields_file}")
        write_subfields_file(winter_fields_file, winter_subfields_file)

    winter_northern_fields = winter_fields[
        (winter_fields["Dec"] > -40) & (winter_fields["Dec"] < 60)
    ].reset_index(drop=True)

    pipeline = IRRefBuildPipeline(night="references", selected_configurations="default")

    res, errorstack = [], []
    for ind, _ in winter_northern_fields.iterrows():
        cent_ra, cent_dec = (
            winter_northern_fields.loc[ind]["RA"],
            winter_northern_fields.loc[ind]["Dec"],
        )

        split_image_batch = dummy_split_image_batch_generator(
            cent_ra=cent_ra,
            cent_dec=cent_dec,
            fieldid=winter_northern_fields.loc[ind]["ID"],
            nx=nx,
            ny=ny,
            full_ra_size_deg=full_ra_size_deg,
            full_dec_size_deg=full_dec_size_deg,
        )

        dataset = Dataset([ImageBatch(x) for x in split_image_batch])

        res, errorstack = pipeline.reduce_images(
            dataset=dataset,
            catch_all_errors=True,
        )

    return res, errorstack


if __name__ == "__main__":
    winter_fields = pd.read_csv(winter_fields_file, delim_whitespace=True)

    logger = logging.getLogger("winterdrp")
    handler = logging.StreamHandler(sys.stdout)
    formatter = logging.Formatter(
        "%(name)s [l %(lineno)d] - %(levelname)s - %(message)s"
    )
    handler.setFormatter(formatter)
    logger.addHandler(handler)
    logger.setLevel(logging.DEBUG)

    run_winter_reference_build_pipeline(winter_fields)
