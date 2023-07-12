"""
Module to run the IR reference building pipeline on WINTER images
"""
import argparse
import logging
import sys
from pathlib import Path

import numpy as np
from astropy.io import fits
from wintertoo.data import winter_fields

from mirar.data import Dataset, Image, ImageBatch
from mirar.data.utils import plot_fits_image
from mirar.paths import BASE_NAME_KEY, core_fields, get_output_dir
from mirar.pipelines.winter.winter_pipeline import WINTERPipeline
from mirar.processors.candidates.utils import get_corners_ra_dec_from_header
from mirar.processors.split import SUB_ID_KEY

logger = logging.getLogger(__name__)

winter_subfields_file = get_output_dir(
    dir_root="cache",
    sub_dir="winter",
).joinpath("WINTER_subfields.txt")
winter_subfields_file.parent.mkdir(parents=True, exist_ok=True)


def dummy_split_image_batch_generator(
    cent_ra: float,
    cent_dec: float,
    fieldid: int,
    nx: int = 3,
    ny: int = 4,
    full_ra_size_deg: float = 1.0,
    full_dec_size_deg: float = 1.2,
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
            hdu.header[SUB_ID_KEY] = subdet
            hdu.header[BASE_NAME_KEY] = f"field{fieldid}_subdet{subdet}.fits"

            image = Image(header=hdu.header, data=data)

            all_images.append(image)
            del hdu

    image_batch = ImageBatch(all_images)

    return image_batch


def write_subfields_file(
    subfields_filename: Path,
    full_ra_size_deg: float = 1.0,
    full_dec_size_deg: float = 1.2,
    nx: int = 3,
    ny: int = 4,
):
    """
    Write a subfields file for the specified fields file
    Args:
        subfields_filename: Path to the subfields file
        full_ra_size_deg: RA-Size of the field in degrees
        full_dec_size_deg: Dec-Size of the field in degrees
        nx: Number of sub-fields in the x direction
        ny: Number of sub-fields in the y direction

    Returns:
        None
    """
    subimg_half_ra_deg = full_ra_size_deg / (2 * nx)
    subimg_half_dec_deg = full_dec_size_deg / (2 * ny)

    with open(subfields_filename, "w", encoding="utf8") as subfield_f:
        subfield_f.write(
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

                with open(subfields_filename, "a", encoding="utf8") as subfield_f:
                    subfield_f.write(
                        f"{int(fieldid)}, {subdetid}, "
                        f"{ra}, {dec}, "
                        f"{ra0_0}, {dec0_0}, {ra0_1}, {dec0_1}, "
                        f"{ra1_0}, {dec1_0}, {ra1_1}, {dec1_1}\n"
                    )


def run_winter_reference_build_pipeline(
    nx: int = 3,
    ny: int = 4,
    full_ra_size_deg: float = 1.0,
    full_dec_size_deg: float = 1.2,
    field_id: int | None = None,
    subdet_id: int | None = None,
    catch_all_errors: bool = True,
):
    """
    Run the reference build pipeline on the winter fields
    Args:
        nx: Number of sub-fields in the x direction
        ny: Number of sub-fields in the y direction
        full_ra_size_deg: Full right ascension size of the field in degrees
        full_dec_size_deg: Full declination size of the field in degrees
        field_id: Run only for this fieldid (for debugging)
        subdet_id: Run only for this subdetid (for debugging)
        catch_all_errors: Catch all errors and continue
    Returns:

    """
    if not winter_subfields_file.exists():
        logger.info(f"Writing subfields file to {winter_subfields_file}")
        write_subfields_file(winter_subfields_file)

    winter_northern_fields = winter_fields[
        (winter_fields["Dec"] > -40) & (winter_fields["Dec"] < 60)
    ].reset_index(drop=True)

    if field_id is not None:
        winter_northern_fields = winter_northern_fields[
            winter_northern_fields["ID"] == field_id
        ].reset_index(drop=True)

    pipeline = WINTERPipeline(night="references", selected_configurations="refbuild")

    res, errorstack = [], []
    for ind, _ in winter_northern_fields.iterrows():
        cent_ra, cent_dec = (
            winter_northern_fields.loc[ind]["RA"],
            winter_northern_fields.loc[ind]["Dec"],
        )
        split_image_batch = dummy_split_image_batch_generator(
            cent_ra=float(cent_ra),
            cent_dec=float(cent_dec),
            fieldid=int(winter_northern_fields.loc[ind]["ID"]),
            nx=nx,
            ny=ny,
            full_ra_size_deg=full_ra_size_deg,
            full_dec_size_deg=full_dec_size_deg,
        )

        subdetids = np.array([x.header[SUB_ID_KEY] for x in split_image_batch])
        if subdet_id is not None:
            split_image_batch = [
                split_image_batch[np.where(subdetids == subdet_id)[0][0]]
            ]
        subdetids = np.array([x.header[SUB_ID_KEY] for x in split_image_batch])

        dataset = Dataset([ImageBatch(x) for x in split_image_batch])

        res, errorstack = pipeline.reduce_images(
            dataset=dataset,
            catch_all_errors=catch_all_errors,
        )

        if len(res) == 0:
            logger.error(
                f"Something went wrong for field {ind}, " f"returned res with len 0"
            )
            continue

        if len(errorstack.failed_images) < len(split_image_batch):
            plots_dir = get_output_dir(
                dir_root="plots",
                sub_dir="winter/references",
            )
            if not plots_dir.exists():
                plots_dir.mkdir(parents=True)

            for res_image in res[0]:
                subdet_id = np.where(subdetids == int(res_image.header[SUB_ID_KEY]))[0][
                    0
                ]
                split_image = split_image_batch[subdet_id]
                corner_wcs_coords = get_corners_ra_dec_from_header(split_image.header)

                plot_savepath = plots_dir / res_image.header[BASE_NAME_KEY].replace(
                    ".fits", ".png"
                )
                if not plot_savepath.exists():
                    plot_fits_image(
                        res_image,
                        savedir=plots_dir,
                        regions_wcs_coords=corner_wcs_coords,
                    )

    return res, errorstack


if __name__ == "__main__":
    logger = logging.getLogger("mirar")
    handler = logging.StreamHandler(sys.stdout)
    formatter = logging.Formatter(
        "%(asctime)s %(name)s [l %(lineno)d] - %(levelname)s - %(message)s"
    )
    handler.setFormatter(formatter)
    logger.addHandler(handler)
    logger.setLevel(logging.DEBUG)

    parser = argparse.ArgumentParser()
    parser.add_argument("-fieldid", type=int, default=None)
    parser.add_argument("-subdetid", type=int, default=None)

    args = parser.parse_args()

    run_winter_reference_build_pipeline(
        field_id=args.fieldid,
        subdet_id=args.subdetid,
    )
