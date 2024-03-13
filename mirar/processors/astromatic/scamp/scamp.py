"""
Module to run Scamp
"""

import logging
import os
import shutil
from collections.abc import Callable
from pathlib import Path

import numpy as np
from astropy.io import fits

from mirar.catalog.base_catalog import BaseCatalog
from mirar.data import Image, ImageBatch
from mirar.paths import (
    BASE_NAME_KEY,
    copy_temp_file,
    get_astrometry_keys,
    get_output_dir,
    get_untemp_path,
)
from mirar.processors.astromatic.sextractor.sextractor import (
    SEXTRACTOR_HEADER_KEY,
    check_sextractor_prerequisite,
)
from mirar.processors.base_processor import BaseImageProcessor
from mirar.utils import execute

logger = logging.getLogger(__name__)

SCAMP_HEADER_KEY = "SCMPHEAD"


def run_scamp(
    scamp_list_path: str | Path,
    scamp_config_path: str | Path,
    ast_ref_cat_path: str | Path,
    output_dir: str | Path,
    timeout_seconds: float = 60.0,
):
    """
    Function to run scamp.
    NOTE : By default, the scamp instance here is only designed to run for astrometry.
    This function thus enforces SOLVE_PHOTOM = N as otherwise Scamp behaves weirdly and
    can output FLXSCALE != 1 in the output header. This can cause incosistencies down
    the line, e.g. with Swarp.
    Args:
        scamp_list_path:
        scamp_config_path:
        ast_ref_cat_path:
        output_dir:
        timeout_seconds:

    Returns:

    """
    scamp_cmd = (
        f"scamp @{scamp_list_path} "
        f"-c {scamp_config_path} "
        f"-ASTREFCAT_NAME {ast_ref_cat_path} "
        f"-VERBOSE_TYPE QUIET -SOLVE_PHOTOM N"
    )

    execute(scamp_cmd, output_dir=output_dir, timeout=np.max([60.0, timeout_seconds]))


def write_scamp_header_to_image(image: Image):
    """
    Function to write the scamp header to the image.
    """
    headerfile = image[SCAMP_HEADER_KEY]
    image_header = image.get_header()
    with open(headerfile, "r") as header_f:  # pylint: disable=unspecified-encoding
        scamp_header_data = header_f.read()

    scamp_header = fits.Header()
    scamp_header = scamp_header.fromstring(scamp_header_data, sep="\n")

    # Remove any existing astrometry keywords
    astrometry_keys = get_astrometry_keys()
    for k in astrometry_keys:
        if k in image_header.keys():
            del image_header[k]

    for k in scamp_header:
        if k in ["HISTORY", "COMMENT"]:
            continue

        image_header[k] = scamp_header[k]
    image_header.add_history(f"{scamp_header['HISTORY']}")
    image.set_header(image_header)
    return image


class Scamp(BaseImageProcessor):
    """
    Class for Scamp Processor
    """

    base_key = "scamp"

    def __init__(
        self,
        ref_catalog_generator: Callable[[fits.Header], BaseCatalog],
        scamp_config_path: str,
        temp_output_sub_dir: str = "scamp",
        cache: bool = False,
        copy_scamp_header_to_image: bool = False,
    ):
        super().__init__()
        self.scamp_config = scamp_config_path
        self.ref_catalog_generator = ref_catalog_generator
        self.temp_output_sub_dir = temp_output_sub_dir
        self.cache = cache
        self.copy_scamp_header_to_image = copy_scamp_header_to_image

    def __str__(self) -> str:
        """
        Description of processor
        Returns:

        """
        return (
            f"Processor to apply Scamp to images and calculate astrometry, "
            f"using the config at {self.scamp_config}."
        )

    def get_scamp_output_dir(self) -> Path:
        """
        Function to get scamp output directory
        Returns:
        Path
        """
        return get_output_dir(self.temp_output_sub_dir, self.night_sub_dir)

    def _apply_to_images(self, batch: ImageBatch) -> ImageBatch:
        basenames = [x[BASE_NAME_KEY] for x in batch]
        sort_inds = np.argsort(basenames)
        batch = ImageBatch([batch[i] for i in sort_inds])

        scamp_output_dir = self.get_scamp_output_dir()
        scamp_output_dir.mkdir(parents=True, exist_ok=True)

        ref_catalog = self.ref_catalog_generator(batch[0])

        ref_cat_path = ref_catalog.write_catalog(batch[0], output_dir=scamp_output_dir)

        scamp_image_list_path = scamp_output_dir.joinpath(
            Path(batch[0][BASE_NAME_KEY]).name + "_scamp_list.txt",
        )

        logger.debug(f"Writing file list to {scamp_image_list_path}")

        temp_files = [scamp_image_list_path, ref_cat_path]

        out_files = []

        with open(scamp_image_list_path, "w", encoding="utf8") as img_list_f:
            for image in batch:
                temp_sextractor_cat_path = copy_temp_file(
                    output_dir=scamp_output_dir, file_path=image[SEXTRACTOR_HEADER_KEY]
                )
                img_list_f.write(f"{temp_sextractor_cat_path}\n")
                temp_files += [temp_sextractor_cat_path]

                out_path = Path(
                    os.path.splitext(temp_sextractor_cat_path)[0]
                ).with_suffix(".head")
                out_files.append(out_path)
        num_files = len(batch)
        run_scamp(
            scamp_list_path=scamp_image_list_path,
            scamp_config_path=self.scamp_config,
            ast_ref_cat_path=ref_cat_path,
            output_dir=scamp_output_dir,
            timeout_seconds=30.0 * num_files,
        )

        if not self.cache:
            for path in temp_files:
                logger.debug(f"Deleting temp file {path}")
                path.unlink()

        assert len(batch) == len(out_files)

        for i, out_path in enumerate(out_files):
            image = batch[i]
            new_out_path = get_untemp_path(out_path)
            shutil.move(out_path, new_out_path)
            image[SCAMP_HEADER_KEY] = str(new_out_path).strip()
            if self.copy_scamp_header_to_image:
                image = write_scamp_header_to_image(image)

            batch[i] = image

        return batch

    def check_prerequisites(
        self,
    ):
        check_sextractor_prerequisite(self)
