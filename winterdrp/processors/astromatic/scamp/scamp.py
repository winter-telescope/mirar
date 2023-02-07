"""
Module for running Scamp
"""
import logging
import os
import shutil
from collections.abc import Callable
from pathlib import Path

import astropy.io.fits
import numpy as np

from winterdrp.catalog.base_catalog import BaseCatalog
from winterdrp.data import ImageBatch
from winterdrp.paths import (
    BASE_NAME_KEY,
    copy_temp_file,
    get_output_dir,
    get_temp_path,
    get_untemp_path,
)
from winterdrp.processors.astromatic.sextractor.sextractor import (
    SEXTRACTOR_HEADER_KEY,
    Sextractor,
)
from winterdrp.processors.base_processor import BaseImageProcessor
from winterdrp.utils import execute

logger = logging.getLogger(__name__)

SCAMP_HEADER_KEY = "SCMPHEAD"


def run_scamp(
    scamp_list_path: str | Path,
    scamp_config_path: str | Path,
    ast_ref_cat_path: str | Path,
    output_dir: str | Path,
    timeout: float,
):
    """
    Function to run scamp command

    :param scamp_list_path: path to text file containing a list of images
    :param scamp_config_path: path to config file for scamp
    :param ast_ref_cat_path: path to astrometry reference catalog
    :param output_dir: directory to output to
    :param timeout: timeout in seconds
    :return: None
    """
    scamp_cmd = (
        f"scamp @{scamp_list_path} "
        f"-c {scamp_config_path} "
        f"-ASTREFCAT_NAME {ast_ref_cat_path} "
        f"-VERBOSE_TYPE QUIET "
    )

    execute(scamp_cmd, output_dir=output_dir, timeout=timeout)


class Scamp(BaseImageProcessor):
    """
    Processor to run Scamp
    """

    base_key = "scamp"

    def __init__(
        self,
        ref_catalog_generator: Callable[[astropy.io.fits.Header], BaseCatalog],
        scamp_config_path: str,
        temp_output_sub_dir: str = "scamp",
        timeout: float = 120.0,
    ):
        super().__init__()
        self.scamp_config = scamp_config_path
        self.ref_catalog_generator = ref_catalog_generator
        self.temp_output_sub_dir = temp_output_sub_dir
        self.timeout = timeout

    def __str__(self) -> str:
        return (
            "Processor to apply Scamp to images, calculating more precise astrometry."
        )

    def get_scamp_output_dir(self) -> Path:
        """
        Get default output directory

        :return: path
        """
        return get_output_dir(self.temp_output_sub_dir, self.night_sub_dir)

    def _apply_to_images(self, batch: ImageBatch) -> ImageBatch:
        scamp_output_dir = self.get_scamp_output_dir()
        scamp_output_dir.mkdir(parents=True, exist_ok=True)

        ref_catalog = self.ref_catalog_generator(batch[0])

        cat_path = copy_temp_file(
            output_dir=scamp_output_dir,
            file_path=ref_catalog.write_catalog(batch[0], output_dir=scamp_output_dir),
        )

        scamp_image_list_path = scamp_output_dir.joinpath(
            Path(batch[0][BASE_NAME_KEY]).name + "_scamp_list.txt",
        )

        logger.debug(f"Writing file list to {scamp_image_list_path}")

        temp_files = [scamp_image_list_path]

        out_files = []

        with open(scamp_image_list_path, "w", encoding="utf8") as txt_file:
            for image in batch:
                temp_cat_path = copy_temp_file(
                    output_dir=scamp_output_dir, file_path=image[SEXTRACTOR_HEADER_KEY]
                )

                temp_img_path = get_temp_path(scamp_output_dir, image[BASE_NAME_KEY])
                self.save_fits(image, temp_img_path)
                temp_mask_path = self.save_weight_image(image, temp_img_path)
                txt_file.write(f"{temp_cat_path}\n")
                temp_files += [temp_cat_path, temp_img_path, temp_mask_path]

                out_files.append(
                    Path(os.path.splitext(temp_cat_path)[0]).with_suffix(".head")
                )

        run_scamp(
            scamp_list_path=scamp_image_list_path,
            scamp_config_path=self.scamp_config,
            ast_ref_cat_path=cat_path,
            output_dir=scamp_output_dir,
            timeout=self.timeout,
        )

        for path in temp_files:
            logger.debug(f"Deleting temp file {path}")
            path.unlink()

        assert len(batch) == len(out_files)

        for i, out_path in enumerate(out_files):
            image = batch[i]
            new_out_path = get_untemp_path(out_path)
            shutil.move(out_path, new_out_path)
            image[SCAMP_HEADER_KEY] = str(new_out_path).strip()
            logger.debug(f"Saved to {new_out_path}")
            batch[i] = image

        return batch

    def check_prerequisites(
        self,
    ):
        check = np.sum([isinstance(x, Sextractor) for x in self.preceding_steps])
        if check < 1:
            err = (
                f"{self.__module__} requires {Sextractor} as a prerequisite. "
                f"However, the following steps were found: {self.preceding_steps}."
            )
            logger.error(err)
            raise ValueError(err)
