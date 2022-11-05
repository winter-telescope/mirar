import logging
import os
import numpy as np
import astropy.io.fits
from winterdrp.processors.base_processor import ProcessorWithCache, BaseImageProcessor
from winterdrp.paths import get_output_dir, copy_temp_file, get_temp_path, get_untemp_path
from winterdrp.utils import execute
from winterdrp.catalog.base_catalog import BaseCatalog
from collections.abc import Callable
from winterdrp.processors.astromatic.sextractor.sextractor import Sextractor, sextractor_header_key
import shutil

logger = logging.getLogger(__name__)

scamp_header_key = "SCMPHEAD"


def run_scamp(
        scamp_list_path: str,
        scamp_config_path: str,
        ast_ref_cat_path: str,
        output_dir: str
):
    scamp_cmd = f"scamp @{scamp_list_path} " \
                f"-c {scamp_config_path} " \
                f"-ASTREFCAT_NAME {ast_ref_cat_path} " \
                f"-VERBOSE_TYPE QUIET "

    execute(scamp_cmd, output_dir=output_dir)


def get_scamp_output_head_path(
        cat_path: str
) -> str:
    return os.path.splitext(cat_path)[0] + ".head"


class Scamp(BaseImageProcessor):

    base_key = "scamp"

    def __init__(
            self,
            ref_catalog_generator: Callable[[astropy.io.fits.Header], BaseCatalog],
            scamp_config_path: str,
            temp_output_sub_dir: str = "scamp",
            *args,
            **kwargs
    ):
        super(Scamp, self).__init__(*args, **kwargs)
        self.scamp_config = scamp_config_path
        self.ref_catalog_generator = ref_catalog_generator
        self.temp_output_sub_dir = temp_output_sub_dir

    def __str__(self) -> str:
        return f"Processor to apply Scamp to images, calculating more precise astrometry."

    def get_scamp_output_dir(self):
        return get_output_dir(self.temp_output_sub_dir, self.night_sub_dir)

    def _apply_to_images(
            self,
            images: list[np.ndarray],
            headers: list[astropy.io.fits.Header],
    ) -> tuple[list[np.ndarray], list[astropy.io.fits.Header]]:

        scamp_output_dir = self.get_scamp_output_dir()

        try:
            os.makedirs(scamp_output_dir)
        except OSError:
            pass

        ref_catalog = self.ref_catalog_generator(headers[0])

        cat_path = copy_temp_file(
            output_dir=scamp_output_dir,
            file_path=ref_catalog.write_catalog(
                headers[0],
                output_dir=scamp_output_dir
            )
        )

        scamp_image_list_path = os.path.join(
            scamp_output_dir,
            os.path.splitext(headers[0]["BASENAME"])[0] + "_scamp_list.txt"
        )

        logger.info(f"Writing file list to {scamp_image_list_path}")

        temp_files = [scamp_image_list_path]

        out_files = []

        with open(scamp_image_list_path, "w") as f:
            for i, data in enumerate(images):
                header = headers[i]

                temp_cat_path = copy_temp_file(
                    output_dir=scamp_output_dir,
                    file_path=header[sextractor_header_key]
                )

                temp_img_path = get_temp_path(scamp_output_dir, header["BASENAME"])
                self.save_fits(data, header, temp_img_path)
                temp_mask_path = self.save_mask(data, header, temp_img_path)
                f.write(f"{temp_cat_path}\n")
                temp_files += [temp_cat_path, temp_img_path, temp_mask_path]

                out_path = get_scamp_output_head_path(temp_cat_path)

                out_files.append(out_path)

        run_scamp(
            scamp_list_path=scamp_image_list_path,
            scamp_config_path=self.scamp_config,
            ast_ref_cat_path=cat_path,
            output_dir=scamp_output_dir,
        )

        for path in temp_files:
            logger.debug(f"Deleting temp file {path}")
            os.remove(path)

        assert len(headers) == len(out_files)

        for i, out_path in enumerate(out_files):
            header = headers[i]
            new_out_path = get_untemp_path(out_path)
            shutil.move(out_path, new_out_path)
            header[scamp_header_key] = new_out_path
            logger.info(f"Saved to {new_out_path}")
            headers[i] = header

        return images, headers

    def check_prerequisites(
            self,
    ):
        check = np.sum([isinstance(x, Sextractor) for x in self.preceding_steps])
        if check < 1:
            err = f"{self.__module__} requires {Sextractor} as a prerequisite. " \
                  f"However, the following steps were found: {self.preceding_steps}."
            logger.error(err)
            raise ValueError


