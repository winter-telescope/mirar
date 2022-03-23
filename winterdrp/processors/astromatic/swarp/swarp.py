import logging
import os
import numpy as np
import astropy.io.fits
from winterdrp.processors.base_processor import BaseProcessor
from winterdrp.paths import get_output_dir, copy_temp_file, get_temp_path, get_untemp_path, get_mask_path
from winterdrp.utils import execute
from winterdrp.catalog.base_catalog import BaseCatalog
from collections.abc import Callable
from winterdrp.processors.astromatic.scamp.scamp import Scamp, scamp_header_key
import shutil

logger = logging.getLogger(__name__)


# def run_swarp(
#         scamp_list_path: str,
#         scamp_config_path: str,
#         ast_ref_cat_path: str,
#         output_dir: str
# ):
#     swarp_cmd = f"scamp @{scamp_list_path} " \
#                 f"-c {scamp_config_path} " \
#                 f"-ASTREFCAT_NAME {ast_ref_cat_path} " \
#                 # f"-VERBOSE_TYPE QUIET "
#
#     execute(swarp_cmd, output_dir=output_dir)

def run_swarp(
        stack_list_path: str,
        swarp_config_path: str,
        out_path: str,
        weight_list_path: str = None,
        weight_out_path: str = None
):  # resample and stack images with swarp
    """
    Resample and stack given images with swarp

    Parameters
    ----------
    stack_list_path : string
        Name of file containing the names of files to be stacked
        One file name per line
    weight_list_path : string
        Name of file containing the names of weight files to be stacked
        One file name per line
    swarp_config_path: str
        Path of Swarp config file
    out_path : string
        Path of stacked output file
    weight_out_path: str
        Path of output weight image
    """

    swarp_command = f'swarp -c {swarp_config_path} ' \
                    f'@{stack_list_path} ' \
                    f'-IMAGEOUT_NAME {out_path}' \
                    f' -RESAMPLE_DIR {os.path.dirname(out_path)} '

    if weight_list_path is not None:
        swarp_command += f'-WEIGHT_IMAGE @{weight_list_path} '

    if weight_out_path is not None:
        swarp_command += f' -WEIGHTOUT_NAME {weight_out_path} '

    print(swarp_command)

    execute(swarp_command)


# def get_scamp_output_head_path(
#         cat_path: str
# ) -> str:
#     return os.path.splitext(cat_path)[0] + ".head"


class Swarp(BaseProcessor):

    base_key = "swarp"

    def __init__(
            self,
            swarp_config_path: str,
            temp_output_sub_dir: str = "swarp",
            *args,
            **kwargs
    ):
        super(Swarp, self).__init__(*args, **kwargs)
        self.swarp_config = swarp_config_path
        self.temp_output_sub_dir = temp_output_sub_dir

    def get_swarp_output_dir(self):
        return get_output_dir(self.temp_output_sub_dir, self.night_sub_dir)

    def _apply_to_images(
            self,
            images: list[np.ndarray],
            headers: list[astropy.io.fits.Header],
    ) -> tuple[list[np.ndarray], list[astropy.io.fits.Header]]:

        swarp_output_dir = self.get_swarp_output_dir()

        try:
            os.makedirs(swarp_output_dir)
        except OSError:
            pass

        swarp_image_list_path = os.path.join(
            swarp_output_dir,
            os.path.splitext(headers[0]["BASENAME"])[0] + "_swarp_img_list.txt"
        )

        swarp_weight_list_path = os.path.join(
            swarp_output_dir,
            os.path.splitext(headers[0]["BASENAME"])[0] + "_swarp_weight_list.txt"
        )

        logger.info(f"Writing file list to {swarp_image_list_path}")

        temp_files = [swarp_image_list_path]

        output_image_path = os.path.join(
            swarp_output_dir,
            os.path.splitext(headers[0]["BASENAME"])[0] + "_stack.fits"
        )

        out_files = []

        with open(swarp_image_list_path, "w") as f, open(swarp_weight_list_path, "w") as g:
            for i, data in enumerate(images):
                header = headers[i]

                print(header[scamp_header_key])

                temp_head_path = copy_temp_file(
                    output_dir=swarp_output_dir,
                    file_path=str(header[scamp_header_key])
                )

                temp_img_path = get_temp_path(swarp_output_dir, header["BASENAME"])
                self.save_fits(data, header, temp_img_path)

                temp_mask_path = self.save_mask(data, header, temp_img_path)

                f.write(f"{temp_img_path}\n")
                g.write(f"{temp_mask_path}\n")
                temp_files += [temp_head_path, temp_img_path, temp_mask_path]

                # out_files.append(out_path)

        run_swarp(
            stack_list_path=swarp_image_list_path,
            swarp_config_path=self.swarp_config,
            out_path=output_image_path,
            weight_list_path=swarp_weight_list_path,
            weight_out_path=output_image_path.replace(".fits", ".weight.fits")
        )

        raise

        # for path in temp_files:
        #     logger.debug(f"Deleting temp file {path}")
        #     os.remove(path)
        #
        # assert len(headers) == len(out_files)
        #
        # for i, out_path in enumerate(out_files):
        #     header = headers[i]
        #     new_out_path = get_untemp_path(out_path)
        #     shutil.move(out_path, new_out_path)
        #     header["HEADPATH"] = new_out_path
        #     logger.info(f"Saved to {new_out_path}")

        return images, headers

    def check_prerequisites(
            self,
    ):
        check = np.sum([isinstance(x, Scamp) for x in self.preceding_steps])
        if check < 1:
            err = f"{self.__module__} requires {Scamp} as a prerequisite. " \
                  f"However, the following steps were found: {self.preceding_steps}."
            logger.error(err)
            raise ValueError

    # @classmethod
    # def single_catalog(
    #         cls,
    #         catalog: BaseCatalog,
    #         *args,
    #         **kwargs
    # ):
    #
    #     def get_catalog(
    #             header: astropy.io.fits.Header
    #     ) -> BaseCatalog:
    #         return catalog
    #
    #     return cls(get_catalog=get_catalog, *args, **kwargs)
