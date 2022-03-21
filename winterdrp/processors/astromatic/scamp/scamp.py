import logging
import os
import numpy as np
import astropy.io.fits
from winterdrp.processors.base_processor import BaseProcessor
from winterdrp.processors.autoastrometry.autoastrometry import run_autoastrometry_single
from winterdrp.paths import output_dir

logger = logging.getLogger(__name__)


class Scamp(BaseProcessor):

    base_key = "scamp"

    def __init__(
            self,
            # write_crosscheck_files: bool = False,
            temp_output_sub_dir: str = "scamp",
            *args,
            **kwargs
    ):
        super(Scamp, self).__init__(*args, **kwargs)
        self.temp_output_sub_dir = temp_output_sub_dir

    def _apply_to_images(
            self,
            images: list[np.ndarray],
            headers: list[astropy.io.fits.Header],
    ) -> tuple[list[np.ndarray], list[astropy.io.fits.Header]]:

        scamp_output_dir = output_dir(self.temp_output_sub_dir, self.night_sub_dir)

        try:
            os.makedirs(scamp_output_dir)
        except OSError:
            pass

        scamp_image_list_path = os.path.join(
            scamp_output_dir,
            os.path.splitext(headers[0]["BASENAME"])[0] + "_scamp_list.txt"
        )

        scamp_img_list = list()

        print(scamp_image_list_path)

        for i, data in enumerate(images):
            header = headers[i]

            temp_path = os.path.join(scamp_output_dir, header["BASENAME"])

            print(temp_path)

            # self.save_fits(data, header, temp_path)

            # run_autoastrometry_single(
            #     img_path=temp_path,
            #     output_dir=scamp_output_dir,
            #     write_crosscheck_files=self.write_crosscheck_files,
            #     overwrite=True
            # )

            # Load up temp path image.header, then delete
            # img, header = self.open_fits(temp_path)
            # images[i] = img
            # headers[i] = header
            # os.remove(temp_path)
            # logger.info(f"Loaded updated header, and deleted temporary file {temp_path}")

        raise

        return images, headers

# def run_scamp(
#         # images: str | list,
#         output_dir: str = ".",
#         # astrefcat_name: str,
#         # images_in_file: bool = False,
#         # config_file: str = os.path.join(scamp_config_dir, "scamp.conf"),
#         keyword_string: str = "",
#         run_local: bool = local_scamp
#
# ):
#     cmd = f"scamp {keyword_string}"
#
#     # if images_in_file:
#     #
#     #     if isinstance(images, list):
#     #         err = "'images_in_file' was set to True, so one file must be provided." \
#     #               f"Instead, 'images' was a list! ({images})."
#     #         logger.error(err)
#     #         raise ValueError(err)
#     #
#     #     cmd += f"@{images}"
#     #
#     # elif isinstance(images, list):
#     #     cmd += ",".join(images)
#     #
#     # else:
#     #     cmd += images
#     #
#     # cmd += f" -c {config_file} -ASTREFCAT_NAME {astrefcat_name}"
#
#     logger.debug(f"Executing '{cmd}'")
#
#     try:
#         execute(cmd, output_dir, local=run_local)
#     except ExecutionError as err:
#         raise ScampExecutionError(err)
