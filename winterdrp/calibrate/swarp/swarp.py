import logging
import os
from winterdrp.calibrate.swarp.configs import swarp_config_dir
from winterdrp.utils import execute, ExecutionError

logger = logging.getLogger(__name__)

local_swarp = True


def execute_swarp(cmd, output_dir="."):
    execute(cmd, output_dir=output_dir, local=local_swarp)


class SwarpExecutionError(ExecutionError):
    pass


def run_swarp(
        images: str | list,
        weight_image: str | list = None,
        resample_dir: str = ".",
        imageout_name: str = "stack.fits",
        weightout_name: str = "stack.weight.fits",
        config_file: str = os.path.join(swarp_config_dir, "config.swarp"),
):
    cmd = "swarp "

    if isinstance(images, list):
        cmd += ",".join(images)
    else:
        cmd += images

    if weight_image is not None:

        cmd += " -WEIGHT_IMAGE "

        if isinstance(weight_image, list):
            cmd += ",".join(weight_image)
        else:
            cmd += weight_image

    cmd += f" -c {config_file} " \
           f"-RESAMPLE_DIR {resample_dir} " \
           f"-WEIGHTOUT_NAME {weightout_name} " \
           f"-IMAGEOUT_NAME {imageout_name} "

    logger.debug(f"Executing '{cmd}'")

    try:
        execute_swarp(cmd, output_dir=resample_dir)
    except ExecutionError as err:
        raise SwarpExecutionError(err)
