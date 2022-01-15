import logging
import os
from winterdrp.calibrate.scamp.configs import scamp_config_dir
from winterdrp.utils import execute, ExecutionError

logger = logging.getLogger(__name__)

local_scamp = False


class ScampExecutionError(ExecutionError):
    pass


def run_scamp(
        # images: str | list,
        output_dir: str = ".",
        # astrefcat_name: str,
        # images_in_file: bool = False,
        # config_file: str = os.path.join(scamp_config_dir, "scamp.conf"),
        keyword_string: str = "",
        run_local: bool = local_scamp

):
    cmd = f"scamp {keyword_string}"

    # if images_in_file:
    #
    #     if isinstance(images, list):
    #         err = "'images_in_file' was set to True, so one file must be provided." \
    #               f"Instead, 'images' was a list! ({images})."
    #         logger.error(err)
    #         raise ValueError(err)
    #
    #     cmd += f"@{images}"
    #
    # elif isinstance(images, list):
    #     cmd += ",".join(images)
    #
    # else:
    #     cmd += images
    #
    # cmd += f" -c {config_file} -ASTREFCAT_NAME {astrefcat_name}"

    logger.debug(f"Executing '{cmd}'")

    try:
        execute(cmd, output_dir, local=run_local)
    except ExecutionError as err:
        raise ScampExecutionError(err)
