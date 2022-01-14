import logging
import os
from winterdrp.calibrate.scamp.configs import scamp_config_dir
from winterdrp.calibrate.sextractor.sourceextractor import execute_sextractor

logger = logging.getLogger(__name__)


def run_scamp(
        images: str | list,
        output_dir: str,
        astrefcat_name: str,
        images_in_file: bool = False,
        config_file: str = os.path.join(scamp_config_dir, "scamp.conf"),

):
    cmd = "scamp "

    if images_in_file:

        if isinstance(images, list):
            err = "'images_in_file' was set to True, so one file must be provided." \
                  f"Instead, 'images' was a list! ({images})."
            logger.error(err)
            raise ValueError(err)

        cmd += f"@{images}"

    elif isinstance(images, list):
        cmd += ",".join(images)

    else:
        cmd += images

    cmd += f" -c {config_file} -ASTREFCAT_NAME {astrefcat_name}"

    # logger.debug(f"Using '{['local', 'docker'][sextractor_cmd == local_sextractor]}' "
    #              f"sextractor installation to run `{cmd}`")

    logger.debug(f"Executing '{cmd}'")
    execute_sextractor(cmd, output_dir)
