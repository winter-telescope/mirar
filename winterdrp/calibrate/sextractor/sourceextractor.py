import os
import logging
from pathlib import Path
from winterdrp.paths import calibration_config_dir
from winterdrp.utils import execute, ExecutionError


logger = logging.getLogger(__name__)

# sextractor_cmd = os.getenv("SEXTRACTOR_CMD")

default_saturation = 1.e10
default_config = os.path.join(calibration_config_dir, 'astrom.sex')


class SextractorError(ExecutionError):
    pass


# Either run sextractor locally or on docker

local_sextractor = True

# Functions to parse commands and generate appropriate sextractor files

def parse_checkimage(
    checkimage_type: str | list = None,
    checkimage_name: str | list = None,
    image: str = None,
):
    """Function to parse the "checkimage" component of Sextractor configuration.

    Parameters
    ----------
    checkimage_type: The 'CHECKIMAGE_TYPE' files for sextractor. The default is None. To quote sextractor,
    available types are: 'NONE, BACKGROUND, BACKGROUND_RMS, MINIBACKGROUND, MINIBACK_RMS, -BACKGROUND,
    FILTERED, OBJECTS, -OBJECTS, SEGMENTATION, or APERTURES'. Multiple arguments should be specified in a list.
    checkimage_name: The name(s) of the checkput images to be output.
    image: The name of the image in question. If specified, the name of each checkimage will include the
    name of the original base image.

    Returns
    -------
    cmd: A string containing the partial sextractor command relating to checkimages. The default is an empty string.
    """
    if isinstance(checkimage_type, str):
        checkimage_type = list(checkimage_type)

    if isinstance(checkimage_name, str):
        checkimage_name = list(checkimage_name)

    cmd = ""

    if checkimage_type is not None:

        cmd = f"-CHECKIMAGE_TYPE {','.join(checkimage_type)} "

        if checkimage_name is not None:
            if not len(checkimage_type) == len(checkimage_name):
                err = f"Number of checkimage types {len(checkimage_type)} does not " \
                      f"match number of checkimage names {len(checkimage_name)}. " \
                      f"These values must be equal. The following types were given: {checkimage_type}, " \
                      f"and the following names were given: {checkimage_name}"
                logger.error(err)
                raise ValueError(err)
            else:
                cmd += f"-CHECKIMAGE_NAME {','.join(checkimage_name)}"

        else:
            if image is not None:
                base_name = f'{os.path.basename(image).split(".")[0]}_'
            else:
                base_name = ""

            cmd += " -CHECKIMAGE_NAME " + ",".join([
                f"{base_name}check_{x.lower()}.fits" for x in checkimage_type
            ])

        cmd += " "

    return cmd


def run_sextractor(
        images: str | list,
        output_dir: str,
        *args,
        **kwargs
):

    if not isinstance(images, list):
        images = [images]

    # Make output directory if it doesn't exist

    try:
        os.makedirs(output_dir)
    except OSError:
        pass

    for img in images:
        run_sextractor_single(img, *args, **kwargs)


def run_sextractor_single(
        img: str,
        output_dir: str,
        catalog_name: str = None,
        config: str = default_config,
        parameters_name: str = None,
        filter_name: str = None,
        starnnw_name: str = None,
        saturation: float = default_saturation,
        weight_image: str = None,
        verbose_type: str = "QUIET",
        checkimage_name: str | list = None,
        checkimage_type: str | list = None,
        gain: float = 0.0,
        run_local: bool = local_sextractor
):
    if parameters_name is None:
        parameters_name = os.path.join(calibration_config_dir, 'astrom.param')
        msg = "No parameters name was passed. Only paths explicitly passed as arguments are used. " \
              f"Using a default parameter file of {filter_name}"
        logger.warning(msg)

    if filter_name is None:
        filter_name = os.path.join(calibration_config_dir, 'default.conv')
        msg = "No filter name was passed. Only paths explicitly passed as arguments are used. " \
              f"Using a default filter file of {filter_name}"
        logger.warning(msg)

    if starnnw_name is None:
        starnnw_name = os.path.join(calibration_config_dir, 'default.nnw')
        msg = "No starnnw name was passed. Only paths explicitly passed as arguments are used. " \
              f"Using a default starnnw file of {filter_name}"
        logger.warning(msg)

    if catalog_name is None:
        image_name = Path(img).stem
        catalog_name = f'{image_name}.cat'

    cmd = f"sex {img} " \
          f"-c {config} " \
          f"-CATALOG_NAME {catalog_name} " \
          f"-PARAMETERS_NAME {parameters_name} " \
          f"-FILTER_NAME {filter_name} " \
          f"-STARNNW_NAME {starnnw_name} " \
          f"-VERBOSE_TYPE {verbose_type} " \
          f"-SATUR_LEVEL {saturation} " \
          f"-GAIN {gain:.3f} "

    cmd += parse_checkimage(
        checkimage_type=checkimage_type,
        checkimage_name=checkimage_name,
        image=img
    )

    if weight_image is None:
        cmd += "-WEIGHT_TYPE None"
    else:
        cmd += f"-WEIGHT_IMAGE {weight_image}"

    try:
        execute(cmd, output_dir, local=run_local)
    except ExecutionError as e:
        raise SextractorError(e)


if __name__ == "__main__":
    run_sextractor(
        "/Users/robertstein/Data/WIRC/20200929/redux/image0240.fits",
        "/Users/robertstein/Data/testersextractor",
        checkimage_type=["BACKGROUND", "BACKGROUND_RMS"]
    )
