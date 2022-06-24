import os
import logging
from pathlib import Path
from winterdrp.processors.astromatic.config import astromatic_config_dir
from winterdrp.utils import execute, ExecutionError

logger = logging.getLogger(__name__)

# sextractor_cmd = os.getenv("SEXTRACTOR_CMD")

default_saturation = 1.e10
default_config_path = os.path.join(astromatic_config_dir, 'astrom.sex')
default_param_path = os.path.join(astromatic_config_dir, 'astrom.param')
default_filter_name = os.path.join(astromatic_config_dir, 'default.conv')
default_starnnw_path = os.path.join(astromatic_config_dir, 'default.nnw')


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
        checkimage_type = [checkimage_type]

    if isinstance(checkimage_name, str):
        checkimage_name = [checkimage_name]

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

    else:

        cmd = f" -CHECKIMAGE_TYPE NONE "
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
        run_sextractor_single(img, output_dir, *args, **kwargs)


def run_sextractor_single(
        img: str,
        output_dir: str,
        catalog_name: str = None,
        config: str = default_config_path,
        parameters_name: str = default_param_path,
        filter_name: str = default_filter_name,
        starnnw_name: str = default_starnnw_path,
        saturation: float = default_saturation,
        weight_image: str = None,
        verbose_type: str = "QUIET",
        checkimage_name: str | list = None,
        checkimage_type: str | list = None,
        gain: float = None,
):
    if catalog_name is None:
        image_name = Path(img).stem
        catalog_name = f'{image_name}.cat'

    cmd = f"sex {img} " \
          f"-c {config} " \
          f"-CATALOG_NAME {catalog_name} " \
          f"-VERBOSE_TYPE {verbose_type} "

    if saturation is not None:
        cmd += f"-SATUR_LEVEL {saturation} "

    if gain is not None:
        cmd += f"-GAIN {gain:.3f} "

    if parameters_name is not None:
        cmd += f"-PARAMETERS_NAME {parameters_name} "

    if filter_name is not None:
        cmd += f"-FILTER_NAME {filter_name} "

    if starnnw_name is not None:
        cmd += f"-STARNNW_NAME {starnnw_name} "

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
        execute(cmd, output_dir)
    except ExecutionError as e:
        raise SextractorError(e)

    return catalog_name


def run_sextractor_dual(
        det_image: str,
        measure_image: str,
        output_dir: str,
        catalog_name: str = None,
        config: str = default_config_path,
        parameters_name: str = default_param_path,
        filter_name: str = default_filter_name,
        starnnw_name: str = default_starnnw_path,
        saturation: float = default_saturation,
        weight_image: str = None,
        verbose_type: str = "QUIET",
        checkimage_name: str | list = None,
        checkimage_type: str | list = None,
        gain: float = None,

):
    if catalog_name is None:
        image_name = Path(measure_image).stem
        catalog_name = f'{image_name}.cat'

    cmd = f"sex {det_image},{measure_image} " \
          f"-c {config} " \
          f"-CATALOG_NAME {catalog_name} " \
          f"-VERBOSE_TYPE {verbose_type} "

    if saturation is not None:
        cmd += f"-SATUR_LEVEL {saturation} "

    if gain is not None:
        cmd += f"-GAIN {gain:.3f} "

    if parameters_name is not None:
        cmd += f"-PARAMETERS_NAME {parameters_name} "

    if filter_name is not None:
        cmd += f"-FILTER_NAME {filter_name} "

    if starnnw_name is not None:
        cmd += f"-STARNNW_NAME {starnnw_name} "

    cmd += parse_checkimage(
        checkimage_type=checkimage_type,
        checkimage_name=checkimage_name,
        image=det_image
    )

    if weight_image is None:
        cmd += "-WEIGHT_TYPE None"
    else:
        cmd += f"-WEIGHT_IMAGE {weight_image}"

    try:
        execute(cmd, output_dir)
    except ExecutionError as e:
        raise SextractorError(e)

    return catalog_name
