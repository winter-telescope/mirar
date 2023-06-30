"""
Module to run source extractor
"""
import logging
import os
from pathlib import Path
from typing import Optional

from mirar.processors.astromatic.config import astromatic_config_dir
from mirar.processors.candidates.utils.regions_writer import write_regions_file
from mirar.utils import ExecutionError, execute
from mirar.utils.ldac_tools import get_table_from_ldac

logger = logging.getLogger(__name__)

# sextractor_cmd = os.getenv("SEXTRACTOR_CMD")

DEFAULT_SATURATION = 10000000000.0
default_config_path = os.path.join(astromatic_config_dir, "astrom.sex")
default_param_path = os.path.join(astromatic_config_dir, "astrom.param")
default_filter_name = os.path.join(astromatic_config_dir, "default.conv")
default_starnnw_path = os.path.join(astromatic_config_dir, "default.nnw")


class SextractorError(ExecutionError):
    """
    Sextractor processing error
    """


# Either run sextractor locally or on docker

LOCAL_SEXTRACTOR = True


# Functions to parse commands and generate appropriate sextractor files


def parse_checkimage(
    checkimage_type: Optional[str | list] = None,
    checkimage_name: Optional[str | list] = None,
    image: Optional[str] = None,
):
    """Function to parse the "checkimage" component of Sextractor configuration.

    Parameters
    ----------
    checkimage_type: The 'CHECKIMAGE_TYPE' files for sextractor. The default is None.
    To quote sextractor,
    available types are: 'NONE, BACKGROUND, BACKGROUND_RMS, MINIBACKGROUND,
    MINIBACK_RMS, -BACKGROUND,
    FILTERED, OBJECTS, -OBJECTS, SEGMENTATION, or APERTURES'. Multiple arguments should
    be specified in a list.
    checkimage_name: The name(s) of the checkput images to be output.
    image: The name of the image in question. If specified, the name of each
    checkimage will include the
    name of the original base image.

    Returns
    -------
    cmd: A string containing the partial sextractor command relating to checkimages.
    The default is an empty string.
    """
    if isinstance(checkimage_type, str):
        checkimage_type = [checkimage_type]

    if isinstance(checkimage_name, str):
        checkimage_name = [checkimage_name]

    if checkimage_type is not None:
        cmd = f"-CHECKIMAGE_TYPE {','.join(checkimage_type)} "

        if checkimage_name is not None:
            if not len(checkimage_type) == len(checkimage_name):
                err = (
                    f"Number of checkimage types {len(checkimage_type)} does not "
                    f"match number of checkimage names {len(checkimage_name)}. "
                    f"These values must be equal. The following types were given: "
                    f"{checkimage_type}, "
                    f"and the following names were given: {checkimage_name}"
                )
                logger.error(err)
                raise ValueError(err)

            cmd += f"-CHECKIMAGE_NAME {','.join(checkimage_name)}"

        else:
            if image is not None:
                base_name = f'{image.split(".fits")[0]}_'
            else:
                base_name = ""

            checkimage_name = [
                f"{base_name}check_{x.lower()}.fits" for x in checkimage_type
            ]
            cmd += " -CHECKIMAGE_NAME " + ",".join(checkimage_name)

        cmd += " "

    else:
        cmd = " -CHECKIMAGE_TYPE NONE "
        checkimage_name = []
    return cmd, checkimage_name


def run_sextractor(images: str | list, output_dir: str, *args, **kwargs):
    """
    Wrapper function to run sextractor
    Args:
        images:
        output_dir:
        *args:
        **kwargs:

    Returns:

    """
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
    catalog_name: Optional[Path] = None,
    config: str = default_config_path,
    parameters_name: str = default_param_path,
    filter_name: str = default_filter_name,
    starnnw_name: str = default_starnnw_path,
    saturation: float = None,
    weight_image: Optional[str] = None,
    verbose_type: str = "QUIET",
    checkimage_name: Optional[str | list] = None,
    checkimage_type: Optional[str | list] = None,
    gain: Optional[float] = None,
    mag_zp: Optional[float] = None,
    write_regions: bool = False,
):
    """
    Function to run sextractor in single mode
    Args:
        img: The image to run sextractor on
        output_dir: The directory to output the catalog to
        catalog_name: The name of the catalog to output.
        config:
        parameters_name:
        filter_name:
        starnnw_name:
        saturation:
        weight_image:
        verbose_type:
        checkimage_name:
        checkimage_type:
        gain: The gain to use for the catalog
        mag_zp: The magnitude zero point to use for the catalog
        write_regions: Whether to write ds9 regions for the objects in the catalog
    Returns:

    """
    if catalog_name is None:
        image_name = Path(img).stem
        catalog_name = f"{image_name}.cat"

    cmd = (
        f"sex {img} "
        f"-c {config} "
        f"-CATALOG_NAME {catalog_name} "
        f"-VERBOSE_TYPE {verbose_type} "
    )

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

    checkimage_cmd, checkimage_name = parse_checkimage(
        checkimage_type=checkimage_type, checkimage_name=checkimage_name, image=img
    )
    cmd += checkimage_cmd

    if weight_image is None:
        cmd += "-WEIGHT_TYPE None"
    else:
        cmd += f"-WEIGHT_IMAGE {weight_image}"

    if mag_zp is not None:
        cmd += f" -MAG_ZEROPOINT {mag_zp}"

    try:
        execute(cmd, output_dir)
    except ExecutionError as exc:
        raise SextractorError(exc) from exc

    if write_regions:
        output_catalog = get_table_from_ldac(catalog_name)

        x_coords = output_catalog["X_IMAGE"]
        y_coords = output_catalog["Y_IMAGE"]

        regions_path = catalog_name.as_posix() + ".reg"

        write_regions_file(
            regions_path=regions_path,
            x_coords=x_coords,
            y_coords=y_coords,
            system="image",
            region_radius=5,
        )
    return catalog_name, checkimage_name


def run_sextractor_dual(
    det_image: str,
    measure_image: str,
    output_dir: str,
    catalog_name: Optional[str] = None,
    config: str = default_config_path,
    parameters_name: str = default_param_path,
    filter_name: str = default_filter_name,
    starnnw_name: str = default_starnnw_path,
    saturation: float = None,
    weight_image: Optional[str] = None,
    verbose_type: str = "QUIET",
    checkimage_name: Optional[str | list] = None,
    checkimage_type: Optional[str | list] = None,
    gain: Optional[float] = None,
    mag_zp: Optional[float] = None,
):
    """
    Run sextractor in the dual mode
    Args:
        det_image:
        measure_image:
        output_dir:
        catalog_name:
        config:
        parameters_name:
        filter_name:
        starnnw_name:
        saturation:
        weight_image:
        verbose_type:
        checkimage_name:
        checkimage_type:
        gain:
        mag_zp:

    Returns:

    """
    if catalog_name is None:
        image_name = Path(measure_image).stem
        catalog_name = f"{image_name}.cat"

    cmd = (
        f"sex {det_image},{measure_image} "
        f"-c {config} "
        f"-CATALOG_NAME {catalog_name} "
        f"-VERBOSE_TYPE {verbose_type} "
    )

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

    if mag_zp is not None:
        cmd += f" -MAG_ZEROPOINT {mag_zp}"

    checkimage_cmd, checkimage_name = parse_checkimage(
        checkimage_type=checkimage_type,
        checkimage_name=checkimage_name,
        image=det_image,
    )

    cmd += checkimage_cmd
    logger.info(weight_image)
    if weight_image is None:
        cmd += "-WEIGHT_TYPE None"
    else:
        cmd += f"-WEIGHT_IMAGE {weight_image}"

    try:
        execute(cmd, output_dir)
    except ExecutionError as err:
        raise SextractorError(err) from err

    return catalog_name, checkimage_name
