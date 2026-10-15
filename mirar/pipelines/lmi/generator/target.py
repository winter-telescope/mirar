"""
Module to annotate target coordinates on images.
"""

import logging

import requests

from mirar.data import ImageBatch
from mirar.paths import OBSCLASS_KEY, TARGET_KEY, TIME_KEY
from mirar.processors.skyportal.client import SkyportalClient

logger = logging.getLogger(__name__)

# One lookup per name is enough, however many dithers and stacks use it
resolved_positions: dict[str, tuple[float, float] | None] = {}


def query_skyportal_position(name: str) -> tuple[float, float] | None:
    """
    Query SkyPortal for the catalogue position of one named source

    A name which cannot be resolved is not an error, since plenty of targets are
    never posted to SkyPortal, so this warns and leaves the caller to fall back.

    :param name: Name of the target, e.g. 'ZTF26abzrlze'
    :return: Right ascension and declination in degrees, or None if not resolved
    """
    try:
        res = SkyportalClient().api("GET", f"sources/{name}")
    except (requests.RequestException, ValueError) as exc:
        logger.warning(f"Could not reach SkyPortal to resolve '{name}': {exc}")
        return None

    if res.status_code != 200:
        logger.warning(
            f"SkyPortal returned no source named '{name}' (status {res.status_code})"
        )
        return None

    try:
        data = res.json()["data"]
        return float(data["ra"]), float(data["dec"])
    except (ValueError, KeyError, TypeError) as exc:
        logger.warning(f"Could not read a position for '{name}' from SkyPortal: {exc}")
        return None


def resolve_target_coordinates(name: str) -> tuple[float, float] | None:
    """
    Get the catalogue position of a target by name, caching the lookup

    :param name: Name of the target, e.g. 'ZTF26abzrlze'
    :return: Right ascension and declination in degrees, or None if not resolved
    """
    if name not in resolved_positions:
        resolved_positions[name] = query_skyportal_position(name)

    return resolved_positions[name]


def annotate_target_coordinates(image_batch: ImageBatch) -> ImageBatch:
    """
    Function to annotate target coordinates on images.

    LMI headers record the telescope pointing, which is dithered and can sit an
    arcminute from the target, so the position is resolved from the target name
    instead. Only the pointing is available if that fails, which is enough to
    group and stack the images, but not to place a forced photometry aperture.

    :param image_batch: ImageBatch object
    :return: ImageBatch object
    """

    times = [image[TIME_KEY] for image in image_batch]
    min_time = min(times)
    first_image = image_batch[times.index(min_time)]

    # In case one of the dithers is mis-named, we'll use the most common name
    names = [x[TARGET_KEY] for x in image_batch]
    most_common_name = max(set(names), key=names.count)

    position = None

    if first_image[OBSCLASS_KEY] == "science":
        position = resolve_target_coordinates(most_common_name)

        if position is None:
            logger.warning(
                f"Falling back to the telescope pointing for '{most_common_name}', "
                f"which is dithered and may not sit on the target"
            )

    if position is None:
        position = (first_image["OBJRA"], first_image["OBJDEC"])

    for image in image_batch:
        image["OBJRA"] = position[0]
        image["OBJDEC"] = position[1]
        image[TARGET_KEY] = most_common_name

    return image_batch
