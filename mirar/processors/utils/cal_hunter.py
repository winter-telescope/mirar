"""
Module for finding archival calibration images if these are missing
"""
import copy
import logging
import os
from collections.abc import Callable
from pathlib import Path

import numpy as np

from mirar.data import Image, ImageBatch
from mirar.errors import ImageNotFoundError
from mirar.io import open_fits
from mirar.processors.utils.image_loader import ImageLoader, load_from_dir
from mirar.processors.utils.image_selector import select_from_images

logger = logging.getLogger(__name__)


class CalRequirement:
    """
    Class to specify particular calibration files that must be present for processing
    """

    target_name: str
    required_field: str
    required_values: str | list[str]

    def __init__(
        self, target_name, required_field: str, required_values: str | list[str]
    ):
        self.target_name = target_name
        self.required_field = required_field
        self.required_values = required_values
        self.success = False
        self.data = {}

    def check_images(self, images: ImageBatch):
        """
        Check a batch of images, to see whether the calibration requirement is met.
        Adds any required images to the cache, then updates
        the check of the self.success

        :param images: ImageBatch
        :return: None
        """

        new_images = select_from_images(
            images, key="TARGET", target_values=self.target_name
        )

        if len(new_images) > 0:
            for value in self.required_values:
                if value not in self.data:
                    sub_images = select_from_images(
                        new_images, key=self.required_field, target_values=value
                    )
                    if len(sub_images) > 0:
                        self.data[value] = sub_images

        self.success = len(self.data) == len(self.required_values)

    def __str__(self):
        return (
            f"<Calibration requirement, checking for '{self.target_name}' images "
            f"with '{self.required_field}'values {self.required_values} >"
        )


def update_requirements(
    requirements: list[CalRequirement],
    images: ImageBatch,
) -> list[CalRequirement]:
    """
    Iteratively check a list of Cal Requirements against an image batch

    :param requirements: CalRequirements to check
    :param images: Images to check
    :return: Updated CalRequirements
    """

    for requirement in requirements:
        if not requirement.success:
            requirement.check_images(images)
        logger.debug(f"{requirement}, {requirement.success}")

    return requirements


def find_required_cals(
    latest_dir: str | Path,
    night: str,
    requirements: list[CalRequirement],
    open_f: Callable[[str], Image] = open_fits,
    images: ImageBatch = ImageBatch(),
    skip_latest_night: bool = False,
) -> ImageBatch:
    """
    Broad function to search for missing calibration files in previous nights

    :param latest_dir: The directory for the raw images
    :param night: The night being processed
    :param requirements: List of calibration requirements
    :param open_f: Function to open raw images
    :param images: Current image list (default: empty)
    :param skip_latest_night: Boolean to skip the directory of night being processed
    :return: Updated image batch
    """

    path = Path(latest_dir)

    logger.debug(f"Searching for archival images for {path}")

    split = latest_dir.split(night)

    root = split[0]

    if len(split) > 1:
        subdir = split[1][1:]
    else:
        subdir = ""

    preceding_dirs = []

    for dir_name in [x for x in Path(root).iterdir() if x.is_dir()]:
        if dir_name.name[0] not in ["."]:
            if len(str(dir_name.name)) == len(str(night)):
                try:
                    if float(dir_name.name) < float(night):
                        preceding_dirs.append(dir_name)
                except ValueError:
                    pass

    if not skip_latest_night:
        preceding_dirs.append(path.parent)

    ordered_nights = sorted(preceding_dirs)[::-1]

    while np.sum([req.success for req in requirements]) != len(requirements):
        if len(ordered_nights) == 0:
            err = (
                "Despite checking all past nights, there are still "
                "missing cal requirements: "
            )

            for req in requirements:
                if not req.success:
                    err += str(req)
            logger.error(err)
            raise ImageNotFoundError(err)

        dir_to_load = ordered_nights[0].joinpath(subdir)

        logger.info(f"Checking night {dir_to_load}")

        ordered_nights = ordered_nights[1:]

        try:
            new_images = load_from_dir(str(dir_to_load), open_f=open_f)
            requirements = update_requirements(requirements, new_images)

        except ImageNotFoundError:
            pass

    n_cal = 0

    for requirement in requirements:
        for cal_imgs in requirement.data.values():
            for cal_img in cal_imgs:
                if cal_img not in images.get_batch():
                    images.append(cal_img)
                    n_cal += 1

    if n_cal > 0:
        logger.warning(
            f"Some required calibration images were missing from image set. "
            f"Found {n_cal} additional calibration images from older nights"
        )

    return images


class CalHunter(ImageLoader):
    """
    Processor to find any missing calibration images,
    by searching previous nights of data
    """

    base_key = "calhunt"

    def __init__(
        self, requirements: CalRequirement | list[CalRequirement], *args, **kwargs
    ):
        super().__init__(*args, **kwargs)

        if not isinstance(requirements, list):
            requirements = [requirements]

        self.requirements = requirements

    def __str__(self):
        reqs = [f"{req.target_name.upper()} images" for req in self.requirements]
        return (
            f"Processor to search through archival data to find any missing "
            f"{' and '.join(reqs)}"
        )

    def _apply_to_images(
        self,
        batch: ImageBatch,
    ) -> ImageBatch:
        requirements = copy.deepcopy(self.requirements)
        requirements = update_requirements(requirements, batch)

        latest_dir = os.path.join(
            self.input_img_dir, os.path.join(self.night_sub_dir, self.input_sub_dir)
        )

        updated_batch = find_required_cals(
            latest_dir=latest_dir,
            night=self.night,
            requirements=requirements,
            open_f=self.open_raw_image,
            images=batch,
            skip_latest_night=True,
        )

        return updated_batch
