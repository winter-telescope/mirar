import copy
import os

import astropy.io
import numpy as np
import logging
from winterdrp.processors.utils.image_selector import select_from_images
from winterdrp.processors.utils.image_loader import ImageLoader, load_from_dir, BaseImageProcessor
from winterdrp.errors import ImageNotFoundError, ProcessorError
from winterdrp.io import open_fits
from collections.abc import Callable
from winterdrp.paths import raw_img_sub_dir
from pathlib import Path
from winterdrp.data import ImageBatch, Image

logger = logging.getLogger(__name__)


class CalRequirement:

    def __init__(self, target_name, required_field: str, required_values: str | list[str]):
        self.target_name = target_name
        self.required_field = required_field
        self.required_values = required_values
        self.success = False
        self.data = dict()

    def check_images(self, images):
        new_images = select_from_images(
            images,
            key="TARGET",
            target_values=self.target_name
        )

        if len(new_images) > 0:
            for value in self.required_values:
                if value not in self.data.keys():
                    sub_images = select_from_images(
                        new_images,
                        key=self.required_field,
                        target_values=value
                    )
                    if len(sub_images) > 0:
                        self.data[value] = [sub_images]

        self.success = len(self.data) == len(self.required_values)


def update_requirements(
        requirements: list[CalRequirement],
        images: ImageBatch,
) -> list[CalRequirement]:

    for requirement in requirements:
        if not requirement.success:
            requirement.check_images(images)
        logger.debug(f"{requirement}, {requirement.success}")

    return requirements


def find_required_cals(
        latest_dir: str,
        night: str,
        requirements: list[CalRequirement],
        open_f: Callable = open_fits,
        images: ImageBatch = None,
        skip_latest_night: bool = False
) -> ImageBatch:

    path = Path(latest_dir)
    logger.debug(path)
    split = latest_dir.split(night)

    root = split[0]

    if len(split) > 1:
        subdir = split[1][1:]
    else:
        subdir = ""

    preceding_dirs = []

    for x in [x for x in Path(root).iterdir() if x.is_dir()]:
        if x.name[0] not in ["."]:
            if len(str(x.name)) == len(str(night)):
                try:
                    if float(x.name) < float(night):
                        preceding_dirs.append(x)
                except ValueError:
                    pass

    if not skip_latest_night:
        preceding_dirs.append(path.parent)

    ordered_nights = sorted(preceding_dirs)[::-1]

    while np.sum([x.success for x in requirements]) != len(requirements):

        if len(ordered_nights) == 0:
            raise ImageNotFoundError("Ran out of nights!")

        dir_to_load = ordered_nights[0].joinpath(subdir)

        ordered_nights = ordered_nights[1:]

        try:
            logger.info(f"Checking night {dir_to_load}")
            new_images = load_from_dir(
                str(dir_to_load), open_f=open_f
            )
            requirements = update_requirements(requirements, new_images)

        except ImageNotFoundError:
            pass

    n_cal = 0

    for requirement in requirements:
        for key, (cal_imgs, cal_headers) in requirement.data.items():
            for i, cal_header in enumerate(cal_headers):
                if cal_header not in images:
                    images.append(cal_imgs[i])
                    n_cal += 1

    if n_cal > 0:
        logger.warning(f"Some required calibration images were missing from image set. "
                       f"Found {n_cal} additional calibration images from older nights")

    return images


class CalHunter(ImageLoader):

    base_key = "calhunt"

    def __init__(
            self,
            requirements: CalRequirement | list[CalRequirement],
            *args,
            **kwargs
    ):
        super().__init__(*args, **kwargs)

        if not isinstance(requirements, list):
            requirements = [requirements]

        self.requirements = requirements

    def __str__(self):
        reqs = [f"{x.target_name.upper()} images" for x in self.requirements]
        return f"Processor to search through archival data to find any missing {' and '.join(reqs)}"

    def _apply_to_images(
            self,
            batch: ImageBatch,
    ) -> ImageBatch:

        requirements = copy.deepcopy(self.requirements)
        requirements = update_requirements(requirements, batch)

        latest_dir = os.path.join(
            self.input_img_dir,
            os.path.join(self.night_sub_dir, self.input_sub_dir)
        )

        images = find_required_cals(
            latest_dir=latest_dir,
            night=self.night,
            requirements=requirements,
            open_f=self.load_image,
            images=batch,
            skip_latest_night=True
        )
        batch += images

        return batch


