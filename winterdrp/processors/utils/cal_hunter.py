import copy
import os

import astropy.io
import numpy as np
import logging
from winterdrp.processors.utils.image_selector import select_from_images
from winterdrp.processors.utils.image_loader import ImageLoader, load_from_dir, BaseImageProcessor
from winterdrp.errors import ImageNotFoundError

logger = logging.getLogger(__name__)


class CalRequirement:

    def __init__(self, target_name, required_field: str, required_values: str | list[str]):
        self.target_name = target_name
        self.required_field = required_field
        self.required_values = required_values
        self.success = False
        self.data = dict()

    def check_images(self, images, headers):

        new_images, new_headers = select_from_images(
            images, headers,
            header_key="TARGET",
            target_values=self.target_name
        )

        if len(new_images) > 0:
            for value in self.required_values:
                if value not in self.data.keys():
                    sub_images, sub_headers = select_from_images(
                        new_images, new_headers,
                        header_key=self.required_field,
                        target_values=value
                    )
                    if len(sub_images) > 0:
                        self.data[value] = [sub_images, sub_headers]

        self.success = len(self.data) == len(self.required_values)


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

    def _apply_to_images(
            self,
            images: list[np.ndarray],
            headers: list[astropy.io.fits.Header],
    ) -> tuple[list[np.ndarray], list[astropy.io.fits.Header]]:

        requirements = copy.deepcopy(self.requirements)
        requirements = self.update_requirements(requirements, images, headers)

        latest_dir = os.path.join(
            self.input_img_dir,
            os.path.join(self.night_sub_dir, self.input_sub_dir)
        )

        preceding_dirs = []

        for x in os.listdir(os.path.dirname(os.path.dirname(latest_dir))):
            if x[0] not in ["."]:
                if len(str(x)) == len(str(self.night)):
                    try:
                        if float(x) < float(self.night):
                            preceding_dirs.append(x)
                    except ValueError:
                        pass

        ordered_nights = sorted(preceding_dirs)[::-1]

        while np.sum([x.success for x in requirements]) != len(requirements):

            if len(ordered_nights) == 0:
                raise ImageNotFoundError("Ran out of nights!")

            new_latest_night = ordered_nights[0]
            ordered_nights = ordered_nights[1:]

            try:

                logger.info(f"Checking night {new_latest_night}")

                dir_to_load = latest_dir.replace(self.night, new_latest_night)

                new_images, new_headers = load_from_dir(
                    dir_to_load, open_f=self.open_raw_image
                )

                requirements = self.update_requirements(requirements, new_images, new_headers)

            except ImageNotFoundError:
                pass

        n_cal = 0

        for requirement in requirements:
            for key, (cal_imgs, cal_headers) in requirement.data.items():
                for i, cal_header in enumerate(cal_headers):
                    if cal_header not in headers:
                        images.append(cal_imgs[i])
                        headers.append(cal_header)
                        n_cal += 1

        if n_cal > 0:
            logger.warning(f"Some required calibration images were missing from image set. "
                           f"Found {n_cal} additional calibration images from older nights")

        return images, headers

    @staticmethod
    def update_requirements(
            requirements: list[CalRequirement],
            images: list[np.ndarray],
            headers: list[astropy.io.fits.Header]
    ) -> list[CalRequirement]:

        for requirement in requirements:
            if not requirement.success:
                requirement.check_images(images, headers)

        return requirements
