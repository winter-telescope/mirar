import logging
import os

import astropy.io.fits
import numpy as np
import copy
from winterdrp.paths import saturate_key


logger = logging.getLogger(__name__)


core_fields = ["OBSCLASS", "TARGET", "UTCTIME"]


class Pipeline:

    pipelines = {}
    name = None

    @property
    def pipeline_configurations(self):
        raise NotImplementedError()

    @property
    def gain(self):
        raise NotImplementedError()

    @property
    def non_linear_level(self):
        raise NotImplementedError()

    def __init__(
            self,
            pipeline_configuration: str = None,
            night: int | str = "",
    ):

        self.night_sub_dir = os.path.join(self.name, night)

        self.processors = self.load_pipeline_configuration(pipeline_configuration)

        self.configure_processors(sub_dir=self.night_sub_dir)

        for i, (processor) in enumerate(self.processors):

            logger.debug(f"Initialising processor {processor.__class__}")
            processor.set_preceding_steps(previous_steps=self.processors[:i])
            processor.check_prerequisites()

        logger.debug("Pipeline initialisation complete.")

    @classmethod
    def __init_subclass__(cls, **kwargs):
        super().__init_subclass__(**kwargs)
        if cls.name in cls.pipelines.keys():
            err = f"Pipeline name '{cls.name}' is already found in the pipeline registered keys. " \
                  f"The Pipeline class variable 'name' must be unique!"
            logger.error(err)
            raise ValueError(err)
        cls.pipelines[cls.name] = cls

    def configure_processors(
            self,
            sub_dir: str = ""
    ):
        for processor in self.processors:
            processor.set_night(night_sub_dir=sub_dir)

    def open_raw_image(
            self,
            path: str
    ) -> tuple[np.array, astropy.io.fits.Header]:

        data, header = self.load_raw_image(path)

        for key in core_fields:
            if key not in header.keys():
                err = f"Essential key {key} not found in header. " \
                      f"Please add this field first. Available fields are: {list(header.keys())}"
                logger.error(err)
                raise KeyError(err)

        return data.astype(np.float64), header

    @staticmethod
    def load_raw_image(
            path: str
    ) -> tuple[np.array, astropy.io.fits.Header]:
        raise NotImplementedError

    @staticmethod
    def download_raw_images_for_night(
            night: str | int
    ):
        raise NotImplemented

    def load_pipeline_configuration(
            self,
            configuration: str | list = None,
    ):
        if isinstance(configuration, str | None):
            return copy.copy(self.pipeline_configurations[configuration])
        else:
            return copy.copy(configuration)

    def reduce_images(
            self,
            batches: list[list[list[np.ndarray], list[astropy.io.fits.header]]],
    ):

        for i, processor in enumerate(self.processors):
            logger.debug(f"Applying '{processor.__class__}' processor to {len(batches)} batches. "
                         f"(Step {i+1}/{len(self.processors)})")
            batches, failures = processor.apply(
                batches
            )

        return batches

    def set_saturation(
            self,
            header: astropy.io.fits.Header
    ) -> astropy.io.fits.Header:
        # update the SATURATE keyword in the header for subsequent sextractor runs
        co_add_head = header['COADDS']
        num_co_adds = int(co_add_head)
        saturation_level = self.non_linear_level * num_co_adds
        if "SKMEDSUB" in header.keys():
            saturation_level -= header['SKMEDSUB']
        header.append((saturate_key, saturation_level, 'Saturation level'), end=True)
        return header
