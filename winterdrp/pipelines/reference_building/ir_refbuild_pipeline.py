import os
from typing import Callable

import pandas as pd

from winterdrp.catalog import Gaia2Mass
from winterdrp.data import Image, ImageBatch
from winterdrp.pipelines.base_pipeline import Pipeline
from winterdrp.pipelines.wirc.load_wirc_image import load_raw_wirc_image
from winterdrp.pipelines.wirc.wirc_files import sextractor_astrometry_config
from winterdrp.processors.astromatic.sextractor.sextractor import Sextractor
from winterdrp.processors.astromatic.swarp.swarp import Swarp
from winterdrp.processors.base_processor import BaseImageProcessor
from winterdrp.processors.photcal import PhotCalibrator
from winterdrp.processors.utils import ImageLoader, ImageSaver
from winterdrp.references import BaseReferenceGenerator
from winterdrp.references.ukirt import UKIRTRef

refbuild_dir = os.path.dirname(__file__)

astromatic_config_dir = os.path.join(refbuild_dir, "config/files")
swarp_config_path = os.path.join(astromatic_config_dir, "swarp.config")


def winter_reference_generator(image: Image):
    filtername = image["FILTER"]
    # TODO if in_ukirt and in_vista, different processing
    return UKIRTRef(
        filter_name=filtername,
        swarp_resampler=winter_reference_image_resampler,
        phot_calibrator=winter_reference_phot_calibrator,
        num_query_points=4,
    )


def dummy_image_generator():
    winter_fields = pd.read_csv(
        "~/winter/gwemopt_sims/input/WINTER_fields.txt", delim_whitespace=True
    )


def winter_reference_image_resampler(**kwargs) -> Swarp:
    """
    Generates a resampler for reference images

    :param kwargs: kwargs
    :return: Swarp processor
    """
    return Swarp(swarp_config_path=swarp_config_path, cache=True, **kwargs)


def wirc_photometric_catalog_generator(image: Image) -> Gaia2Mass:
    """
    Function to crossmatch WIRC to GAIA/2mass for photometry

    :param image: Image
    :return: catalogue
    """
    filter_name = image["FILTER"]
    return Gaia2Mass(
        min_mag=10, max_mag=20, search_radius_arcmin=30, filter_name=filter_name
    )


def winter_reference_phot_calibrator(image: Image, **kwargs) -> PhotCalibrator:
    """
    Generates a resampler for reference images

    :param kwargs: kwargs
    :return: Swarp processor
    """
    return PhotCalibrator(
        ref_catalog_generator=wirc_photometric_catalog_generator, **kwargs
    )


class MakeDummyImages(BaseImageProcessor):
    pass


class GetReferenceImage(BaseImageProcessor):
    base_key = "refimg_returner"

    def __init__(
        self,
        ref_image_generator: Callable[..., BaseReferenceGenerator],
    ):
        super().__init__()
        self.ref_image_generator = ref_image_generator

    def _apply_to_images(
        self,
        batch: ImageBatch,
    ) -> ImageBatch:
        ref_batch = ImageBatch()
        for image in batch:
            ref_generator = self.ref_image_generator(image)
            ref_image_hdu = ref_generator.get_reference(image)
            ref_image = Image(data=ref_image_hdu.data, header=ref_image_hdu.header)
            ref_batch.append(ref_image)

        return ref_batch


class IRRefBuildPipeline(Pipeline):
    name = "ir_refbuild"

    refbuild = [
        ImageLoader(load_image=load_raw_wirc_image),
        GetReferenceImage(
            ref_image_generator=winter_reference_generator,
        ),
        Sextractor(
            output_sub_dir="phot",
            **sextractor_astrometry_config,
            write_regions_bool=True
        ),
        PhotCalibrator(ref_catalog_generator=wirc_photometric_catalog_generator),
        ImageSaver(output_dir_name="stacked_ref"),
    ]

    all_pipeline_configurations = {"default": refbuild}
