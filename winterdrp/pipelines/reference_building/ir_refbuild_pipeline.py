import logging
import os
from typing import Callable

import pandas as pd

from winterdrp.catalog import Gaia2Mass
from winterdrp.data import Image, ImageBatch
from winterdrp.paths import get_output_dir
from winterdrp.pipelines.base_pipeline import Pipeline
from winterdrp.pipelines.reference_building.db_models import RefComponents, RefStacks
from winterdrp.pipelines.wirc.load_wirc_image import load_raw_wirc_image
from winterdrp.pipelines.wirc.wirc_files import sextractor_astrometry_config
from winterdrp.processors.astromatic.sextractor.sextractor import Sextractor
from winterdrp.processors.astromatic.swarp.swarp import Swarp
from winterdrp.processors.base_processor import BaseImageProcessor
from winterdrp.processors.photcal import PhotCalibrator
from winterdrp.processors.sqldatabase.database_exporter import DatabaseImageExporter
from winterdrp.processors.utils import ImageDebatcher, ImageLoader, ImageSaver
from winterdrp.references import BaseReferenceGenerator
from winterdrp.references.ukirt import UKIRTRef

logger = logging.getLogger(__name__)

refbuild_dir = os.path.dirname(__file__)

astromatic_config_dir = os.path.join(refbuild_dir, "config/")
swarp_config_path = os.path.join(astromatic_config_dir, "config.swarp")


def winter_reference_generator(image: Image):
    filtername = image["FILTER"]
    # TODO check if exists in DB
    # TODO if in_ukirt and in_vista, different processing
    return UKIRTRef(
        filter_name=filtername,
        swarp_resampler=winter_reference_image_resampler,
        sextractor_generator=ref_sextractor,
        phot_calibrator_generator=winter_reference_phot_calibrator,
        num_query_points=4,
        components_table=RefComponents,
        write_to_db=True,
        write_db_table=RefStacks,
    )


def winter_reference_image_resampler(**kwargs) -> Swarp:
    """
    Generates a resampler for reference images

    :param kwargs: kwargs
    :return: Swarp processor
    """
    logger.info(kwargs)
    return Swarp(
        swarp_config_path=swarp_config_path, subtract_bkg=True, cache=True, **kwargs
    )


def wirc_photometric_catalog_generator(image: Image) -> Gaia2Mass:
    """
    Function to crossmatch WIRC to GAIA/2mass for photometry

    :param image: Image
    :return: catalogue
    """
    filter_name = image["FILTER"]
    return Gaia2Mass(
        min_mag=10,
        max_mag=20,
        search_radius_arcmin=10,
        filter_name=filter_name,
        snr_threshold=20,
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


def ref_sextractor(image: Image):
    return Sextractor(
        output_sub_dir="phot",
        **sextractor_astrometry_config,
        write_regions_bool=True,
        cache=True
    )


def ref_phot_calibrator(image: Image):
    return PhotCalibrator(
        ref_catalog_generator=wirc_photometric_catalog_generator,
        write_regions=True,
        x_lower_limit=0,
        fwhm_threshold_arcsec=3,
    )


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
            ref_image_path = ref_generator.write_reference(
                image,
                output_dir=get_output_dir(
                    dir_root="mock/stacked_ref", sub_dir="ir_refbuild"
                ).as_posix(),
            )

            ref_image = self.open_fits(ref_image_path)

            ref_batch.append(ref_image)

        return ref_batch


class IRRefBuildPipeline(Pipeline):
    name = "ir_refbuild"

    refbuild = [
        # ImageLoader(load_image=load_raw_wirc_image),
        ImageDebatcher(),
        GetReferenceImage(
            ref_image_generator=winter_reference_generator,
        ),
        # Sextractor(
        #     output_sub_dir="phot",
        #     **sextractor_astrometry_config,
        #     write_regions_bool=True,
        #     cache=True
        # ),
        # PhotCalibrator(
        #     ref_catalog_generator=wirc_photometric_catalog_generator,
        #     write_regions=True,
        #     x_lower_limit=0,
        #     fwhm_threshold_arcsec=3,
        # ),
        ImageSaver(output_dir_name="stacked_ref"),
        # DatabaseImageExporter(
        #     db_table=RefStacks,
        #     duplicate_protocol="replace",
        #     q3c_bool=False,
        # ),
    ]

    all_pipeline_configurations = {"default": refbuild}
