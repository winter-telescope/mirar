from winterdrp.processors.dark import DarkCalibrator
from winterdrp.processors.flat import SkyFlatCalibrator
from winterdrp.processors.sky import NightSkyMedianCalibrator
from winterdrp.processors.mask import MaskPixels
from winterdrp.processors.utils import ImageSaver
from winterdrp.pipelines.wirc.wirc_files import wirc_mask_path, sextractor_astrometry_config, scamp_fp_path, \
    swarp_sp_path
from winterdrp.processors.autoastrometry import AutoAstrometry
from winterdrp.processors.astromatic import Sextractor, Scamp, Swarp
from winterdrp.processors.photcal import PhotCalibrator
from winterdrp.processors.utils.image_selector import ImageSelector, ImageBatcher, ImageDebatcher
from winterdrp.processors.utils import ImageLoader
from winterdrp.processors.csvlog import CSVLog
from winterdrp.pipelines.wirc.load_wirc_image import load_raw_wirc_image
from winterdrp.pipelines.wirc.generator import wirc_astrometric_catalog_generator, wirc_photometric_catalog_generator

test_pipeline = [
    ImageLoader(
        input_sub_dir="raw",
        load_image=load_raw_wirc_image
    ),
    CSVLog(
        export_keys=["OBJECT", "FILTER", "UTSHUT", "EXPTIME", "COADDS", "OBSTYPE", "OBSCLASS"]
    ),
    MaskPixels(mask_path=wirc_mask_path),
    ImageSelector(("exptime", "45.0")),
    DarkCalibrator(),
    ImageDebatcher(),
    ImageSelector(("obsclass", "science")),
    ImageBatcher(split_key="filter"),
    SkyFlatCalibrator(),
    NightSkyMedianCalibrator(),
    AutoAstrometry(catalog="tmc"),
    Sextractor(
        output_sub_dir="postprocess",
        **sextractor_astrometry_config
    ),
    Scamp(
        ref_catalog_generator=wirc_astrometric_catalog_generator,
        scamp_config_path=scamp_fp_path,
    ),
    Swarp(swarp_config_path=swarp_sp_path),
    Sextractor(
        output_sub_dir="final_sextractor",
        **sextractor_astrometry_config
    ),
    ImageSaver(output_dir_name="final"),
    PhotCalibrator(ref_catalog_generator=wirc_photometric_catalog_generator)
]