from winterdrp.pipelines.base_pipeline import Pipeline
from winterdrp.references.wirc import WIRCRef
from winterdrp.processors.astromatic.swarp.swarp import Swarp
from winterdrp.processors.astromatic.sextractor.sextractor import Sextractor
from winterdrp.processors.astromatic.psfex import PSFex
from winterdrp.processors.reference import Reference


def wirc_reference_image_generator():
    return WIRCRef()


def wirc_reference_image_resampler():
    return Swarp(swarp_config_path='xx')


def wirc_reference_sextractor():
    return Sextractor(config_path='xx',
                      parameter_path='xx',
                      filter_path='xx',
                      starnnw_path='xx'
                      )


def wirc_reference_psfex():
    return Sextractor(config_path='xx',
                      parameter_path='xx',
                      filter_path='xx',
                      starnnw_path='xx'
                      )


class WircImsubPipeline(Pipeline):

    pipeline_configurations = {
        None: [
            Reference(),
            Swarp(),
            Sextractor(),
            PSFex(),
            ImPrepare(),
            Zogy(),
            CandidateDetect()

            ]
    }
    '''
            MaskPixels(mask_path=wirc_mask_path),
            DarkCalibrator(),
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
            PhotCalibrator(ref_catalog_generator=wirc_photometric_catalog_generator),
        ]
    }
    '''
    pass
