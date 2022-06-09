from winterdrp.pipelines.base_pipeline import Pipeline
from winterdrp.references.wirc import WIRCRef

def wirc_reference_image_generator():
    return WIRCRef()


class WircImsubPipeline(Pipeline):

    '''
    pipeline_configurations = {
        None: [
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

