from winterdrp.pipelines.base_pipeline import Pipeline
from winterdrp.references.wirc import WIRCRef
from winterdrp.processors.astromatic.swarp.swarp import Swarp
from winterdrp.processors.astromatic.sextractor.sextractor import Sextractor
from winterdrp.processors.astromatic.psfex import PSFex
from winterdrp.processors.reference import Reference
from winterdrp.processors.zogy.zogy import ZOGY, ZOGYPrepare
from winterdrp.processors.candidates import DetectCandidates
import numpy as np
from astropy.io import fits
import os
from astropy.time import Time
import logging

logger = logging.getLogger(__name__)


def wirc_reference_image_generator(
        header: fits.header,
        images_directory: str = "/Users/viraj/wirc_imsub",
):
    object_name = header['OBJECT']
    filter_name = header['FILTER']
    return WIRCRef(
        object_name=object_name,
        filter_name=filter_name,
        images_directory_path=images_directory
    )


def wirc_reference_image_resampler(pixsize,
                                   x_imgpixsize,
                                   y_imgpixsize,
                                   center_ra,
                                   center_dec,
                                   propogate_headerlist,
                                   temp_output_sub_dir,
                                   night_sub_dir,
                                   include_scamp,
                                   combine):
    logger.debug(f'Night sub dir is {night_sub_dir}')
    return Swarp(swarp_config_path='~/wirc_imsub/config/config.swarp',
                 pixsize=pixsize,
                 x_imgpixsize=x_imgpixsize,
                 y_imgpixsize=y_imgpixsize,
                 center_ra=center_ra,
                 center_dec=center_dec,
                 propogate_headerlist=propogate_headerlist,
                 temp_output_sub_dir=temp_output_sub_dir,
                 night_sub_dir=night_sub_dir,
                 include_scamp=include_scamp,
                 combine=combine
                 )


def wirc_reference_sextractor(output_sub_dir, gain):
    return Sextractor(config_path='~/wirc_imsub/config/photomCat.sex',
                      parameter_path='~/wirc_imsub/config/photom.param',
                      filter_path='~/wirc_imsub/config/default.conv',
                      starnnw_path='~/wirc_imsub/config/default.nnw',
                      gain=gain,
                      output_sub_dir=output_sub_dir
                      )


def wirc_reference_psfex(output_sub_dir, norm_fits):
    return PSFex(config_path='~/wirc_imsub/config/photom.psfex',
                 output_sub_dir=output_sub_dir,
                 norm_fits=norm_fits
                 )


def detect_candidates_sextractor():
    pass


class WircImsubPipeline(Pipeline):
    name = "wirc_imsub"

    header_keys = [
        "UTSHUT",
        'OBJECT',
        "FILTER",
        "EXPTIME"
    ]
    batch_split_keys = ["UTSHUT"]

    pipeline_configurations = {
        None: [
            Reference(
                ref_image_generator=wirc_reference_image_generator,
                ref_swarp_resampler=wirc_reference_image_resampler,
                ref_sextractor=wirc_reference_sextractor,
                ref_psfex=wirc_reference_psfex
            ),
            # Swarp(),
            Sextractor(config_path='~/wirc_imsub/config/photomCat.sex',
                       parameter_path='~/wirc_imsub/config/photom.param',
                       filter_path='~/wirc_imsub/config/default.conv',
                       starnnw_path='~/wirc_imsub/config/default.nnw',
                       output_sub_dir='subtract'),
            PSFex(config_path='~/wirc_imsub/config/photom.psfex',
                  output_sub_dir="subtract",
                  norm_fits=True),
            ZOGYPrepare(output_sub_dir="subtract"),
            ZOGY(output_sub_dir="subtract"),
            DetectCandidates(output_sub_dir="subtract",
                             cand_det_sextractor_config='~/wirc_imsub/config/photomCat.sex',
                             cand_det_sextractor_nnw='~/wirc_imsub/config/default.nnw',
                             cand_det_sextractor_filter='~/wirc_imsub/config/default.conv',
                             cand_det_sextractor_params='~/wirc_imsub/config/Scorr.param')
        ]
    }

    def load_raw_image(
            self,
            path: str
    ) -> tuple[np.array, fits.Header]:
        with fits.open(path) as img:
            data = img[0].data
            header = img[0].header
            header["FILTER"] = header["AFT"].split("__")[0]
            header["OBSCLASS"] = ["calibration", "science"][header["OBSTYPE"] == "object"]
            header["CALSTEPS"] = ""
            header["BASENAME"] = os.path.basename(path)
            header["TARGET"] = header["OBJECT"].lower()
            header["UTCTIME"] = header["UTSHUT"]
            header["MJD-OBS"] = Time(header['UTSHUT']).mjd
            # header.append(('GAIN', self.gain, 'Gain in electrons / ADU'), end=True)
            # header = self.set_saturation(header)
            if not 'COADDS' in header.keys():
                logger.debug('Setting COADDS to 0')
                header['COADDS'] = 0
            if not 'CALSTEPS' in header.keys():
                logger.debug('Setting CALSTEPS to blank')
                header['CALSTEPS'] = ''
        return data, header

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
