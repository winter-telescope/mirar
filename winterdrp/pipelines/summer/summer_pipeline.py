import os
import astropy.io.fits
import numpy as np
from astropy.coordinates import SkyCoord
import astropy.units as u
from astropy.io.fits import HDUList
from winterdrp.pipelines.base_pipeline import Pipeline

from winterdrp.processors.bias import BiasCalibrator
from winterdrp.processors.flat import FlatCalibrator
from winterdrp.processors.mask import MaskPixels
from winterdrp.processors.utils import ImageSaver
from winterdrp.processors.autoastrometry import AutoAstrometry
from winterdrp.processors.astromatic import Sextractor, Scamp, Swarp
from winterdrp.catalog import Gaia2Mass
from winterdrp.pipelines.summer.summer_files import summer_mask_path, summer_weight_path, sextractor_astrometry_config, scamp_path, \
    swarp_path
from winterdrp.pipelines.summer.calibration import select_bias, select_flats_archival

summer_flats_dir = os.path.join(os.path.dirname(os.path.abspath(__file__)))
summer_gain = 1.0

def summer_astrometric_catalog_generator(
        header: astropy.io.fits.Header
):
    return Gaia2Mass(min_mag=10, max_mag=20, search_radius_arcmin=30)


class SummerPipeline(Pipeline):

    name = "summer"

    astrometry_cal = ("GAIA", 10., 20.)
    photometry_cal = {
        "J": ()
    }

    # Set up elements to use

    header_keys = [
        "UTC",
        'FIELDID',
        "FILTERID",
        "EXPTIME",
        "OBSTYPE"
    ]

    batch_split_keys = ["RAWIMAGEPATH"]

    #pipeline_configurations = {
    #    None: [
    #        (BiasCalibrator, select_bias),
    #        (FlatCalibrator, select_flats_archival),
    #        (ImageSaver, "preprocess"),
    #        (Sextractor, "pass1"),
    #        # "stack",
    #        # "dither"
    #    ]
    #}

    pipeline_configurations = {
        None: [
            MaskPixels(mask_path=summer_mask_path),
            BiasCalibrator(),
            FlatCalibrator(),
            #ImageSaver(output_dir_name="testa"),
            AutoAstrometry(pa=0, inv=True, pixel_scale=0.466),
            ImageSaver(output_dir_name="testb"),
            Sextractor(
                 output_sub_dir="postprocess",
                 weight_image=summer_weight_path,
                checkimage_name=None,
                checkimage_type=None,
                 **sextractor_astrometry_config
             ),
            #ImageSaver(output_dir_name="testc"),
            Scamp(
                 ref_catalog_generator=summer_astrometric_catalog_generator,
                 scamp_config_path=scamp_path,
             ),
            #ImageSaver(output_dir_name="testd"),
            Swarp(swarp_config_path=swarp_path,imgpixsize=2400),
            #ImageSaver(output_dir_name="latest"),
            #PhotCalibrator(ref_catalog_generator=summer_photometric_catalog_generator),
            #PhotCalibrator(ref_catalog_generator=summer_backup_photometric_catalog_generator,redo=False),
        ]
    }

    @staticmethod
    def reformat_raw_data(
            img: HDUList,
            path: str
    ) -> [np.array, astropy.io.fits.Header]:
        header = img[0].header
        header["OBSCLASS"] = ["calibration", "science"][header["OBSTYPE"] == "SCIENCE"]
        #print(header['OBSCLASS'])
        header['UTCTIME'] = header['UTCSHUT']
        header['TARGET'] = header['OBSTYPE'].lower()
        #header['TARGET'] = header['FIELDID']
        crd = SkyCoord(ra=img[0].header['RA'], dec=img[0].header['DEC'], unit=(u.deg, u.deg))
        header['RA'] = crd.ra.deg
        header['DEC'] = crd.dec.deg
        header['CRVAL1'] = header['RA']
        header['CRVAL2'] = header['DEC']
        tel_crd = SkyCoord(ra=img[0].header['TELRA'], dec=img[0].header['TELDEC'], unit=(u.deg, u.deg))
        header['TELRA'] = tel_crd.ra.deg
        header['TELDEC'] = tel_crd.dec.deg
        #filters = {'4': 'OPEN', '3': 'r', '1': 'u'}
        header['BZERO'] = 0

        #print(img[0].data.shape)
        img[0].data = img[0].data * 1.0
        #img[0].data[2048, :] = np.nan

        if 'other' in header['FILTERID']:
            header['FILTERID'] = 'r'

        if "CALSTEPS" not in header.keys():
            header["CALSTEPS"] = ""
        header["BASENAME"] = os.path.basename(path)
        header.append(('GAIN', summer_gain, 'Gain in electrons / ADU'), end=True)

        img[0].header = header
        return img[0].data, img[0].header

    # def apply_reduction(self, raw_image_list):
    #     return
