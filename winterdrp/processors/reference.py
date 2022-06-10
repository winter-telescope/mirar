import logging
import os.path

import numpy as np
from winterdrp.processors.base_processor import BaseProcessor
from astropy.io import fits
from winterdrp.paths import get_output_dir, copy_temp_file, base_name_key, sextractor_header_key, latest_mask_save_key
from winterdrp.references.base_reference_generator import BaseReferenceGenerator
from collections.abc import Callable
from winterdrp.processors.astromatic.swarp.swarp import Swarp
from winterdrp.processors.astromatic.sextractor.sextractor import Sextractor
from winterdrp.processors.astromatic.psfex.psfex import PSFex
from winterdrp.io import open_fits
from astropy.wcs import WCS

logger = logging.getLogger(__name__)


class Reference(BaseProcessor):

    def __init__(self,
                 ref_image_generator: Callable[[fits.Header], BaseReferenceGenerator],
                 ref_swarp_resampler: Callable[Swarp],
                 ref_sextractor: Callable[Sextractor],
                 ref_psfex: Callable[PSFex],
                 temp_output_subtract_dir: str = "subtract",
                 *args,
                 **kwargs
                 ):
        super(Reference, self).__init__(*args, **kwargs)
        self.ref_image_generator = ref_image_generator
        self.ref_swarp_resampler = ref_swarp_resampler
        self.ref_sextractor = ref_sextractor
        self.ref_psfex = ref_psfex
        self.temp_output_subtract_dir = temp_output_subtract_dir

    def get_sub_output_dir(self):
        return get_output_dir(self.temp_output_subtract_dir, self.night_sub_dir)

    def _apply_to_images(
            self,
            images: list[np.ndarray],
            headers: list[fits.Header],
    ) -> tuple[list[np.ndarray], list[fits.Header]]:
        for header in headers:
            ref_image = self.ref_image_generator(header)
            ref_image_path = ref_image.write_reference(header=header)
            ref_data, ref_header = open_fits(ref_image_path)

            if not (base_name_key in ref_header.keys()):
                ref_data[base_name_key] = os.path.basename(ref_image_path)

            ref_data[base_name_key] = ref_data[base_name_key] + '_ref'

            sci_x_cent, sci_y_cent = header['NAXIS1'] / 2, header['NAXIS2'] / 2
            wcs = WCS(header)
            [sci_ra_cent, sci_dec_cent] = wcs.all_pix2world(sci_x_cent, sci_y_cent, 1)
            sci_pixscale = np.abs(header['CD1_1']) * 3600
            sci_x_imgsize, sci_y_imgsize = header['NAXIS1'], header['NAXIS2']
            sci_gain = header['GAIN']

            propogate_headerlist = ['TMC_ZP', 'TMC_ZPSD']
            ref_resampler = self.ref_swarp_resampler(pixsize=sci_pixscale,
                                                     x_imgpixsize=sci_x_imgsize,
                                                     y_imgpixsize=sci_y_imgsize,
                                                     center_ra=sci_ra_cent,
                                                     center_dec=sci_dec_cent,
                                                     propogate_headerlist=propogate_headerlist,
                                                     temp_output_sub_dir=self.temp_output_subtract_dir
                                                     )

            resampled_ref_image, resampled_ref_header = ref_resampler.apply(images=[ref_data],
                                                                            headers=[ref_header])

            ref_resamp_gain = resampled_ref_header['GAIN']
            ref_sextractor = self.ref_sextractor(output_sub_dir=self.temp_output_subtract_dir,
                                                 gain=ref_resamp_gain
                                                 )

            resampled_ref_sex_image, resampled_ref_sex_header \
                = ref_sextractor.apply(images=[resampled_ref_image], headers=[resampled_ref_header])

            ref_psfex = self.ref_psfex(output_sub_dir=self.temp_output_subtract_dir, norm_fits=True)

            resampled_ref_sex_image, resampled_ref_sex_psf_header = ref_psfex.apply(resampled_ref_sex_image,
                                                                                    resampled_ref_sex_header)

            return resampled_ref_sex_image, resampled_ref_sex_psf_header
