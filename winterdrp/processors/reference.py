import logging
import os.path

import numpy as np
from winterdrp.processors.base_processor import BaseImageProcessor
from winterdrp.paths import get_output_dir,  base_name_key, norm_psfex_header_key, raw_img_key, ref_img_key
from winterdrp.references.base_reference_generator import BaseReferenceGenerator
from collections.abc import Callable
from winterdrp.processors.astromatic.swarp.swarp import Swarp
from winterdrp.processors.astromatic.sextractor.sextractor import Sextractor
from winterdrp.processors.astromatic.psfex.psfex import PSFex
from winterdrp.io import open_fits, save_to_path
from astropy.wcs import WCS
from winterdrp.data import ImageBatch, Image

logger = logging.getLogger(__name__)


class Reference(BaseImageProcessor):
    base_key = "REFPREP"

    def __init__(self,
                 ref_image_generator: Callable[..., BaseReferenceGenerator],
                 ref_swarp_resampler: Callable[..., Swarp],
                 ref_sextractor: Callable[..., Sextractor],
                 ref_psfex: Callable[..., PSFex],
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

    @staticmethod
    def get_image_header_params(image: Image):
        header = image.get_header()

        wcs = WCS(header)

        image_x_cent = header['NAXIS1'] / 2
        image_y_cent = header['NAXIS2'] / 2
        x_imgsize = header['NAXIS1']
        y_imgsize = header['NAXIS2']

        [ra_cent, dec_cent] = wcs.all_pix2world(image_x_cent, image_y_cent, 1)

        pixscale = np.abs(header['CD1_1']) * 3600
        gain = header['GAIN']

        return image_x_cent, image_y_cent, ra_cent, dec_cent, pixscale, x_imgsize, y_imgsize, gain

    def _apply_to_images(
            self,
            batch: ImageBatch,
    ) -> ImageBatch:

        try:
            os.makedirs(self.get_sub_output_dir())
        except OSError:
            pass

        new_batch = ImageBatch()

        for ind, image in enumerate(batch):

            ref_image = self.ref_image_generator(image)
            ref_image_path = ref_image.write_reference(header=image, output_dir=self.get_sub_output_dir())

            ref_data, ref_header = open_fits(ref_image_path)

            if not (base_name_key in ref_header.keys()):
                logger.debug(os.path.basename(ref_image_path))
                ref_header[base_name_key] = os.path.basename(ref_image_path)

            if not (raw_img_key in ref_header.keys()):
                ref_header[raw_img_key] = ref_image_path

            ref_image = Image(ref_data, ref_header)

            # ref_header[base_name_key] = ref_header[base_name_key] + '_ref'
            ref_gain = ref_header['GAIN']

            sci_x_cent, sci_y_cent = image['NAXIS1'] / 2, image['NAXIS2'] / 2

            header = image.get_header()
            wcs = WCS(header)

            [sci_ra_cent, sci_dec_cent] = wcs.all_pix2world(sci_x_cent, sci_y_cent, 1)
            sci_pixscale = np.abs(image['CD1_1']) * 3600
            sci_x_imgsize, sci_y_imgsize = image['NAXIS1'], image['NAXIS2']
            sci_gain = image['GAIN']

            propogate_headerlist = ['TMC_ZP', 'TMC_ZPSD']
            ref_resampler = self.ref_swarp_resampler(
                pixscale=sci_pixscale,
                x_imgpixsize=sci_x_imgsize,
                y_imgpixsize=sci_y_imgsize,
                center_ra=sci_ra_cent,
                center_dec=sci_dec_cent,
                propogate_headerlist=propogate_headerlist,
                temp_output_sub_dir=self.temp_output_subtract_dir,
                night_sub_dir=self.night_sub_dir,
                include_scamp=False,
                combine=False,
                gain=ref_gain,
                subtract_bkg=True
            )

            ref_resampler.set_night(night_sub_dir=self.night_sub_dir)
            resampled_ref_img_batch = ref_resampler.apply(ImageBatch([ref_image]))

            save_to_path(resampled_ref_img_batch,
                         path=os.path.join(self.get_sub_output_dir(), resampled_ref_img_batch[0][base_name_key]))

            ref_resamp_x_cent, ref_resamp_y_cent, ref_resamp_ra_cent, ref_resamp_dec_cent, \
            ref_resamp_pixscale, ref_resamp_x_imgsize, ref_resamp_y_imgsize, \
            ref_resamp_gain = self.get_image_header_params(resampled_ref_header)

            sci_gain = image['GAIN']
            sci_resampler = self.ref_swarp_resampler(
                pixscale=ref_resamp_pixscale,
                x_imgpixsize=ref_resamp_x_imgsize,
                y_imgpixsize=ref_resamp_y_imgsize,
                center_ra=ref_resamp_ra_cent,
                center_dec=ref_resamp_dec_cent,
                propogate_headerlist=propogate_headerlist,
                temp_output_sub_dir=self.temp_output_subtract_dir,
                night_sub_dir=self.night_sub_dir,
                include_scamp=False,
                combine=False,
                gain=sci_gain,
                subtract_bkg=True
            )
            sci_resampler.set_night(night_sub_dir=self.night_sub_dir)

            resampled_sci_image_batch = sci_resampler.apply(ImageBatch([image]))

            ref_sextractor = self.ref_sextractor(
                output_sub_dir=self.temp_output_subtract_dir,
                gain=ref_resamp_gain
            )
            ref_sextractor.set_night(night_sub_dir=self.night_sub_dir)

            resampled_ref_sextractor_image_batch = ref_sextractor.apply(resampled_sci_image_batch)

            save_to_path(image=resampled_ref_sextractor_image_batch[0],
                         path=os.path.join(self.get_sub_output_dir()))
            logger.info(
                f"Saved reference image to {os.path.join(self.get_sub_output_dir(), resampled_ref_sex_header[base_name_key])}")

            sci_resamp_x_cent, sci_resamp_y_cent, sci_resamp_ra_cent, sci_resamp_dec_cent, \
            sci_resamp_pixscale, sci_resamp_x_imgsize, sci_resamp_y_imgsize, \
            sci_resamp_gain = self.get_image_header_params(resampled_sci_header)

            sci_sextractor = self.ref_sextractor(
                output_sub_dir=self.temp_output_subtract_dir,
                gain=sci_resamp_gain
            )

            [[resampled_sci_sex_image], [resampled_sci_sex_header]] \
                = ref_sextractor.apply([[resampled_sci_image], [resampled_sci_header]])
            save_to_path(data=resampled_sci_sex_image, header=resampled_sci_sex_header,
                         path=os.path.join(self.get_sub_output_dir(), resampled_sci_sex_header[base_name_key]))

            ref_psfex = self.ref_psfex(output_sub_dir=self.temp_output_subtract_dir, norm_fits=True)

            [[resampled_ref_sextractor_image_batch], [resampled_ref_sex_psf_header]] = ref_psfex.apply(
                [[resampled_ref_sextractor_image_batch], [resampled_ref_sex_header]])

            resampled_sci_sex_header["REFPSF"] = resampled_ref_sex_psf_header[norm_psfex_header_key]
            resampled_sci_sex_header[ref_img_key] = resampled_ref_sex_psf_header[base_name_key]

            new_batch.append(resampled_sci_sex_image)
        return new_batch
