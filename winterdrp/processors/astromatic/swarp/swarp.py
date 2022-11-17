import logging
import os
import numpy as np
import astropy.io.fits
from winterdrp.processors.base_processor import BaseImageProcessor
from winterdrp.paths import get_output_dir, copy_temp_file, get_temp_path, base_name_key, latest_mask_save_key, \
    raw_img_key
from winterdrp.utils import execute
from winterdrp.processors.astromatic.scamp.scamp import Scamp, scamp_header_key
from astropy.wcs import WCS
from winterdrp.data import ImageBatch

logger = logging.getLogger(__name__)


def run_swarp(
        stack_list_path: str,
        swarp_config_path: str,
        out_path: str,
        weight_list_path: str = None,
        weight_out_path: str = None,
        pixscale: float = None,
        x_imgpixsize: float = None,
        y_imgpixsize: float = None,
        propogate_headerlist: list = None,
        center_ra: float = None,
        center_dec: float = None,
        combine: bool = True,
        gain: float = None,
        subtract_bkg: bool = False
):  # resample and stack images with swarp
    """
    Resample and stack given images with swarp
    Parameters
    ----------
    stack_list_path : string
        Name of file containing the names of files to be stacked
        One file name per line
    weight_list_path : string
        Name of file containing the names of weight files to be stacked
        One file name per line
    swarp_config_path: str
        Path of Swarp config file
    out_path : string
        Path of stacked output file
    weight_out_path: str
        Path of output weight image
    """

    swarp_command = f'swarp -c {swarp_config_path} ' \
                    f'@{stack_list_path} ' \
                    f'-IMAGEOUT_NAME {out_path} ' \
                    f'-RESAMPLE Y -RESAMPLE_DIR {os.path.dirname(out_path)} '

    if subtract_bkg:
        swarp_command += f'-SUBTRACT_BACK Y '
    else:
        swarp_command += f'-SUBTRACT_BACK N '
    if combine:
        swarp_command += f'-COMBINE Y -COMBINE_TYPE MEDIAN '
    else:
        swarp_command += f'-COMBINE N '

    if weight_list_path is not None:
        swarp_command += f' -WEIGHT_TYPE MAP_WEIGHT -WEIGHT_IMAGE @{weight_list_path} '

    if weight_out_path is not None:
        swarp_command += f' -WEIGHTOUT_NAME {weight_out_path}'

    if pixscale is not None:
        swarp_command += f' -PIXELSCALE_TYPE MANUAL -PIXEL_SCALE {pixscale}'

    if propogate_headerlist is not None:
        swarp_command += f' -COPY_KEYWORDS '
        for keyword in propogate_headerlist:
            swarp_command += f'{keyword},'

        # remove final comma
        swarp_command = swarp_command[:-1]

    if np.logical_and(center_ra is not None, center_dec is not None):
        swarp_command += f' -CENTER {center_ra},{center_dec}'

    if x_imgpixsize is not None:
        swarp_command += f' -IMAGE_SIZE {x_imgpixsize}'
        if y_imgpixsize is not None:
            swarp_command += f',{y_imgpixsize}'

    if gain is not None:
        swarp_command += f' -GAIN {gain}'
    execute(swarp_command)


class Swarp(BaseImageProcessor):
    base_key = "swarp"

    def __init__(
            self,
            swarp_config_path: str,
            temp_output_sub_dir: str = "swarp",
            pixscale: float = None,
            x_imgpixsize: float = None,
            y_imgpixsize: float = None,
            propogate_headerlist: list = None,
            center_ra: float = None,
            center_dec: float = None,
            gain: float = None,
            include_scamp: bool = True,
            combine: bool = False,
            cache: bool = False,
            subtract_bkg: bool = False,
            *args,
            **kwargs
    ):
        super(Swarp, self).__init__(*args, **kwargs)
        self.swarp_config = swarp_config_path
        self.temp_output_sub_dir = temp_output_sub_dir
        self.pixscale = pixscale
        self.propogate_headerlist = propogate_headerlist
        self.center_ra = center_ra
        self.center_dec = center_dec
        self.x_imgpixsize = x_imgpixsize
        self.y_imgpixsize = y_imgpixsize
        self.include_scamp = include_scamp
        self.combine = combine
        self.cache = cache
        self.gain = gain
        self.subtract_bkg = subtract_bkg

    def __str__(self) -> str:
        return f"Processor to apply swarp to images, stacking them together."

    def get_swarp_output_dir(self):
        return get_output_dir(self.temp_output_sub_dir, self.night_sub_dir)

    def _apply_to_images(
            self,
            batch: ImageBatch,
    ) -> ImageBatch:

        swarp_output_dir = self.get_swarp_output_dir()

        try:
            os.makedirs(swarp_output_dir)
        except OSError:
            pass

        swarp_image_list_path = os.path.join(
            swarp_output_dir,
            os.path.splitext(batch[0][base_name_key])[0] + "_swarp_img_list.txt"
        )

        swarp_weight_list_path = os.path.join(
            swarp_output_dir,
            os.path.splitext(batch[0][base_name_key])[0] + "_swarp_weight_list.txt"
        )
        logger.debug(f"Writing file list to {swarp_image_list_path}")

        temp_files = [swarp_image_list_path, swarp_weight_list_path]

        # If swarp is run with combine -N option, it just sets the output name as inpname+.resamp.fits
        if self.combine:
            output_image_path = os.path.join(
                swarp_output_dir,
                os.path.splitext(batch[0][base_name_key])[0] + "_stack.fits"
            )
        else:
            output_image_path = os.path.join(
                swarp_output_dir,
                os.path.splitext(batch[0][base_name_key])[0] + ".resamp.fits"
            )
        logger.debug(f"Saving to {output_image_path}")

        all_pixscales = []
        all_imgpixsizes = []
        all_ras = []
        all_decs = []
        with open(swarp_image_list_path, "w") as f, open(swarp_weight_list_path, "w") as g:
            for image in batch:

                pixscale_to_use = None
                x_imgpixsize_to_use = None
                y_imgpixsize_to_use = None
                center_ra_to_use = None
                center_dec_to_use = None

                if not (
                        (self.pixscale is None) |
                        (self.x_imgpixsize is None) |
                        (self.y_imgpixsize is None) |
                        (self.center_ra is None) |
                        (self.center_dec is None)
                ):
                    pixscale_to_use = self.pixscale
                    x_imgpixsize_to_use = self.x_imgpixsize
                    y_imgpixsize_to_use = self.y_imgpixsize
                    center_ra_to_use = self.center_ra
                    center_dec_to_use = self.center_dec
                else:
                    w = WCS(image.get_header())

                    cd11 = image['CD1_1']
                    cd21 = image['CD2_1']
                    cd12 = image['CD1_2']
                    cd22 = image['CD2_2']

                    nxpix = image['NAXIS1']
                    nypix = image['NAXIS2']

                    image_x_cen = nxpix / 2
                    image_y_cen = nypix / 2

                    [ra, dec] = w.all_pix2world(image_x_cen, image_y_cen, 1)

                    xscale = np.sqrt(cd11 ** 2 + cd21 ** 2)
                    yscale = np.sqrt(cd12 ** 2 + cd22 ** 2)

                    pixscale = xscale * 3600
                    imgpixsize = max(nxpix, nypix)

                    all_pixscales.append(pixscale)
                    all_imgpixsizes.append(imgpixsize)
                    all_ras.append(ra)
                    all_decs.append(dec)

                logger.debug(f"{all_ras}, {all_decs}")

                if self.include_scamp:
                    temp_head_path = copy_temp_file(
                        output_dir=swarp_output_dir,
                        file_path=str(image[scamp_header_key]).replace('\n', '')
                    )

                temp_img_path = get_temp_path(swarp_output_dir, image[base_name_key])

                self.save_fits(image, temp_img_path)

                temp_mask_path = self.save_mask(image, temp_img_path)

                f.write(f"{temp_img_path}\n")
                g.write(f"{temp_mask_path}\n")

                if self.include_scamp:
                    temp_files += [temp_head_path, temp_img_path, temp_mask_path]
                else:
                    temp_files += [temp_img_path, temp_mask_path]

        if pixscale_to_use is None:
            pixscale_to_use = np.max(all_pixscales)
        if x_imgpixsize_to_use is None:
            x_imgpixsize_to_use = np.max(all_imgpixsizes)
        if y_imgpixsize_to_use is None:
            y_imgpixsize_to_use = np.max(all_imgpixsizes)
        if center_ra_to_use is None:
            center_ra_to_use = np.median(all_ras)
        if center_dec_to_use is None:
            center_dec_to_use = np.median(all_decs)

        logger.debug(f"{self.center_ra}, {center_ra_to_use}")

        output_image_weight_path = output_image_path.replace(".fits", ".weight.fits")

        run_swarp(
            stack_list_path=swarp_image_list_path,
            swarp_config_path=self.swarp_config,
            out_path=output_image_path,
            weight_list_path=swarp_weight_list_path,
            weight_out_path=output_image_weight_path,
            pixscale=pixscale_to_use,
            x_imgpixsize=x_imgpixsize_to_use,
            y_imgpixsize=y_imgpixsize_to_use,
            propogate_headerlist=self.propogate_headerlist,
            center_ra=center_ra_to_use,
            center_dec=center_dec_to_use,
            combine=self.combine,
            gain=self.gain,
            subtract_bkg=self.subtract_bkg
        )

        # Check if output image exists if combine is no
        if not self.combine:

            temp_output_image_path = get_temp_path(
                swarp_output_dir,
                os.path.splitext(image[0][base_name_key])[0] + ".resamp.fits"
            )
            temp_output_image_weight_path = temp_output_image_path.replace(".fits", ".weight.fits")

            if os.path.exists(temp_output_image_path):
                os.rename(temp_output_image_path, output_image_path)
                os.rename(temp_output_image_weight_path, output_image_weight_path)
            else:
                err = f'Swarp seems to have misbehaved, and not made the correct output file {output_image_path}'
                logger.error(err)
                raise ValueError

        new_image = self.open_fits(output_image_path)

        for key in batch[0]:
            if np.sum([x[key] == batch[0][key] for x in batch]) == len(batch):
                if key not in new_image:
                    new_image[key] = batch[0][key]

        new_image["COADDS"] = np.sum([x["COADDS"] for x in batch])

        if not self.cache:
            for temp_file in temp_files:
                os.remove(temp_file)
                logger.debug(f"Deleted temporary file {temp_file}")

        new_image[raw_img_key] = ",".join([x[raw_img_key] for x in batch])
        new_image[base_name_key] = os.path.basename(output_image_path)
        new_image[latest_mask_save_key] = os.path.basename(output_image_weight_path)
        return ImageBatch([new_image])

    def check_prerequisites(
            self,
    ):
        check = np.sum([isinstance(x, Scamp) for x in self.preceding_steps])
        if check < 1:
            err = f"{self.__module__} requires {Scamp} as a prerequisite. " \
                  f"However, the following steps were found: {self.preceding_steps}."
            logger.error(err)
            raise ValueError
