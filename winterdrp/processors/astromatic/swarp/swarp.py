import logging
import os
import numpy as np
import astropy.io.fits
from winterdrp.processors.base_processor import BaseProcessor
from winterdrp.paths import get_output_dir, copy_temp_file, get_temp_path, base_name_key
from winterdrp.utils import execute
from winterdrp.processors.astromatic.scamp.scamp import Scamp, scamp_header_key
import shutil
from astropy.wcs import WCS

logger = logging.getLogger(__name__)


# def run_swarp(
#         scamp_list_path: str,
#         scamp_config_path: str,
#         ast_ref_cat_path: str,
#         output_dir: str
# ):
#     swarp_cmd = f"scamp @{scamp_list_path} " \
#                 f"-c {scamp_config_path} " \
#                 f"-ASTREFCAT_NAME {ast_ref_cat_path} " \
#                 # f"-VERBOSE_TYPE QUIET "
#
#     execute(swarp_cmd, output_dir=output_dir)

def run_swarp(
        stack_list_path: str,
        swarp_config_path: str,
        out_path: str,
        weight_list_path: str = None,
        weight_out_path: str = None,
        pixscale: float=None,
        imgpixsize: float=None,
        propogate_headerlist: list=None,
        center_ra: float=None,
        center_dec: float=None,
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
                    f'-RESAMPLE Y -RESAMPLE_DIR {os.path.dirname(out_path)} '\
                    f'-COMBINE Y -COMBINE_TYPE MEDIAN '\
                    f'-SUBTRACT_BACK N'

    if weight_list_path is not None:
        swarp_command += f' -WEIGHT_TYPE MAP_WEIGHT -WEIGHT_IMAGE @{weight_list_path} '

    if weight_out_path is not None:
        swarp_command += f' -WEIGHTOUT_NAME {weight_out_path}'

    if pixscale is not None:
        swarp_command += f' -PIXEL_SCALE {pixscale}'

    if propogate_headerlist is not None:
        swarp_command += f' -COPY_KEYWORDS ' + propogate_headerlist

    if np.logical_and(center_ra is not None, center_dec is not None):
        swarp_command += f' -CENTER {center_ra},{center_dec}'
    print(swarp_command)

    execute(swarp_command)


# def get_scamp_output_head_path(
#         cat_path: str
# ) -> str:
#     return os.path.splitext(cat_path)[0] + ".head"


class Swarp(BaseProcessor):

    base_key = "swarp"

    def __init__(
            self,
            swarp_config_path: str,
            temp_output_sub_dir: str = "swarp",
            pixscale: float = None,
            imgpixsize: float=None,
            propogate_headerlist: list=None,
            center_ra: float = None,
            center_dec: float = None,
            *args,
            **kwargs
    ):
        super(Swarp, self).__init__(*args, **kwargs)
        self.swarp_config = swarp_config_path
        self.temp_output_sub_dir = temp_output_sub_dir
        self.pixscale = pixscale
        self.imgpixsize = imgpixsize
        self.propogate_headerlist = propogate_headerlist
        self.center_ra = center_ra
        self.center_dec = center_dec

    def get_swarp_output_dir(self):
        return get_output_dir(self.temp_output_sub_dir, self.night_sub_dir)

    def _apply_to_images(
            self,
            images: list[np.ndarray],
            headers: list[astropy.io.fits.Header],
    ) -> tuple[list[np.ndarray], list[astropy.io.fits.Header]]:

        swarp_output_dir = self.get_swarp_output_dir()

        try:
            os.makedirs(swarp_output_dir)
        except OSError:
            pass

        swarp_image_list_path = os.path.join(
            swarp_output_dir,
            os.path.splitext(headers[0]["BASENAME"])[0] + "_swarp_img_list.txt"
        )

        swarp_weight_list_path = os.path.join(
            swarp_output_dir,
            os.path.splitext(headers[0]["BASENAME"])[0] + "_swarp_weight_list.txt"
        )

        logger.info(f"Writing file list to {swarp_image_list_path}")

        temp_files = [swarp_image_list_path]

        output_image_path = os.path.join(
            swarp_output_dir,
            os.path.splitext(headers[0]["BASENAME"])[0] + "_stack.fits"
        )

        out_files = []

        all_pixscales = []
        all_imgpixsizes = []
        all_ras = []
        all_decs = []
        with open(swarp_image_list_path, "w") as f, open(swarp_weight_list_path, "w") as g:
            for i, data in enumerate(images):
                header = headers[i]

                if not np.logical_or(self.pixscale is None, self.imgpixsize is None):
                    pass
                else:
                    w = WCS(header)
                    cd11 = header['CD1_1']
                    cd21 = header['CD2_1']
                    cd12 = header['CD1_2']
                    cd22 = header['CD2_2']
                    nxpix = header['NAXIS1']
                    nypix = header['NAXIS2']
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

                print(str(header[scamp_header_key]).replace('\n',''))

                temp_head_path = copy_temp_file(
                    output_dir=swarp_output_dir,
                    file_path=str(header[scamp_header_key]).replace('\n', '')
                )

                temp_img_path = get_temp_path(swarp_output_dir, header["BASENAME"])
                self.save_fits(data, header, temp_img_path)
                temp_mask_path = self.save_mask(data, header, temp_img_path)

                f.write(f"{temp_img_path}\n")
                g.write(f"{temp_mask_path}\n")
                temp_files += [temp_head_path, temp_img_path, temp_mask_path]

                # out_files.append(out_path)

        if self.pixscale is None:
            self.pixscale = np.max(pixscale)
        if self.imgpixsize is None:
            self.imgpixsize = np.max(imgpixsize)
        if self.center_ra is None:
            self.center_ra = np.median(all_ras)
        if self.center_dec is None:
            self.center_dec = np.median(all_decs)

        run_swarp(
            stack_list_path=swarp_image_list_path,
            swarp_config_path=self.swarp_config,
            out_path=output_image_path,
            weight_list_path=swarp_weight_list_path,
            weight_out_path=output_image_path.replace(".fits", ".weight.fits"),
            pixscale=self.pixscale,
            imgpixsize=self.imgpixsize,
            propogate_headerlist=self.propogate_headerlist,
            center_ra=self.center_ra,
            center_dec=self.center_dec
        )

        image, new_header = self.open_fits(output_image_path)

        for key in headers[0]:
            if np.sum([x[key] == headers[0][key] for x in headers]) == len(headers):
                if key not in new_header:
                    new_header[key] = headers[0][key]

        new_header["COADDS"] = np.sum([x["COADDS"] for x in headers])
        # print(new_header["COADDS"], len(headers))
        # print(header["SATURATION"])
        # # Saturation
        #
        # raise

        # fix
        # headers
        # here.....load
        # up
        # comps

        new_header[base_name_key] = os.path.basename(output_image_path)

        return [image], [new_header]

    def check_prerequisites(
            self,
    ):
        check = np.sum([isinstance(x, Scamp) for x in self.preceding_steps])
        if check < 1:
            err = f"{self.__module__} requires {Scamp} as a prerequisite. " \
                  f"However, the following steps were found: {self.preceding_steps}."
            logger.error(err)
            raise ValueError

    # @classmethod
    # def single_catalog(
    #         cls,
    #         catalog: BaseCatalog,
    #         *args,
    #         **kwargs
    # ):
    #
    #     def get_catalog(
    #             header: astropy.io.fits.Header
    #     ) -> BaseCatalog:
    #         return catalog
    #
    #     return cls(get_catalog=get_catalog, *args, **kwargs)
