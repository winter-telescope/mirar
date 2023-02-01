import os
from pathlib import Path

from astropy.wcs import WCS

from winterdrp.paths import BASE_NAME_KEY, UNC_IMG_KEY, get_output_dir
from winterdrp.processors.base_processor import (
    BaseDataframeProcessor,
    BaseImageProcessor,
    BaseProcessor,
)
from winterdrp.processors.photometry.utils import (
    aper_photometry,
    get_rms_image,
    make_cutouts,
    make_psf_shifted_array,
    psf_photometry,
)


class BasePhotometry:
    """
    Parent class to run photometry on cutouts
    """

    def __init__(self):
        pass

    def perform_photometry(self, image_cutout, unc_image_cutout):
        raise NotImplementedError


class PSFPhotometry(BasePhotometry):
    """
    Class to run PSF photometry on cutouts
    """

    def __init__(self, psf_filename: str):
        super().__init__()
        self.psf_filename = psf_filename

    def perform_photometry(self, image_cutout, unc_image_cutout):
        psfmodels = make_psf_shifted_array(
            psf_filename=self.psf_filename, cutout_size_psf_phot=image_cutout.shape[0]
        )

        flux, fluxunc, minchi2, xshift, yshift = psf_photometry(
            image_cutout, unc_image_cutout, psfmodels
        )
        return flux, fluxunc, minchi2, xshift, yshift


class AperturePhotometry(BasePhotometry):
    """
    Class to run aperture photometry on cutouts
    """

    def __init__(
        self,
        aper_diameters: float | list[float] = 10.0,
        bkg_in_diameters: float | list[float] = 25.0,
        bkg_out_diameters: float | list[float] = 40.0,
    ):
        if not isinstance(aper_diameters, list):
            aper_diameters = [aper_diameters]
        if not isinstance(bkg_in_diameters, list):
            bkg_in_diameters = [bkg_in_diameters]
        if not isinstance(bkg_out_diameters, list):
            bkg_out_diameters = [bkg_out_diameters]
        super().__init__()

        self.aper_diameters = aper_diameters
        self.bkg_in_diameters = bkg_in_diameters
        self.bkg_out_diameters = bkg_out_diameters

    def perform_photometry(self, image_cutout, unc_image_cutout):
        fluxes, fluxuncs = [], []
        for ind, aper_diam in enumerate(self.aper_diameters):
            flux, fluxunc = aper_photometry(
                image_cutout,
                unc_image_cutout,
                aper_diam,
                self.bkg_in_diameters[ind],
                self.bkg_out_diameters[ind],
            )
            fluxes.append(flux)
            fluxuncs.append(fluxunc)
        return fluxes, fluxuncs


class BasePhotometryProcessor(BaseProcessor):
    """
    Parent processor to run photometry
    """

    def __init__(self, phot_cutout_size: int = 20, output_sub_dir: str = "photometry"):
        super().__init__()
        self.phot_cutout_size = phot_cutout_size
        self.output_sub_dir = output_sub_dir

    def get_sub_output_dir(self) -> Path:
        return Path(get_output_dir(self.output_sub_dir, self.night_sub_dir))

    def get_path(self, name: str | Path) -> Path:
        """
        Get output path for a file of name

        :param name: Name of file
        :return: Full output path
        """
        return self.get_sub_output_dir().joinpath(name)

    def get_filenames(self, data_item):
        raise NotImplementedError

    def get_physical_coordinates(self, data_item):
        raise NotImplementedError

    def generate_cutouts(self, data_item):
        imagename, unc_imagename = self.get_filenames(data_item)
        x, y = self.get_physical_coordinates(data_item)
        image_cutout, unc_image_cutout = make_cutouts(
            image_paths=[imagename, unc_imagename],
            position=(x, y),
            half_size=self.phot_cutout_size,
        )
        return image_cutout, unc_image_cutout


class BaseImagePhotometry(BasePhotometryProcessor, BaseImageProcessor):
    """
    Processor to run photometry on an image
    """

    def __init__(
        self,
        target_ra_key: str = "TARGRA",
        target_dec_key: str = "TARGDEC",
        *args,
        **kwargs
    ):
        super().__init__(*args, **kwargs)
        self.target_ra_key = target_ra_key
        self.target_dec_key = target_dec_key

    def get_filenames(self, image):
        imagename = image.header[BASE_NAME_KEY]
        unc_filename = image.header[UNC_IMG_KEY]
        if not os.path.exists(unc_filename):
            rms_image = get_rms_image(image)
            unc_filename = image[BASE_NAME_KEY] + ".unc"
            self.save_fits(rms_image, path=os.path.join(self.get_path(unc_filename)))
        return imagename, unc_filename

    def get_physical_coordinates(self, image):
        ra, dec = image[self.target_ra_key], image[self.target_dec_key]
        wcs = WCS(image.header)
        x, y = wcs.all_world2pix(ra, dec, 1)
        return x, y


class BaseCandidatePhotometry(BasePhotometryProcessor, BaseDataframeProcessor):
    """
    Processor to run photometry on a candidates table
    """

    def __init__(
        self,
        image_colname="diffimname",
        unc_image_colname="diffuncname",
        psf_file_colname="diffpsfname",
        x_colname="xpeak",
        y_colname="ypeak",
        *args,
        **kwargs
    ):
        super().__init__(*args, **kwargs)
        self.image_colname = image_colname
        self.unc_image_colname = unc_image_colname
        self.psf_file_colname = psf_file_colname
        self.x_colname = x_colname
        self.y_colname = y_colname

    def get_filenames(self, row):
        imagename = row[self.image_colname]
        unc_filename = row[self.unc_image_colname]
        return imagename, unc_filename

    def get_physical_coordinates(self, row):
        x, y = row[self.x_colname], row[self.y_colname]
        return x, y
