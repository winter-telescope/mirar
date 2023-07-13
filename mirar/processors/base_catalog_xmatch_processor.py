"""
Base processor for any functionality requiring cross-matching a
reference catalog to a Sextractor catalog
"""
from pathlib import Path
from typing import Callable

import numpy as np
from astropy import units as u
from astropy.coordinates import Angle, SkyCoord
from astropy.table import Table

from mirar.catalog import BaseCatalog
from mirar.data import Image
from mirar.errors import ProcessorError
from mirar.paths import BASE_NAME_KEY, copy_temp_file, get_output_dir, get_output_path
from mirar.processors import BaseImageProcessor
from mirar.processors.astromatic import Sextractor
from mirar.processors.astromatic.sextractor.sextractor import SEXTRACTOR_HEADER_KEY
from mirar.processors.base_processor import PrerequisiteError, logger
from mirar.processors.candidates.utils.regions_writer import write_regions_file
from mirar.utils.ldac_tools import get_table_from_ldac


def default_image_sextractor_catalog_purifier(
    catalog: Table,
    image: Image,
    edge_width_pixels: float = 100.0,
    fwhm_threshold_arcsec: float = 4.0,
) -> Table:
    """
    Default function to purify the photometric image catalog
    """
    x_lower_limit = edge_width_pixels
    x_upper_limit = image.get_data().shape[1] - edge_width_pixels
    y_lower_limit = edge_width_pixels
    y_upper_limit = image.get_data().shape[0] - edge_width_pixels

    clean_mask = (
        (catalog["FLAGS"] == 0)
        & (catalog["FWHM_WORLD"] < fwhm_threshold_arcsec / 3600.0)
        & (catalog["X_IMAGE"] > x_lower_limit)
        & (catalog["X_IMAGE"] < x_upper_limit)
        & (catalog["Y_IMAGE"] > y_lower_limit)
        & (catalog["Y_IMAGE"] < y_upper_limit)
    )

    return catalog[clean_mask]


class BaseProcessorWithCrossMatch(BaseImageProcessor):
    """
    Photometric calibrator processor
    """

    base_key = "astrometrystatswriter"

    def __init__(
        self,
        ref_catalog_generator: Callable[[Image], BaseCatalog],
        crossmatch_radius_arcsec: float,
        required_parameters: list[str],
        sextractor_catalog_purifier: Callable[[Table, Image], Table],
        write_regions: bool = False,
        cache: bool = False,
        temp_output_sub_dir: str = "astrom_stats",
    ):
        super().__init__()
        self.ref_catalog_generator = ref_catalog_generator
        self.crossmatch_radius_arcsec = crossmatch_radius_arcsec
        self.temp_output_sub_dir = temp_output_sub_dir
        self.sextractor_catalog_purifier = sextractor_catalog_purifier
        self.cache = cache
        self.required_parameters = required_parameters

        self.write_regions = write_regions

    def __str__(self) -> str:
        return (
            "Base processor for any functionality requiring "
            "cross-matching a reference catalog to a Sextractor catalog"
        )

    def setup_catalogs(self, image):
        """
        Setup the reference catalog and image catalog
        """
        ref_catalog = self.ref_catalog_generator(image)
        output_dir = get_output_dir(
            dir_root=self.temp_output_sub_dir,
            sub_dir=self.night_sub_dir,
        )
        ref_cat_path = ref_catalog.write_catalog(image, output_dir=output_dir)
        temp_cat_path = copy_temp_file(
            output_dir=Path(output_dir), file_path=image[SEXTRACTOR_HEADER_KEY]
        )

        ref_cat = get_table_from_ldac(ref_cat_path)
        img_cat = get_table_from_ldac(temp_cat_path)

        if self.write_regions:
            self.write_regions_files(image=image, ref_cat=ref_cat, img_cat=img_cat)

        cleaned_img_cat = self.sextractor_catalog_purifier(img_cat, image)

        if not self.cache:
            temp_cat_path.unlink()
            logger.debug(f"Deleted temporary file {temp_cat_path}")

        return ref_cat, img_cat, cleaned_img_cat

    def xmatch_catalogs(
        self, ref_cat: Table, image_cat: Table, crossmatch_radius_arcsec: float
    ) -> (Table, Table, Angle):
        """
        Cross-match the reference catalog to the image catalog
        """
        ref_coords = SkyCoord(ra=ref_cat["ra"], dec=ref_cat["dec"], unit=(u.deg, u.deg))

        img_coords = SkyCoord(
            ra=image_cat["ALPHAWIN_J2000"],
            dec=image_cat["DELTAWIN_J2000"],
            unit=(u.deg, u.deg),
        )

        logger.debug(
            f"Cross-matching {len(ref_cat)} sources in catalog to {len(img_coords)} "
            f"image with radius {crossmatch_radius_arcsec} arcsec."
        )

        idx, d2d, _ = ref_coords.match_to_catalog_sky(img_coords)
        match_mask = d2d < crossmatch_radius_arcsec * u.arcsec
        matched_ref_cat = ref_cat[match_mask]
        matched_img_cat = image_cat[idx[match_mask]]
        logger.debug(
            f"Cross-matched {len(matched_img_cat)} sources from catalog to the image."
        )

        return matched_img_cat, matched_ref_cat, d2d[match_mask]

    def write_regions_files(self, image: Image, ref_cat: Table, img_cat: Table):
        """
        Write regions files for the reference catalog and the image catalog
        """
        ref_coords = SkyCoord(ra=ref_cat["ra"], dec=ref_cat["dec"], unit=(u.deg, u.deg))

        img_coords = SkyCoord(
            ra=img_cat["ALPHAWIN_J2000"],
            dec=img_cat["DELTAWIN_J2000"],
            unit=(u.deg, u.deg),
        )

        clean_img_coords = SkyCoord(
            ra=img_cat["ALPHAWIN_J2000"],
            dec=img_cat["DELTAWIN_J2000"],
            unit=(u.deg, u.deg),
        )

        ref_regions_path = get_output_path(
            base_name=image.header[BASE_NAME_KEY] + "ref.reg",
            dir_root=self.temp_output_sub_dir,
            sub_dir=self.night_sub_dir,
        )
        cleaned_img_regions_path = get_output_path(
            base_name=image.header[BASE_NAME_KEY] + "cleaned_img.reg",
            dir_root=self.temp_output_sub_dir,
            sub_dir=self.night_sub_dir,
        )
        img_regions_path = get_output_path(
            base_name=image.header[BASE_NAME_KEY] + "img.reg",
            dir_root=self.temp_output_sub_dir,
            sub_dir=self.night_sub_dir,
        )

        write_regions_file(
            regions_path=ref_regions_path,
            x_coords=ref_coords.ra.deg,
            y_coords=ref_coords.dec.deg,
            system="wcs",
            region_radius=2.0 / 3600,
        )
        write_regions_file(
            regions_path=cleaned_img_regions_path,
            x_coords=clean_img_coords.ra.deg,
            y_coords=clean_img_coords.dec.deg,
            system="wcs",
            region_radius=2.0 / 3600,
        )
        write_regions_file(
            regions_path=img_regions_path,
            x_coords=img_coords.ra.deg,
            y_coords=img_coords.dec.deg,
            system="wcs",
            region_radius=2.0 / 3600,
        )

    def get_sextractor_module(self) -> Sextractor:
        """
        Get the Sextractor module from the preceding steps
        """
        mask = [isinstance(x, Sextractor) for x in self.preceding_steps]
        return np.array(self.preceding_steps)[mask][-1]

    def check_prerequisites(
        self,
    ):
        mask = [isinstance(x, Sextractor) for x in self.preceding_steps]
        if np.sum(mask) < 1:
            err = (
                f"{self.__module__} requires {Sextractor} as a prerequisite. "
                f"However, the following steps were found: {self.preceding_steps}."
            )
            logger.error(err)
            raise PrerequisiteError(err)

        sextractor_param_path = self.get_sextractor_module().parameters_name

        logger.debug(f"Checking file {sextractor_param_path}")

        with open(sextractor_param_path, "rb") as param_file:
            sextractor_params = [
                x.strip().decode() for x in param_file.readlines() if len(x.strip()) > 0
            ]
            sextractor_params = [
                x.split("(")[0] for x in sextractor_params if x[0] not in ["#"]
            ]

        for param in self.required_parameters:
            if param not in sextractor_params:
                err = (
                    f"Missing parameter: {self.__module__} requires {param} to run, "
                    f"but this parameter was not found in sextractor config file "
                    f"'{sextractor_param_path}' . "
                    f"Please add the parameter to this list!"
                )
                logger.error(err)
                raise PrerequisiteError(err)

    def get_sextractor_apertures(self) -> list[float]:
        """
        Function to extract sextractor aperture sizes from config file
        Returns:
        """
        sextractor_config_path = self.get_sextractor_module().config

        with open(sextractor_config_path, "rb") as sextractor_config_file:
            aperture_lines = [
                x.decode()
                for x in sextractor_config_file.readlines()
                if np.logical_and(b"PHOT_APERTURES" in x, x.decode()[0] != "#")
            ]

        if len(aperture_lines) > 1:
            err = (
                f"The config file {sextractor_config_path} has "
                f"multiple entries for PHOT_APERTURES."
            )
            logger.error(err)
            raise ProcessorError(err)

        line = aperture_lines[0].replace("PHOT_APERTURES", " ").split("#")[0]

        return [float(x) for x in line.split(",") if x not in [""]]
