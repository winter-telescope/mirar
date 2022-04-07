import astropy.table
import os
import logging
import astropy.io.fits
from winterdrp.utils.ldac_tools import save_table_as_ldac
from winterdrp.paths import base_name_key

logger = logging.getLogger(__name__)


class BaseCatalog:

    @property
    def abbreviation(self):
        raise NotImplementedError()

    def __init__(
            self,
            search_radius_arcmin: float,
            min_mag: float,
            max_mag: float,
    ):
        self.search_radius_arcmin = search_radius_arcmin
        self.min_mag = min_mag
        self.max_mag = max_mag

    @staticmethod
    def get_catalog(
            ra_deg: float,
            dec_deg: float,
            search_radius_arcmin: float,
            min_mag: float,
            max_mag: float
    ) -> astropy.table.Table:
        raise NotImplementedError()

    def write_catalog(
            self,
            header: astropy.io.fits.Header,
            output_dir: str
    ) -> str:
        ra_deg = header['CRVAL1']
        dec_deg = header['CRVAL2']

        base_name = os.path.basename(header[base_name_key])

        cat = self.get_catalog(
            ra_deg=ra_deg,
            dec_deg=dec_deg,
            search_radius_arcmin=self.search_radius_arcmin,
            min_mag=self.min_mag,
            max_mag=self.max_mag
        )

        output_path = self.get_output_path(output_dir, base_name) + ".ldac"

        if os.path.exists(output_path):
            os.remove(output_path)

        logger.info(f"Saving catalog to {output_path}")

        save_table_as_ldac(cat, output_path)

        return output_path

    def get_output_path(
            self,
            output_dir: str,
            base_name: str
    ) -> str:
        cat_base_name = os.path.splitext(base_name)[0] + f".{self.abbreviation}.cat"
        return os.path.join(output_dir, cat_base_name)

    def get_catalog_from_header(
            self,
            header: astropy.io.fits.header
            ) -> astropy.table:
        ra_deg = header['CRVAL1']
        dec_deg = header['CRVAL2']

        cat = self.get_catalog(
            ra_deg=ra_deg,
            dec_deg=dec_deg,
            search_radius_arcmin=self.search_radius_arcmin,
            min_mag=self.min_mag,
            max_mag=self.max_mag
        )

        return cat
