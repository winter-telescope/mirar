import astropy.table
import os
import logging
import astropy.io.fits
import pandas as pd

from winterdrp.utils.ldac_tools import save_table_as_ldac
from winterdrp.paths import base_name_key
from penquins import Kowalski
import numpy as np

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
            filter_name: str
    ):
        self.search_radius_arcmin = search_radius_arcmin
        self.min_mag = min_mag
        self.max_mag = max_mag
        self.filter_name = filter_name

    @staticmethod
    def get_catalog(
            ra_deg: float,
            dec_deg: float
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
            dec_deg=dec_deg
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
            dec_deg=dec_deg
        )

        return cat


class BaseXMatchCatalog:

    @property
    def abbreviation(self):
        raise NotImplementedError

    @property
    def catalog_name(self):
        raise NotImplementedError

    @property
    def projection(self):
        raise NotImplementedError

    @property
    def column_names(self):
        raise NotImplementedError

    @property
    def column_dtypes(self):
        raise NotImplementedError

    def __init__(self,
                 search_radius_arcsec: float = 10,
                 num_sources: int = 1,
                 *args,
                 **kwargs
                 ):
        super(BaseXMatchCatalog, self).__init__(*args, **kwargs)
        self.search_radius_arcsec = search_radius_arcsec
        self.num_sources = num_sources

    def query(self, coords: dict) -> dict:
        raise NotImplementedError


class BaseKowalskiXMatch(BaseXMatchCatalog):

    def __init__(self,
                 kowalski: Kowalski,
                 max_time_ms: float = 10000,
                 *args,
                 **kwargs):
        super(BaseKowalskiXMatch, self).__init__(*args,**kwargs)
        self.kowalski = kowalski
        self.max_time_ms = max_time_ms

    def near_query_kowalski(self, coords: dict) -> dict:
        query = {
            "query_type": "near",
            "query": {
                "max_distance": self.search_radius_arcsec,
                "distance_units": "arcsec",
                "radec": coords,
                "catalogs": {
                    f"{self.catalog_name}": {
                        "filter": {},
                        "projection": self.projection,
                    }
                }
            },
            "kwargs": {
                "max_time_ms": self.max_time_ms,
                "limit": self.num_sources,
            },
        }
        logger.info(f'Kowalski is {self.kowalski}')
        response = self.kowalski.query(query=query)
        data = response.get("data")
        return data[self.catalog_name]

    def query(self, coords) -> dict:
        logger.info('Querying kowalski')
        data = self.near_query_kowalski(coords)
        return data
