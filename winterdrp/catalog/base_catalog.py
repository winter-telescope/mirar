import logging
import os
from abc import ABC

import astropy.io.fits
import astropy.table
import astropy.units as u
from astropy.coordinates import SkyCoord
from astropy.table import Table
from astroquery.vizier import Vizier
from penquins import Kowalski

from winterdrp.paths import base_name_key
from winterdrp.utils.ldac_tools import save_table_as_ldac

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
        filter_name: str,
    ):
        self.search_radius_arcmin = search_radius_arcmin
        self.min_mag = min_mag
        self.max_mag = max_mag
        self.filter_name = filter_name

    @staticmethod
    def get_catalog(ra_deg: float, dec_deg: float) -> astropy.table.Table:
        raise NotImplementedError()

    def write_catalog(self, header: astropy.io.fits.Header, output_dir: str) -> str:
        ra_deg = header["CRVAL1"]
        dec_deg = header["CRVAL2"]

        base_name = os.path.basename(header[base_name_key])

        cat = self.get_catalog(ra_deg=ra_deg, dec_deg=dec_deg)

        output_path = self.get_output_path(output_dir, base_name) + ".ldac"

        if os.path.exists(output_path):
            os.remove(output_path)

        logger.info(f"Saving catalog to {output_path}")

        save_table_as_ldac(cat, output_path)

        return output_path

    def get_output_path(self, output_dir: str, base_name: str) -> str:
        cat_base_name = os.path.splitext(base_name)[0] + f".{self.abbreviation}.cat"
        return os.path.join(output_dir, cat_base_name)

    def get_catalog_from_header(self, header: astropy.io.fits.header) -> astropy.table:
        ra_deg = header["CRVAL1"]
        dec_deg = header["CRVAL2"]

        cat = self.get_catalog(ra_deg=ra_deg, dec_deg=dec_deg)

        return cat


class VizierCatalog(BaseCatalog, ABC):
    @property
    def catalog_vizier_code(self):
        raise NotImplementedError()

    @property
    def ra_key(self):
        raise NotImplementedError()

    @property
    def dec_key(self):
        raise NotImplementedError()

    def __init__(
        self,
        search_radius_arcmin: float,
        min_mag: float,
        max_mag: float,
        filter_name: str,
        snr_threshold: float = 3.0,
    ):
        super().__init__(search_radius_arcmin, min_mag, max_mag, filter_name)
        self.snr_threshold = snr_threshold

    def get_mag_key(self):
        return f"{self.filter_name}mag"

    def get_mag_error_key(self):
        return f"e_{self.get_mag_key()}"

    def get_catalog(self, ra_deg: float, dec_deg: float) -> astropy.table.Table:

        logger.info(
            f"Querying {self.abbreviation} catalog around RA {ra_deg:.4f}, "
            f"Dec {dec_deg:.4f} with a radius of {self.search_radius_arcmin:.4f} arcmin"
        )

        v = Vizier(
            columns=["*"],
            column_filters={
                f"{self.get_mag_key()}": f"< {self.max_mag}",
                f"{self.get_mag_error_key()}": "<%.3f" % (1.086 / self.snr_threshold),
            },
            row_limit=-1,
        )

        query = v.query_region(
            SkyCoord(ra=ra_deg, dec=dec_deg, unit=(u.deg, u.deg)),
            radius=str(self.search_radius_arcmin) + "m",
            catalog=self.catalog_vizier_code,
            cache=False,
        )

        if len(query) == 0:
            logger.info(f"No matches found in the given radius in {self.abbreviation}")
            t = Table()
            self.check_coverage(ra_deg, dec_deg)

        else:
            t = query[0]
            t["ra"] = t[self.ra_key]
            t["dec"] = t[self.dec_key]
            t["magnitude"] = t[self.get_mag_key()]
            logger.info(
                f"{len(t)} matches found in the given radius in {self.abbreviation}"
            )
            t.meta["description"] = ""
        return t

    @staticmethod
    def check_coverage(ra_deg: float, dec_deg: float):
        pass


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

    def __init__(
        self, search_radius_arcsec: float = 10, num_sources: int = 1, *args, **kwargs
    ):
        super(BaseXMatchCatalog, self).__init__(*args, **kwargs)
        self.search_radius_arcsec = search_radius_arcsec
        self.num_sources = num_sources

    def query(self, coords: dict) -> dict:
        raise NotImplementedError


class BaseKowalskiXMatch(BaseXMatchCatalog):
    def __init__(
        self, kowalski: Kowalski = None, max_time_ms: float = 10000, *args, **kwargs
    ):
        super(BaseKowalskiXMatch, self).__init__(*args, **kwargs)
        self.max_time_ms = max_time_ms
        self.kowalski = kowalski

    def get_kowalski(self) -> Kowalski:
        protocol, host, port = "https", "kowalski.caltech.edu", 443

        token_kowalski = os.environ.get("kowalski_token")

        if token_kowalski is not None:
            logger.debug("Using kowalski token")

            k = Kowalski(token=token_kowalski, protocol=protocol, host=host, port=port)

        else:

            username_kowalski = os.environ.get("kowalski_user")
            password_kowalski = os.environ.get("kowalski_pwd")

            if username_kowalski is None:
                err = "Kowalski username not provided, please run export kowalski_user=<user>"
                logger.error(err)
                raise ValueError
            if password_kowalski is None:
                err = "Kowalski password not provided, please run export KOWALSKI_PWD=<user>"
                logger.error(err)
                raise ValueError

            k = Kowalski(
                username=username_kowalski,
                password=password_kowalski,
                protocol=protocol,
                host=host,
                port=port,
            )

        connection_ok = k.ping()
        logger.info(f"Connection OK?: {connection_ok}")

        return k

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
                },
            },
            "kwargs": {
                "max_time_ms": self.max_time_ms,
                "limit": self.num_sources,
            },
        }
        logger.info(f"Kowalski is {self.kowalski}")
        response = self.kowalski.query(query=query)
        data = response.get("data")
        return data[self.catalog_name]

    def query(self, coords) -> dict:
        if self.kowalski is None:
            self.kowalski = self.get_kowalski()
        logger.info("Querying kowalski")
        data = self.near_query_kowalski(coords)
        return data
