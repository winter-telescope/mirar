"""
Module for obtaining a Gaia/2Mass catalog
"""
import logging
from typing import Optional

import astropy.table
import astropy.units as u
from astropy.coordinates import SkyCoord
from astroquery.gaia import Gaia

from mirar.catalog.base_catalog import BaseCatalog
from mirar.utils.ldac_tools import get_table_from_ldac

logger = logging.getLogger(__name__)


class Gaia2Mass(BaseCatalog):
    """
    Crossmatched Gaia/2Mass catalog
    """

    abbreviation = "tmass"

    def __init__(
        self,
        *args,
        filter_name: str = "j",
        ph_qual_cut: bool = False,
        snr_threshold: float = 5,
        trim: bool = False,
        image_catalog_path: Optional[str] = None,
        **kwargs,
    ):
        super().__init__(*args, filter_name=filter_name, **kwargs)
        self.ph_qual_cut = ph_qual_cut
        self.trim = trim
        self.image_catalog_path = image_catalog_path
        self.snr_threshold = snr_threshold

        logger.debug(f"Sextractor catalog path is {self.image_catalog_path}")

    def get_catalog(
        self,
        ra_deg: float,
        dec_deg: float,
    ) -> astropy.table.Table:
        logger.info(
            f"Querying 2MASS - Gaia cross-match around RA {ra_deg:.4f}, "
            f"Dec {dec_deg:.4f} with a radius of {self.search_radius_arcmin:.4f} arcmin"
        )

        cmd = (
            f"SELECT * FROM gaiadr2.gaia_source AS g, "
            f"gaiadr2.tmass_best_neighbour AS tbest, "
            f"gaiadr1.tmass_original_valid AS tmass "
            f"WHERE g.source_id = tbest.source_id "
            f"AND tbest.tmass_oid = tmass.tmass_oid "
            f"AND CONTAINS(POINT('ICRS', g.ra, g.dec), "
            f"CIRCLE('ICRS', {ra_deg:.4f}, {dec_deg:.4f}, "
            f"{self.search_radius_arcmin / 60:.4f}))=1 "
            f"AND tmass.{self.filter_name}_m > {self.min_mag:.2f} "
            f"AND tmass.{self.filter_name}_m < {self.max_mag:.2f} "
            f"AND tbest.number_of_mates=0 "
            f"AND tbest.number_of_neighbours=1"
        )

        if self.ph_qual_cut:
            cmd += "AND tmass.ph_qual='AAA';"
        else:
            cmd += ";"

        job = Gaia.launch_job_async(cmd, dump_to_file=False)
        src_list = job.get_results()
        src_list["ph_qual"] = src_list["ph_qual"].astype(str)
        src_list["ra_errdeg"] = src_list["ra_error"] / 3.6e6
        src_list["dec_errdeg"] = src_list["dec_error"] / 3.6e6
        src_list["FLAGS"] = 0
        src_list["magnitude"] = src_list[f"{self.filter_name.lower()}_m"]
        src_list["magnitude_err"] = src_list[f"{self.filter_name.lower()}_msigcom"]
        logger.info(f"Found {len(src_list)} sources in Gaia")

        src_list = src_list[src_list["magnitude_err"] < 1.086 / self.snr_threshold]
        if self.trim:
            if self.image_catalog_path is None:
                logger.error(
                    "Gaia catalog trimming requested but "
                    "no sextractor catalog path specified."
                )
                raise ValueError

            image_catalog = get_table_from_ldac(self.image_catalog_path)
            src_list = self.trim_catalog(src_list, image_catalog)
            logger.info(f"Trimmed to {len(src_list)} sources in Gaia")

        return src_list

    @staticmethod
    def trim_catalog(
        ref_catalog: astropy.table.Table, image_catalog: astropy.table.Table
    ) -> astropy.table.Table:
        """
        Trim a reference catalog by only taking ref sources within 2 arcseconds of
        image sources

        :param ref_catalog: reference catalog
        :param image_catalog: image catalog
        :return: trimmed ref catalog
        """
        ref_coords = SkyCoord(
            ra=ref_catalog["ra"], dec=ref_catalog["dec"], unit=(u.deg, u.deg)
        )
        image_coords = SkyCoord(
            ra=image_catalog["ALPHAWIN_J2000"],
            dec=image_catalog["DELTAWIN_J2000"],
            unit=(u.deg, u.deg),
        )
        idx, d2d, _ = image_coords.match_to_catalog_sky(ref_coords)
        match_mask = d2d < 2 * u.arcsec
        matched_catalog = ref_catalog[idx[match_mask]]
        return matched_catalog
