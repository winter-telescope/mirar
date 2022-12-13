import logging
from typing import Optional

import astropy.table
import astropy.units as u
from astropy.coordinates import SkyCoord
from astroquery.gaia import Gaia

from winterdrp.catalog.base_catalog import BaseCatalog
from winterdrp.utils.ldac_tools import get_table_from_ldac

logger = logging.getLogger(__name__)


class Gaia2Mass(BaseCatalog):
    abbreviation = "tmass"

    def __init__(
        self,
        search_radius_arcmin: float,
        min_mag: float,
        max_mag: float,
        filter_name: str = "j",
        ph_qual_cut: bool = False,
        trim: bool = False,
        image_catalog_path: Optional[str] = None,
    ):
        super().__init__(search_radius_arcmin, min_mag, max_mag, filter_name)
        self.ph_qual_cut = ph_qual_cut
        self.trim = trim
        self.image_catalog_path = image_catalog_path

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
            f"CIRCLE('ICRS', {ra_deg:.4f}, {dec_deg:.4f}, {self.search_radius_arcmin / 60:.4f}))=1 "
            f"AND tmass.{self.filter_name}_m > {self.min_mag:.2f} "
            f"AND tmass.{self.filter_name}_m < {self.max_mag:.2f} "
            f"AND tbest.number_of_mates=0 "
            f"AND tbest.number_of_neighbours=1"
        )

        if self.ph_qual_cut:
            cmd += f"AND tmass.ph_qual='AAA';"
        else:
            cmd += ";"

        job = Gaia.launch_job_async(cmd, dump_to_file=False)
        t = job.get_results()
        t["ph_qual"] = t["ph_qual"].astype(str)
        t["ra_errdeg"] = t["ra_error"] / 3.6e6
        t["dec_errdeg"] = t["dec_error"] / 3.6e6
        t["FLAGS"] = 0
        t["magnitude"] = t[f"{self.filter_name.lower()}_m"]

        logger.info(f"Found {len(t)} sources in Gaia")
        if self.trim:
            if self.image_catalog_path is None:
                logger.error(
                    "Gaia catalog trimming requested but no sextractor catalog path specified."
                )
                raise ValueError
            else:
                image_catalog = get_table_from_ldac(self.image_catalog_path)
                t = self.trim_catalog(t, image_catalog)
                logger.info(f"Trimmed to {len(t)} sources in Gaia")

        return t

    @staticmethod
    def trim_catalog(ref_catalog, image_catalog):
        ref_coords = SkyCoord(
            ra=ref_catalog["ra"], dec=ref_catalog["dec"], unit=(u.deg, u.deg)
        )
        image_coords = SkyCoord(
            ra=image_catalog["ALPHAWIN_J2000"],
            dec=image_catalog["DELTAWIN_J2000"],
            unit=(u.deg, u.deg),
        )
        idx, d2d, d3d = image_coords.match_to_catalog_sky(ref_coords)
        match_mask = d2d < 2 * u.arcsec
        matched_catalog = ref_catalog[idx[match_mask]]
        return matched_catalog
