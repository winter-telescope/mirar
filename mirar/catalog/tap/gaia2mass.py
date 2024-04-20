"""
Module for obtaining a Gaia/2Mass catalog
"""

import logging

import astropy.table
from astroquery.gaia import Gaia

from mirar.catalog.base.base_gaia import BaseGaia2Mass

logger = logging.getLogger(__name__)

# Disable astroquery warnings
logging.getLogger("astroquery").setLevel(logging.WARNING)


class Gaia2MassTAP(BaseGaia2Mass):
    """
    Crossmatched Gaia/2Mass catalog
    """

    def get_source_table(
        self,
        ra_deg: float,
        dec_deg: float,
    ) -> astropy.table.Table:
        logger.debug(
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
            f"AND tbest.number_of_neighbours=1;"
        )

        job = Gaia.launch_job_async(cmd, dump_to_file=False)
        src_list = job.get_results()

        src_list = self.convert_to_ab_mag(src_list)

        src_list["ph_qual"] = src_list["ph_qual"].astype(str)
        src_list["ra_errdeg"] = src_list["ra_error"] / 3.6e6
        src_list["dec_errdeg"] = src_list["dec_error"] / 3.6e6
        src_list["FLAGS"] = 0
        src_list["magnitude"] = src_list[f"{self.filter_name.lower()}_m"]
        src_list["magnitude_err"] = src_list[f"{self.filter_name.lower()}_msigcom"]
        logger.debug(f"Found {len(src_list)} sources in Gaia")
        return src_list
