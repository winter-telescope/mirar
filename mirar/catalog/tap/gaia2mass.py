"""
Module for obtaining a Gaia/2Mass catalog
"""

import logging

import astropy.table
from astroquery.gaia import Gaia
from astroquery.utils.tap.core import TapPlus

from mirar.catalog.base.base_catalog import DEFAULT_SNR_THRESHOLD
from mirar.catalog.base.base_gaia import BaseGaia2Mass

logger = logging.getLogger(__name__)

# Disable astroquery warnings
logging.getLogger("astroquery").setLevel(logging.WARNING)

# URL for Gaia backup
ARI_URL = "https://gaia.ari.uni-heidelberg.de/tap"

gaia_ari = TapPlus(url=ARI_URL)


class Gaia2MassTAP(BaseGaia2Mass):
    """
    Crossmatched Gaia/2Mass catalog
    """

    def __init__(self, *args, snr_threshold: float = DEFAULT_SNR_THRESHOLD, **kwargs):
        super().__init__(*args, **kwargs)
        self.snr_threshold = snr_threshold

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
            f"SELECT * FROM gaiadr3.gaia_source AS g "
            f"JOIN gaiadr3.tmass_psc_xsc_best_neighbour AS xmatch USING (source_id) "
            f"JOIN gaiadr3.tmass_psc_xsc_join AS xjoin "
            f"  ON xmatch.original_ext_source_id = xjoin.original_psc_source_id "
            f"JOIN gaiadr1.tmass_original_valid AS tmass "
            f"  ON xjoin.original_psc_source_id = tmass.designation "
            f"WHERE CONTAINS(POINT('ICRS', g.ra, g.dec), "
            f"CIRCLE('ICRS', {ra_deg:.4f}, {dec_deg:.4f}, "
            f"{self.search_radius_arcmin / 60.:.4f}))=1 "
            f"AND tmass.{self.filter_name.lower()[0]}_m > {self.min_mag:.2f} "
            f"AND tmass.{self.filter_name.lower()[0]}_m < {self.max_mag:.2f} "
            f"AND tmass.{self.filter_name.lower()[0]}_msigcom < "
            f"{1.086 / self.snr_threshold:.3f} "
            f"AND xmatch.number_of_mates=0 "
            f"AND xmatch.number_of_neighbours=1"
            f";"
        )

        job = Gaia.launch_job_async(cmd, dump_to_file=False)
        src_list = job.get_results()

        src_list["k_m"] = src_list["ks_m"]
        src_list["k_msigcom"] = src_list["ks_msigcom"]
        src_list = self.convert_to_ab_mag(src_list)
        src_list["ph_qual"] = src_list["ph_qual"].astype(str)
        src_list["ra_errdeg"] = src_list["ra_error"] / 3.6e6
        src_list["dec_errdeg"] = src_list["dec_error"] / 3.6e6

        src_list["FLAGS"] = 0
        src_list["magnitude"] = src_list[f"{self.filter_name.lower()[0]}_m"]
        src_list["magnitude_err"] = src_list[f"{self.filter_name.lower()[0]}_msigcom"]

        logger.debug(f"Found {len(src_list)} sources in Gaia")
        return src_list


class Gaia2MassARI(BaseGaia2Mass):
    """
    Crossmatched Gaia/2Mass catalog
    """

    def __init__(self, *args, snr_threshold: float = DEFAULT_SNR_THRESHOLD, **kwargs):
        super().__init__(*args, **kwargs)
        self.snr_threshold = snr_threshold

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
            f"SELECT g.source_id, g.ra_error, g.dec_error, g.ra, g.dec, g.ruwe, "
            f"g.phot_rp_mean_mag, g.phot_bp_mean_mag, g.phot_g_mean_mag, "
            f"tbest.*, tmass.* FROM gaiadr3.gaia_source AS g, "
            f"gaiadr3.tmass_psc_xsc_best_neighbour AS tbest, "
            f"extcat.twomass AS tmass "
            f"WHERE g.source_id = tbest.source_id "
            f"AND tbest.original_ext_source_id = tmass.mainid "
            f"AND CONTAINS(POINT('ICRS', g.ra, g.dec), "
            f"CIRCLE('ICRS', {ra_deg:.4f}, {dec_deg:.4f}, "
            f"{self.search_radius_arcmin / 60.:.4f}))=1 "
            f"AND tmass.{self.filter_name}mag > {self.min_mag:.2f} "
            f"AND tmass.{self.filter_name}mag < {self.max_mag:.2f} "
            f"AND tmass.e_{self.filter_name}mag < {1.086 / self.snr_threshold: .3f} "
            f"AND tbest.number_of_mates=0 "
            f"AND tbest.number_of_neighbours=1;"
        )
        job = gaia_ari.launch_job_async(cmd)
        src_list = job.get_results()

        src_list.rename_column("jmag", "j_m")
        src_list.rename_column("hmag", "h_m")
        src_list.rename_column("kmag", "k_m")
        src_list.rename_column("e_jmag", "j_msigcom")
        src_list.rename_column("e_hmag", "h_msigcom")
        src_list.rename_column("e_kmag", "k_msigcom")

        src_list = self.convert_to_ab_mag(src_list)

        src_list["ph_qual"] = src_list["qflg"].astype(str)
        src_list["ra_errdeg"] = src_list["ra_error"] / 3.6e6
        src_list["dec_errdeg"] = src_list["dec_error"] / 3.6e6

        src_list["FLAGS"] = 0

        src_list["magnitude"] = src_list[f"{self.filter_name.lower()}_m"]
        src_list["magnitude_err"] = src_list[f"{self.filter_name.lower()}_msigcom"]

        logger.debug(f"Found {len(src_list)} sources in Gaia")
        return src_list
