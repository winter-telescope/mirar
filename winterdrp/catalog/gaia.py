import logging
import astropy.table
from astroquery.gaia import Gaia
from winterdrp.catalog.base_catalog import BaseCatalog

logger = logging.getLogger(__name__)


class Gaia2Mass(BaseCatalog):

    abbreviation = "tmass"

    def __init__(
            self,
            search_radius_arcmin: float,
            min_mag: float,
            max_mag: float,
            ph_qual_cut: bool = False
    ):
        super().__init__(search_radius_arcmin, min_mag, max_mag)
        self.ph_qual_cut = ph_qual_cut

    def get_catalog(
            self,
            ra_deg: float,
            dec_deg: float,
            search_radius_arcmin: float,
            min_mag: float,
            max_mag: float,
    ) -> astropy.table.Table:

        logger.info(
            f'Querying 2MASS - Gaia cross-match around RA {ra_deg:.4f}, '
            f'Dec {dec_deg:.4f} with a radius of {search_radius_arcmin:.4f} arcmin'
        )

        cmd = f"SELECT * FROM gaiadr2.gaia_source AS g, " \
              f"gaiadr2.tmass_best_neighbour AS tbest, " \
              f"gaiadr1.tmass_original_valid AS tmass " \
              f"WHERE g.source_id = tbest.source_id " \
              f"AND tbest.tmass_oid = tmass.tmass_oid " \
              f"AND CONTAINS(POINT('ICRS', g.ra, g.dec), " \
              f"CIRCLE('ICRS', {ra_deg:.4f}, {dec_deg:.4f}, {search_radius_arcmin / 60:.4f}))=1 " \
              f"AND tmass.j_m > {min_mag:.2f} " \
              f"AND tmass.j_m < {max_mag:.2f} " \
              f"AND tbest.number_of_mates=0 " \
              f"AND tbest.number_of_neighbours=1"

        if self.ph_qual_cut:
            cmd += f"AND tmass.ph_qual='AAA';"
        else:
            cmd += ";"

        job = Gaia.launch_job_async(
            cmd,
            dump_to_file=False
        )
        t = job.get_results()
        t['ph_qual'] = t['ph_qual'].astype(str)
        t['ra_errdeg'] = t['ra_error'] / 3.6e6
        t['dec_errdeg'] = t['dec_error'] / 3.6e6
        t['FLAGS'] = 0
        return t


# def make_gaia_catalog(
#         ra_deg: float,
#         dec_deg: float,
#         tm_cat_name: str,
#         write_ldac: bool = False
# ):
#
#     t = get_catalog(
#         ra_deg=ra_deg,
#         dec_deg=dec_deg,
#         search_radius_arcmin=mk.catalog_search_radius_arcmin,
#         min_mag=mk.catalog_min_mag,
#         max_mag=mk.catalog_max_mag,
#         ph_qual_cut=False
#     )
#
#     t['ph_qual'] = t['ph_qual'].astype(str)
#     t['ra_errdeg'] = t['ra_error'] / 3.6e6
#     t['dec_errdeg'] = t['dec_error'] / 3.6e6
#     t['FLAGS'] = 0
#
#     if write_ldac:
#
#         output_path = tm_cat_name + '.ldac'
#
#         if os.path.exists(output_path):
#             os.remove(output_path)
#
#         logger.info(f"Saving catalog to {output_path}")
#
#         aw.save_table_as_ldac(t, output_path)
