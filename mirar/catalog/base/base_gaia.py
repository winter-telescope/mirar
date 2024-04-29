"""
Module for obtaining a Gaia/2Mass catalog
"""

import logging
from abc import ABC
from typing import Optional

import astropy.table
import astropy.units as u
import numpy as np
from astropy.coordinates import SkyCoord

from mirar.catalog.base.base_catalog import BaseCatalog
from mirar.utils.ldac_tools import get_table_from_ldac

logger = logging.getLogger(__name__)

# 2MASS values from https://iopscience.iop.org/article/10.1086/376474
zeromag_2mass = {"j": 1594.0 * u.Jansky, "h": 1024.0 * u.Jansky, "k": 666.8 * u.Jansky}

offsets_2mass = {key: zm.to("mag(AB)").value for key, zm in zeromag_2mass.items()}


class BaseGaia2Mass(BaseCatalog, ABC):
    """
    Base Gaia/2Mass catalog
    """

    abbreviation = "tmass"

    def __init__(
        self,
        *args,
        filter_name: str = "j",
        snr_threshold: float = 5,
        trim: bool = False,
        image_catalog_path: Optional[str] = None,
        acceptable_j_ph_quals: str | list[str] = None,
        acceptable_h_ph_quals: str | list[str] = None,
        acceptable_k_ph_quals: str | list[str] = None,
        **kwargs,
    ):
        filter_name = filter_name.lower()

        super().__init__(*args, filter_name=filter_name, **kwargs)

        assert (
            filter_name in offsets_2mass.keys()
        ), f"Filter name must be one of {offsets_2mass.keys()}, not '{filter_name}'"

        self.trim = trim
        self.image_catalog_path = image_catalog_path
        self.snr_threshold = snr_threshold

        if isinstance(acceptable_j_ph_quals, str):
            acceptable_j_ph_quals = [acceptable_j_ph_quals]
        if isinstance(acceptable_h_ph_quals, str):
            acceptable_h_ph_quals = [acceptable_h_ph_quals]
        if isinstance(acceptable_k_ph_quals, str):
            acceptable_k_ph_quals = [acceptable_k_ph_quals]

        self.acceptable_ph_quals = {
            "j": acceptable_j_ph_quals,
            "h": acceptable_h_ph_quals,
            "k": acceptable_k_ph_quals,
        }

        if self.acceptable_ph_quals[self.filter_name] is None:
            self.acceptable_ph_quals[self.filter_name] = ["A"]

        for filt, val in self.acceptable_ph_quals.items():
            if val is None:
                self.acceptable_ph_quals[filt] = ["A", "B", "C"]

    def convert_to_ab_mag(self, src_list: astropy.table.Table) -> astropy.table.Table:
        """
        Convert 2MASS magnitudes to AB magnitudes

        :param src_list: Source list
        :return: Source list with AB magnitudes
        """
        # Convert to AB magnitudes
        for key in ["j", "h", "k"]:
            offset = offsets_2mass[self.filter_name.lower()]
            src_list[f"{key}_m"] += offset
            logger.debug(f"Adding {offset:.2f} to convert from 2MASS to AB magnitudes")
        return src_list

    def get_catalog(
        self,
        ra_deg: float,
        dec_deg: float,
    ) -> astropy.table.Table:

        src_list = self.get_source_table(ra_deg, dec_deg)

        j_phquals = [x[0] for x in src_list["ph_qual"]]
        h_phquals = [x[1] for x in src_list["ph_qual"]]
        k_phquals = [x[2] for x in src_list["ph_qual"]]

        j_phmask = np.array([x in self.acceptable_ph_quals["j"] for x in j_phquals])
        h_phmask = np.array([x in self.acceptable_ph_quals["h"] for x in h_phquals])
        k_phmask = np.array([x in self.acceptable_ph_quals["k"] for x in k_phquals])

        phmask = j_phmask & h_phmask & k_phmask

        src_list = src_list[phmask]
        src_list = src_list[src_list["magnitude_err"] < (1.086 / self.snr_threshold)]
        if self.trim:
            if self.image_catalog_path is None:
                err = (
                    "Gaia catalog trimming requested "
                    "but no sextractor catalog path specified."
                )
                logger.error(err)
                raise ValueError(err)

            image_catalog = get_table_from_ldac(self.image_catalog_path)
            src_list = self.trim_catalog(src_list, image_catalog)
            logger.debug(f"Trimmed to {len(src_list)} sources in Gaia")

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

    def get_source_table(self, ra_deg: float, dec_deg: float) -> astropy.table.Table:
        """
        Get the source table for a given position

        :param ra_deg: RA in degrees
        :param dec_deg: Dec in degrees
        :return: Table
        """
        raise NotImplementedError
