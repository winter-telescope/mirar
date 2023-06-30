"""
Module for downloading a reference image to perform calibration,
and extracting sources from that reference image
"""
import logging
import urllib.error
import urllib.request
from typing import Optional

import astropy.units as u
import numpy as np
from astropy.coordinates import SkyCoord
from astroquery.vizier import Vizier

from mirar.processors.astrometry.autoastrometry.errors import (
    AstrometryReferenceError,
    AstrometryURLError,
)
from mirar.processors.astrometry.autoastrometry.sources import BaseSource, compare_mag
from mirar.processors.astrometry.autoastrometry.utils import dec_str_2_deg, ra_str_2_deg

logger = logging.getLogger(__name__)


def get_catalog(
    catalog: str,
    ra: float,
    dec: float,
    box_size_arcsec: float,
    min_mag: float = 8.0,
    max_mag: Optional[float] = None,
    max_pm: float = 60.0,
) -> list[BaseSource]:
    """
    Get a reference catalog around ra/dec, with radius

    :param catalog: Catalog to use ('ub2'/'sdss'/'tmc')
    :param ra: ra
    :param dec: declination
    :param box_size_arcsec: box to use
    :param min_mag: min mag of sources
    :param max_mag: max mag of sources
    :param max_pm: max pm
    :return: list of reference sources
    """
    # Get catalog from USNO

    if max_mag is None:
        if catalog == "ub2":
            max_mag = 21.0  # 19.5
        elif catalog == "sdss":
            max_mag = 22.0
        elif catalog == "tmc":
            max_mag = 20.0
        else:
            err = "Catalog not recognised. Please select 'ub2', 'sdss' or tmc'."
            logger.error(err)
            raise ValueError(err)

    ra_col = 1
    dec_col = 2

    if catalog == "tmc":
        mag_col = 3
    else:
        mag_col = 6

    pm_ra_col = 10
    pm_dec_col = 11
    query_url = (
        f"http://tdc-www.harvard.edu/cgi-bin/scat?catalog={catalog}"
        f"&ra={ra}&dec={dec}&system=J2000&rad={-box_size_arcsec}"
        f"&sort=mag&epoch=2000.00000&nstar=6400"
    )

    with urllib.request.urlopen(query_url) as cat:
        cat_lines = cat.readlines()

    if len(cat_lines) > 6400 - 20:
        logger.warning(
            "Reached maximum catalog query size. Gaps may be "
            "present in the catalog, leading to a poor solution "
            "or no solution. Decrease the search radius."
        )

    cat_list = []

    for line in cat_lines:
        inline = line.strip()

        if len(inline) <= 2:
            continue

        if inline[0:2] == "#:":
            inline_arg = inline[2:].split(",")
            ra_col = int(inline_arg[0]) - 1
            dec_col = int(inline_arg[1]) - 1
            if len(inline_arg) > 2:
                mag_col = int(inline_arg[2]) - 1
            continue

        if (int(inline[0]) < ord("0") or int(inline[0]) > ord("9")) and str(
            inline[0]
        ) != ".":
            continue  # this may be too overzealous about
        if (int(inline[1]) < ord("0") or int(inline[1]) > ord("9")) and str(
            inline[1]
        ) != ".":
            continue  # removing comments...

        inline_arg_byte = inline.split()
        inline_arg = [str(bytes(a), "utf-8") for a in inline_arg_byte]
        n_arg = len(inline_arg)

        if inline_arg[ra_col].find(":") == -1:
            ra = float(inline_arg[ra_col])
        else:
            ra = ra_str_2_deg(inline_arg[ra_col])

        if inline_arg[dec_col].find(":") == -1:
            dec = float(inline_arg[dec_col])
        else:
            dec = dec_str_2_deg(inline_arg[dec_col])

        if n_arg > mag_col >= 0:
            try:
                mag = float(inline_arg[mag_col])
            except ValueError:
                mag = float(inline_arg[mag_col][0:-2])
        else:
            mag = max_mag

        # if user_cat is False and n_arg > pm_ra_col and n_arg > pm_dec_col:
        if n_arg > pm_ra_col and n_arg > pm_dec_col:
            pm_ra = float(inline_arg[pm_ra_col])
            pm_dec = float(inline_arg[pm_dec_col])
        else:
            pm_ra = pm_dec = 0

        if mag > max_mag:
            continue
        if mag < min_mag:
            continue

        if abs(pm_ra) > max_pm or abs(pm_dec) > max_pm:
            continue

        source = BaseSource(ra, dec, mag)  # process the line into an object
        cat_list.append(source)

    cat_list.sort(key=compare_mag)

    return cat_list


def get_catalog_astroquery(
    catalog: str,
    ra: float,
    dec: float,
    box_size_arcsec: float,
    min_mag: float = 8.0,
    max_mag: Optional[float] = None,
    max_pm: float = 60.0,
) -> list[BaseSource]:
    """
    Get a reference catalogue using astroquery, around ra/dec

    :param catalog: Catalog to use
    :param ra: ra
    :param dec: dec
    :param box_size_arcsec: radius of box
    :param min_mag: min mag of stars
    :param max_mag: max mag of stars
    :param max_pm: max pm
    :return: list of reference stars
    """

    ra_col_key = dec_col_key = mag_col_key = catalog_str = pm_ra_key = pm_dec_key = ""

    if catalog == "sdss":
        ra_col_key = "RA_ICRS"
        dec_col_key = "DE_ICRS"
        mag_col_key = "rmag"
        pm_ra_key = "pmRA"
        pm_dec_key = "pmDE"
        max_mag = 22
        catalog_str = "V/147/sdss12"

    if catalog == "usno":
        ra_col_key = "RAJ2000"
        dec_col_key = "DEJ2000"
        mag_col_key = "R1mag"
        pm_ra_key = "pmRA"
        pm_dec_key = "pmDE"
        max_mag = 21
        catalog_str = "I/284/out"

    if catalog == "tmc":
        ra_col_key = "RAJ2000"
        dec_col_key = "DEJ2000"
        mag_col_key = "Jmag"
        pm_ra_key = "pmRA"
        pm_dec_key = "pmDE"
        max_mag = 17
        catalog_str = (
            "I/317/sample"  # using the PPMXL catalog of 2MASS + proper motions
        )

    crd = SkyCoord(ra=ra, dec=dec, unit=(u.deg, u.deg))
    vizier_cat = Vizier(
        columns=[ra_col_key, dec_col_key, mag_col_key, pm_ra_key, pm_dec_key],
        column_filters={f"{mag_col_key}": f"<{max_mag}"},
    )
    vizier_cat.ROW_LIMIT = -1
    logger.info(
        f"Querying {catalog} around {ra},{dec} and a radius "
        f"{int(box_size_arcsec / 60)} arcminutes "
        f"using - {vizier_cat.columns}, {vizier_cat.column_filters}"
    )
    # pylint: disable=no-member
    result = vizier_cat.query_region(
        crd, width=f"{int(box_size_arcsec / 60)}m", catalog=catalog_str
    )
    if len(result) == 0:
        return []
    table = result[0]

    n_cat = len(table)
    cat_density = n_cat / (2 * box_size_arcsec / 60.0) ** 2
    logger.debug(f"{n_cat} good catalog objects.")
    logger.debug(f"Source density of {cat_density} /arcmin^2")

    mask = (
        (table[mag_col_key] < max_mag)
        & (table[mag_col_key] > min_mag)
        & (np.abs(table[pm_ra_key]) < max_pm)
        & (np.abs(table[pm_dec_key]) < max_pm)
    )
    cat = table[mask]
    cat_list = []

    for src in cat:
        cat_list.append(BaseSource(src[ra_col_key], src[dec_col_key], src[mag_col_key]))

    cat_list.sort(key=compare_mag)
    return cat_list


def get_ref_sources_from_catalog_astroquery(
    catalog: str, center_ra: float, center_dec: float, box_size_arcsec: float
) -> tuple[list[BaseSource], int, float]:
    """
    Get reference sources from an astropquery catalogue

    :param catalog: catalog to use
    :param center_ra: ra
    :param center_dec: dec
    :param box_size_arcsec: radius
    :return: ref catalog, n_cat, cat_density
    """
    ref_src_list = []
    if catalog is None:
        try:
            trycats = ["sdss", "usno", "tmc"]
            for trycat in trycats:
                ref_src_list = get_catalog_astroquery(
                    catalog=trycat,
                    ra=center_ra,
                    dec=center_dec,
                    box_size_arcsec=box_size_arcsec,
                )
                if len(ref_src_list) > 15:
                    break

        except urllib.error.URLError as exc:
            err = "No catalog is available.  Check your internet connection."
            logger.error(err)
            raise AstrometryURLError(err) from exc

    n_cat = len(ref_src_list)
    cat_density = n_cat / (2 * box_size_arcsec / 60.0) ** 2
    logger.debug(f"{n_cat} good catalog objects.")
    logger.debug(f"Source density of {cat_density} /arcmin^2")

    if n_cat == 0:
        err = (
            "No objects found in catalog. "
            "The web query failed, all stars were excluded "
            "by the FHWM clip, or the image is too small.  "
            "Check input parameters or your internet connection."
        )
        logger.error(err)
        raise AstrometryReferenceError(err)

    if 0 < n_cat < 5:
        err = (
            f"Only {n_cat} catalog objects in the search zone. "
            f"Increase the magnitude threshold or box size."
        )
        logger.error(err)
        raise AstrometryReferenceError(err)

    return ref_src_list, n_cat, cat_density


def get_ref_sources_from_catalog(
    catalog: str,
    center_ra: float,
    center_dec: float,
    box_size_arcsec: float,
) -> tuple[list[BaseSource], int, float]:
    """
    Get reference sources from a catalogue

    :param catalog: catalog to use
    :param center_ra: ra
    :param center_dec: dec
    :param box_size_arcsec: radius
    :return: ref catalog, n_cat, cat_density
    """
    # If no catalog specified, check availability of SDSS
    ref_src_list, n_cat, cat_density = [], 0, 0
    logger.debug(f"catalog is {catalog}")
    if catalog is None:
        trycats = ["ub2", "tmc", "sdss"]
        logger.debug(f"No catalog specified, so defaulting to {trycats}")
    else:
        trycats = [catalog]

    try:
        for trycat in trycats:
            testqueryurl = (
                f"http://tdc-www.harvard.edu/cgi-bin/"
                f"scat?catalog={trycat}&ra={center_ra}"
                f"&dec={center_dec}&system=J2000&rad=-90"
            )
            logger.debug(f"Trying {testqueryurl}")
            with urllib.request.urlopen(testqueryurl, timeout=30) as check:
                checklines = check.readlines()
            logger.debug(f"Found {len(checklines)}")
            if len(checklines) > 0:
                catalog = trycat
                logger.debug(f"Using catalog {catalog}")

                ref_src_list = get_catalog(
                    catalog=catalog,
                    ra=center_ra,
                    dec=center_dec,
                    box_size_arcsec=box_size_arcsec,
                )

                n_cat = len(ref_src_list)
                cat_density = n_cat / (2 * box_size_arcsec / 60.0) ** 2
                logger.debug(f"{n_cat} good catalog objects.")
                logger.debug(f"Source density of {cat_density} /arcmin^2")

                if n_cat > 5:
                    break

    except urllib.error.URLError as exc:
        err = "No catalog is available.  Check your internet connection."
        logger.error(err)
        raise AstrometryURLError(err) from exc

    if n_cat == 0:
        err = (
            "No objects found in catalog. The web query failed, "
            "all stars were excluded by the FHWM clip, or the image is too small.  "
            "Check input parameters or your internet connection."
        )
        logger.error(err)
        raise AstrometryReferenceError(err)

    if 0 < n_cat < 5:
        err = (
            f"Only {n_cat} catalog objects in the search zone. "
            f"Increase the magnitude threshold or box size."
        )
        logger.error(err)
        raise AstrometryReferenceError(err)
    # Load in reference star catalog
    logger.debug(f"Number of catalog objects found is {n_cat}")
    return ref_src_list, n_cat, cat_density
