#!/bin/env python
"""
autoastrometry.py - a fast astrometric solver

author: Daniel Perley (dperley@astro.caltech.edu)
last significant modifications 2012-04-23

4/23: program can actually be overwhelmed by too many good matches (too high maxrad).
need to fix this.

Modified by Kishalay De (kde@astro.caltech.edu) for removing dependency
on deprecated pyfits and making it compatible
with astropy headers and python 3.6 (June 11, 2018)

Modified/dissected by Robert Stein (rdstein@caltech.edu) in 2021/2022
for incorporation into mirar
Also converted to python 3.10
"""

import logging
import math
import os
import sys
from pathlib import Path
from typing import Optional

import numpy as np
from astropy.io import fits

from mirar.paths import base_output_dir
from mirar.processors.astromatic.sextractor.settings import (
    write_config_file,
    write_param_file,
)
from mirar.processors.astromatic.sextractor.sourceextractor import DEFAULT_SATURATION
from mirar.processors.astrometry.autoastrometry.crossmatch import (
    crosscheck_source_lists,
    distance_match,
)
from mirar.processors.astrometry.autoastrometry.detect import (
    DEFAULT_MAX_FWHM,
    DEFAULT_MIN_FWHM,
    get_img_src_list,
)
from mirar.processors.astrometry.autoastrometry.errors import (
    AstrometryCrossmatchError,
    AstrometryReferenceError,
    AstrometrySourceError,
    AstrometryURLError,
)
from mirar.processors.astrometry.autoastrometry.io import (
    export_src_lists,
    parse_header,
    write_region_file,
    write_text_file,
)
from mirar.processors.astrometry.autoastrometry.reference import (
    get_ref_sources_from_catalog,
    get_ref_sources_from_catalog_astroquery,
)
from mirar.processors.astrometry.autoastrometry.utils import median, stdev

logger = logging.getLogger(__name__)

# these defaults should generally not be altered.
DEFAULT_TOLERANCE = 0.01
DEFAULT_PA_TOLERANCE = 1.4


############################################
def autoastrometry(
    filename: str | Path,
    pixel_scale: Optional[float] = None,
    pa: Optional[float] = None,
    inv: bool = False,
    unc_pa: Optional[float] = None,
    user_ra_deg: Optional[float] = None,
    user_dec_deg: Optional[float] = None,
    max_ellip: float = 0.5,
    box_size_arcsec: Optional[float] = None,
    max_rad: Optional[float] = None,
    tolerance: float = DEFAULT_TOLERANCE,
    catalog: Optional[str] = None,
    overwrite: bool = True,
    outfile: str = "",
    output_dir: str = base_output_dir,
    temp_file: Optional[str] = None,
    saturation: float = DEFAULT_SATURATION,
    no_rot: bool = False,
    min_fwhm: float = DEFAULT_MIN_FWHM,
    max_fwhm: float = DEFAULT_MAX_FWHM,
    write_crosscheck_files: bool = False,
):
    """

    Parameters
    ----------
    filename: Path of file
    pixel_scale: The pixel scale in arcsec/pix.  Must be within ~1%. By default: ???
    pa: The position angle in degrees.  Not usually needed.
    unc_pa: Uncertainty of the position angle (degrees)
    inv: Reverse(=positive) parity.
    user_ra_deg: RA in deg
    user_dec_deg: Dec in deg
    max_ellip: Maximum elliptical something?
    box_size_arcsec: Half-width of box for reference catalog query (arcsec)
    max_rad: Maximum distance to look for star pairs.
    tolerance: Amount of slack allowed in match determination
    catalog: Catalog to use (ub2, tmc, sdss, or file)
    overwrite: Overwrite output files
    outfile: Output file
    output_dir: Directory for output file
    temp_file: Temporary file
    saturation: Saturation level; do not use stars exceeding.
    no_rot: Some kind of bool
    min_fwhm: Minimum fwhm
    max_fwhm: Maximum fwhm
    write_crosscheck_files: Bool for whether to write region and other crosscheck files

    Returns
    -------

    """

    if temp_file is None:
        temp_file = f"temp_{os.path.basename(filename)}"

    # temp_path = os.path.join(os.path.dirname(filename), temp_file)
    temp_path = os.path.join(output_dir, temp_file)

    if overwrite:
        temp_path = filename

    base_output_path = os.path.join(
        output_dir, ".".join(os.path.basename(filename).split(".")[:-1])
    )

    # Block A
    # Parse header

    nx_pix, ny_pix, cd11, cd12, cd21, cd22, crpix1, crpix2, cra, cdec = parse_header(
        file_path=filename,
        temp_path=temp_path,
        pixel_scale=pixel_scale,
        pa=pa,
        inv=inv,
        user_ra_deg=user_ra_deg,
        user_dec_deg=user_dec_deg,
    )

    # Block B
    # Sextract stars to produce image star catalog

    img_src_list = get_img_src_list(
        img_path=temp_path,
        nx_pix=nx_pix,
        ny_pix=ny_pix,
        border=3,
        corner=12,
        min_fwhm=min_fwhm,
        max_fwhm=max_fwhm,
        max_ellip=max_ellip,
        saturation=saturation,
        base_output_path=base_output_path,
        write_crosscheck_files=write_crosscheck_files,
    )

    # Block C

    sci_ext = 0

    with fits.open(temp_path) as hdu:
        header = hdu[sci_ext].header

    # Get image info from header (even if we put it there in the first place)
    if cd11 * cd22 < 0 or cd12 * cd21 > 0:
        parity = -1
    else:
        parity = 1

    x_scale = math.sqrt(cd11**2 + cd21**2)
    y_scale = math.sqrt(cd12**2 + cd22**2)
    init_pa = -parity * np.arctan2(cd21 * y_scale, cd22 * x_scale) * 180 / math.pi
    x_scale = abs(x_scale)
    y_scale = abs(y_scale)
    field_width = max(x_scale * nx_pix, y_scale * ny_pix) * 3600.0
    area_sq_deg = x_scale * nx_pix * y_scale * ny_pix
    area_sq_min = area_sq_deg * 3600.0
    center_x = nx_pix / 2
    center_y = ny_pix / 2
    center_dx = center_x - crpix1
    center_dy = center_y - crpix2
    center_ra = (
        cra
        - center_dx * x_scale * math.cos(init_pa * math.pi / 180.0)
        + center_dy * y_scale * math.sin(init_pa * math.pi / 180.0)
    )
    center_dec = (
        cdec
        + parity * center_dx * x_scale * math.sin(-init_pa * math.pi / 180.0)
        + center_dy * y_scale * math.cos(init_pa * math.pi / 180.0)
    )

    if box_size_arcsec is None:
        box_size_arcsec = field_width

    # this has only been checked for a PA of zero.

    logger.debug(
        f"Initial WCS info: \n"
        f'   pixel scale:     x={x_scale * 3600:.4f}"/pix,   '
        f'y={y_scale * 3600:.4f}"/pix \n'
        f"   position angle: PA={init_pa:.2f}"
    )

    if parity == 1:
        logger.debug("   normal parity")
    if parity == -1:
        logger.debug("   inverse parity")

    logger.debug(f"   center:        RA={center_ra:.6f}, dec={center_dec:.6f}")

    n_img = len(img_src_list)

    if n_img < 4:
        err = (
            f"Only {n_img} good stars were found in the image.  "
            f"The image is too small or shallow, "
            f"the detection threshold is set too high, "
            f"or stars and cosmic rays are being confused."
        )
        logger.error(err)
        write_text_file("det.init.txt", img_src_list)
        write_region_file("det.im.reg", img_src_list, "red", "img")
        raise AstrometrySourceError(err)

    img_density = len(img_src_list) / area_sq_min
    logger.debug(f"Source img_density of {img_density:.4f} /arcmin^2")

    # Block D

    try:
        ref_src_list, n_ref, ref_density = get_ref_sources_from_catalog(
            catalog=catalog,
            center_ra=center_ra,
            center_dec=center_dec,
            box_size_arcsec=box_size_arcsec,
        )
    except (
        TimeoutError,
        AstrometryURLError,
        AstrometrySourceError,
        AstrometryReferenceError,
    ):
        ref_src_list, n_ref, ref_density = get_ref_sources_from_catalog_astroquery(
            catalog=catalog,
            center_ra=center_ra,
            center_dec=center_dec,
            box_size_arcsec=box_size_arcsec,
        )
    # Block E

    (
        img_src_list,
        n_img,
        img_density,
        ref_src_list,
        n_ref,
        ref_density,
    ) = crosscheck_source_lists(
        img_src_list=img_src_list,
        n_img=n_img,
        img_density=img_density,
        ref_src_list=ref_src_list,
        n_ref=n_ref,
        ref_density=ref_density,
        box_size_arcsec=box_size_arcsec,
        area_sq_min=area_sq_min,
    )

    # Block F

    if write_crosscheck_files:
        export_src_lists(
            img_src_list=img_src_list,
            ref_src_list=ref_src_list,
            base_output_path=base_output_path,
        )

    # The catalogs have now been completed.

    # Block G

    # Now start getting into the actual astrometry.

    min_rad = 5.0
    if max_rad is None:
        # numcomp ~ 15 [look at 15 closest objects in sparse dataset]
        max_rad = 60 * (15 / (math.pi * min(img_density, ref_density))) ** 0.5
        max_rad = max(max_rad, 60.0)
        if max_rad == 60.0:
            min_rad = (
                10.0  # in theory could scale this up further to reduce #comparisons
            )
        max_rad = min(max_rad, field_width * 3.0 / 4)

        # note that img_density is per arcmin^2, while the radii are in arcsec,
        # hence the conversion factor.

    circ_density = img_density * min(
        [area_sq_min, (math.pi * (max_rad / 60.0) ** 2 - math.pi * (min_rad / 60) ** 2)]
    )
    circ_ref_density = ref_density * (
        math.pi * (max_rad / 60.0) ** 2 - math.pi * (min_rad / 60) ** 2
    )
    n_ref_src_per_image = ref_density * area_sq_min

    logger.debug(
        "After trimming: "
        f"{len(img_src_list)} detected objects "
        f"({img_density:.2f}/arcmin^2, "
        f"{circ_density:.1f}/searchzone)"
    )
    logger.debug(
        f"{len(ref_src_list)} catalog objects "
        f"({ref_density:.2f}/arcmin^2, "
        f"{circ_ref_density:.1f}/searchzone)"
    )

    pa_tolerance = DEFAULT_PA_TOLERANCE

    # RS: WTF is x**1 for?
    expect_false_trios = (
        n_img
        * n_ref
        * circ_density**2
        * circ_ref_density**2
        * tolerance**2
        * (pa_tolerance / 360.0) ** 1
    )

    # fraction of stars in image that are also in catalog - a guess
    overlap_first_guess = 0.3 * min(1.0, ref_density / img_density)
    # but how many matches >3 and >4?  some annoying binomial thing
    true_matches_per_star = circ_density * overlap_first_guess

    req_matches = 3
    if expect_false_trios > 30 and true_matches_per_star >= 4:
        req_matches = 4
    # should check that this will actually work for the catalog, too.
    if n_ref_src_per_image <= 6 or n_img <= 6:
        req_matches = 2
    if n_ref_src_per_image <= 3 or n_img <= 3:
        req_matches = 1
    # for an extremely small or shallow image

    logger.debug(f'Pair comparison search radius: {max_rad:.2f}"')
    logger.debug(f"Using req_matches = {req_matches}")

    primary_match_img, primary_match_ref, mpa = distance_match(
        img_src_list=img_src_list,
        ref_src_list=ref_src_list,
        base_output_path=base_output_path,
        max_rad=max_rad,
        min_rad=min_rad,
        tolerance=tolerance,
        req_match=req_matches,
        pa_tolerance=pa_tolerance,
        unc_pa=unc_pa,
        write_crosscheck_files=write_crosscheck_files,
    )

    n_match = len(primary_match_img)
    if n_match == 0:
        err = (
            " No valid matches found!\n "
            "Possible issues: \n"
            "  - The specified pixel scale (or PA or parity) is incorrect.  "
            "Double-check the input value. \n"
            "  - The field is outside the catalog search region.  "
            "Check header RA/DEC or increase search radius. \n"
            " - The routine is flooded by bad sources.  "
            "Specify or check the input seeing. \n"
            "  - The routine is flagging many real stars.  Check the input seeing. \n"
            " You can display a list of detected/catalog sources "
            "using det.im.reg and cat.wcs.reg. \n"
        )
        logger.error(err)
        raise AstrometryCrossmatchError(err)

    if n_match <= 2:
        logger.warning(
            f"Warning: only {n_match} match(es).  Astrometry may be unreliable."
        )
        logger.warning("   Check the pixel scale and parity and consider re-running.")

    # We now have the PA and a list of stars that are almost certain matches.
    median_pa = median(mpa)  # get average PA from the excellent values
    stdev_pa = stdev(mpa)

    sky_offset_pa = -parity * median_pa
    # This appears to be necessary for the printed value to agree with our
    # normal definition.

    logger.debug("PA offset:")
    logger.debug(f"  dPA = {sky_offset_pa:.3f}  (unc. {stdev_pa:.3f})")

    if no_rot <= 0:
        # Rotate the image to the new, correct PA
        #  NOTE: when CRPIX don't match CRVAL this shifts the center
        #  and screws things up.
        #  I don't understand why they don't always match.
        #  [[I think this was an equinox issue.
        #  should be solved now, but be alert for further problems.]]

        # Rotate....
        rot = median_pa * math.pi / 180
        # ...the image itself
        header["CD1_1"] = math.cos(rot) * cd11 - math.sin(rot) * cd21
        header["CD1_2"] = (
            math.cos(rot) * cd12 - math.sin(rot) * cd22
        )  # a parity issue may be involved here?
        header["CD2_1"] = math.sin(rot) * cd11 + math.cos(rot) * cd21
        header["CD2_2"] = math.sin(rot) * cd12 + math.cos(rot) * cd22
        # ...the coordinates (so we don't have to resex)
        for img_src in img_src_list:
            # do all of them, though this is not necessary
            img_src.rotate(median_pa, cra, cdec)

    else:
        if abs(sky_offset_pa) > 1.0:
            logger.warning(" (WARNING: image appears rotated, may produce bad shift)")
        logger.debug("  Skipping rotation correction ")

    im_ra_offset = []
    im_dec_offset = []
    for i, x in enumerate(primary_match_img):
        im_ra_offset.append(
            img_src_list[x].ra_deg - ref_src_list[primary_match_ref[i]].ra_deg
        )
        im_dec_offset.append(
            img_src_list[x].dec_deg - ref_src_list[primary_match_ref[i]].dec_deg
        )

    ra_offset = -median(im_ra_offset)
    dec_offset = -median(im_dec_offset)
    ra_std = stdev(im_ra_offset) * math.cos(
        cdec * math.pi / 180
    )  # all of these are in degrees
    dec_std = stdev(im_dec_offset)
    std_offset = math.sqrt(ra_std**2 + dec_std**2)

    ra_offset_arcsec = ra_offset * 3600 * math.cos(cdec * math.pi / 180)
    dec_offset_arcsec = dec_offset * 3600
    tot_offset_arcsec = (ra_offset_arcsec**2 + dec_offset**2) ** 0.5
    std_offset_arcsec = std_offset * 3600

    logger.debug("Spatial offset:")

    msg = (
        f'  dra = {ra_offset_arcsec:.2f}",'
        f'  ddec = {dec_offset_arcsec:.2f}"'
        f'  (unc. {std_offset_arcsec:.3f}")'
    )
    logger.debug(msg)

    if std_offset * 3600 > 1.0:
        logger.debug(
            "WARNING: poor solution - some matches may be bad.  Check pixel scale?"
        )

    header["CRVAL1"] = cra + ra_offset
    header["CRVAL2"] = cdec + dec_offset

    logger.debug(f'Updated header {header["CRVAL1"], header["CRVAL2"]}')
    try:
        oldcat = header["ASTR_CAT"]
        header["OLD_CAT"] = (oldcat, "Earlier reference catalog")
    except KeyError:
        pass

    header["ASTR_CAT"] = (catalog, "Reference catalog for autoastrometry")
    header["ASTR_UNC"] = (std_offset_arcsec, "Astrometric scatter vs. catalog (arcsec)")
    header["ASTR_SPA"] = (stdev_pa, "Measured uncertainty in PA (degrees)")
    header["ASTR_DPA"] = (sky_offset_pa, "Change in PA (degrees)")
    header["ASTR_OFF"] = (tot_offset_arcsec, "Change in center position (arcsec)")
    header["ASTR_NUM"] = (len(primary_match_img), "Number of matches")

    if write_crosscheck_files:
        write_text_file(
            file_path=os.path.splitext(base_output_path)[0] + ".det.wcs.txt",
            src_list=img_src_list,
        )

        # Write out a match list to allow doing a formal fit with WCStools.

        match_list_path = os.path.splitext(base_output_path)[0] + ".match.list"

        logger.info(f"Writing match list to {match_list_path}")

        with open(match_list_path, "w", encoding="utf8") as outmatch:
            for i, src_idx in enumerate(primary_match_img):
                ref_idx = primary_match_ref[i]
                outmatch.write(
                    f"{img_src_list[src_idx].x} {img_src_list[src_idx].y} "
                    f"{ref_src_list[ref_idx].ra_deg} {ref_src_list[ref_idx].dec_deg}\n"
                )

    logger.debug(f"Finished deriving astrometry for {filename}")

    # Could repeat with scale adjustment
    # Could then go back to full good catalog and match all sources

    if overwrite:
        outfile = filename
    elif outfile == "":
        slash_pos = filename.rfind("/")
        dir_name = filename[0 : slash_pos + 1]
        fil = filename[slash_pos + 1 :]
        # alternate behavior would always output to current directory
        outfile = f"{dir_name}a{fil}"

    if outfile is not None:
        with fits.open(temp_path) as hdu:
            hdu[sci_ext].header = header

            hdu.writeto(outfile, output_verify="silentfix", overwrite=True)
            logger.info(f"Written updated file to {outfile}")
        logger.debug(
            f"Derived center coordinates of {header['CRVAL1']}, {header['CRVAL2']}."
        )

    return (
        n_match,
        sky_offset_pa,
        stdev_pa,
        ra_offset_arcsec,
        dec_offset_arcsec,
        std_offset_arcsec,
    )


def run_autoastrometry_single(
    img_path: str | Path,
    seeing: Optional[float] = None,
    pixel_scale: Optional[float] = None,
    pa: Optional[float] = None,
    uncpa: Optional[float] = None,
    inv: bool = False,
    user_ra: Optional[float] = None,
    user_dec: Optional[float] = None,
    max_ellip: float = 0.5,
    box_size: Optional[float] = None,
    max_rad: Optional[float] = None,
    tolerance: float = DEFAULT_TOLERANCE,
    catalog: Optional[str] = None,
    overwrite: bool = False,
    outfile: Optional[str] = None,
    output_dir: str = base_output_dir,
    saturation: float = DEFAULT_SATURATION,
    no_rot: bool = False,
    write_crosscheck_files: bool = False,
):
    """Function based on 'autoastrometry.py' by Daniel Perley and Kishalay De.

    This function runs sextractor and automatically performs astrometry.

    Supplying the correct pixel scale (within 1%) and correct parity is critical
    if the image does not already contain this information in the FITS header.
    If you have difficulty solving a field correctly, double-check these values.
    If still having trouble, try opening temp.fits and an archival image of the
    field (from DSS, etc.) and loading the .reg files in DS9. The problem might
    be in the telescope pointing/header info (in this case, increase the boxsize)
    or good matching stars may be thrown away or confused by artifacts (in this
    case, specify a seeing value).  If the PA is known, restricting it can also
    help (try -upa 0.5); by default all orientations are searched.

    Catalog info:
    Leave the catalog field blank will use SDSS if available and USNO otherwise.
    The catalog query uses wcstools (tdc-www.harvard.edu/wcstools).  However, you
    can also use your own catalog file on disk if you prefer using -c [img_path]
    The default format is a text file with the first three columns indicating
    ra, dec, magnitude.  However, you can change the order of the columns by
    adding, e.g.#:1,2,6 to the first line.
    In this case, this would indicate that the RA is in the
    1st column, dec in the 2nd, and magnitude in the 6th.   The mag column can be
    omitted completely, although if the catalog is not the same depth as the
    image this may compromise the search results.

    Parameters
    ----------
    img_path: File to reduce
    seeing: Approximate seeing in pixels for CR/star/galaxy ID'ing.
    pixel_scale: The pixel scale in arcsec/pix.  Must be within ~1%. By default: ???
    pa: The position angle in degrees.  Not usually needed.
    uncpa: Uncertainty of the position angle (degrees)
    inv: Reverse(=positive) parity.
    user_ra: RA in deg
    user_dec: Dec in deg
    max_ellip: Maximum elliptical something?
    box_size: Half-width of box for reference catalog query (arcsec)
    max_rad: Maximum distance to look for star pairs.
    tolerance: Amount of slack allowed in match determination
    catalog: Catalog to use (ub2, tmc, sdss, or file)
    overwrite: Overwrite output files
    outfile: Output file
    output_dir: Output directory
    saturation: Saturation level; do not use stars exceeding.
    no_rot: Some kind of bool involving rotation
    write_crosscheck_files: Bool for whether to write region and other crosscheck files

    Returns
    -------
    """
    if seeing is None:
        min_fwhm = DEFAULT_MIN_FWHM  # 1.5
        max_fwhm = DEFAULT_MAX_FWHM  # 40
    else:
        min_fwhm = 0.7 * seeing
        max_fwhm = 2.0 * seeing

        write_param_file()

    write_config_file()
    logger.debug(f"Outfile is {outfile}")
    fit_info = autoastrometry(
        filename=img_path,
        pixel_scale=pixel_scale,
        pa=pa,
        inv=inv,
        unc_pa=uncpa,
        min_fwhm=min_fwhm,
        max_fwhm=max_fwhm,
        max_ellip=max_ellip,
        box_size_arcsec=box_size,
        max_rad=max_rad,
        user_ra_deg=user_ra,
        user_dec_deg=user_dec,
        tolerance=tolerance,
        catalog=catalog,
        overwrite=overwrite,
        outfile=outfile,
        output_dir=output_dir,
        saturation=saturation,
        no_rot=no_rot,
        write_crosscheck_files=write_crosscheck_files,
    )
    return fit_info


######################################################################


def run_autoastrometry_batch(
    files: str | Path | list[str | Path],
    seeing: Optional[float] = None,
    pixel_scale: Optional[float] = None,
    pa: Optional[float] = None,
    uncpa: Optional[float] = None,
    inv: bool = False,
    user_ra: Optional[float] = None,
    user_dec: Optional[float] = None,
    max_ellip: float = 0.5,
    box_size: Optional[float] = None,
    max_rad: Optional[float] = None,
    tolerance: float = DEFAULT_TOLERANCE,
    catalog: Optional[str] = None,
    overwrite: bool = False,
    outfile: Optional[str] = None,
    output_dir: str = base_output_dir,
    saturation: float = DEFAULT_SATURATION,
    no_rot: bool = False,
    write_crosscheck_files: bool = False,
):
    """Function based on 'autoastrometry.py' by Daniel Perley and Kishalay De.

    This function runs sextractor and automatically performs astrometry.

    Supplying the correct pixel scale (within 1%) and correct parity is critical
    if the image does not already contain this information in the FITS header.
    If you have difficulty solving a field correctly, double-check these values.
    If still having trouble, try opening temp.fits and an archival image of the
    field (from DSS, etc.) and loading the .reg files in DS9. The problem might
    be in the telescope pointing/header info (in this case, increase the boxsize)
    or good matching stars may be thrown away or confused by artifacts (in this
    case, specify a seeing value).  If the PA is known, restricting it can also
    help (try -upa 0.5); by default all orientations are searched.

    Catalog info:
    Leave the catalog field blank will use SDSS if available and USNO otherwise.
    The catalog query uses wcstools (tdc-www.harvard.edu/wcstools).  However, you
    can also use your own catalog file on disk if you prefer using -c [img_path]
    The default format is a text file with the first three columns indicating
    ra, dec, magnitude.  However, you can change the order of the columns by
    adding, e.g.#:1,2,6 to the first line.
    In this case, this would indicate that the RA is in the
    1st column, dec in the 2nd, and magnitude in the 6th.   The mag column can be
    omitted completely, although if the catalog is not the same depth as the
    image this may compromise the search results.

    Parameters
    ----------
    files: Files to reduce
    seeing: Approximate seeing in pixels for CR/star/galaxy ID'ing.
    pixel_scale: The pixel scale in arcsec/pix.  Must be within ~1%. By default: ???
    pa: The position angle in degrees.  Not usually needed.
    uncpa: Uncertainty of the position angle (degrees)
    inv: Reverse(=positive) parity.
    user_ra: RA in deg
    user_dec: Dec in deg
    max_ellip: Maximum elliptical something?
    box_size: Half-width of box for reference catalog query (arcsec)
    max_rad: Maximum distance to look for star pairs.
    tolerance: Amount of slack allowed in match determination
    catalog: Catalog to use (ub2, tmc, sdss, or file)
    overwrite: Overwrite output files
    outfile: Output file
    output_dir: Output directory
    saturation: Saturation level; do not use stars exceeding.
    no_rot: Some kind of bool involving rotation
    write_crosscheck_files: Bool for whether to write region and other crosscheck files

    Returns
    -------

    """

    if isinstance(files, str):
        files = [files]

    if len(files) == 0:
        err = "No files selected!"
        logger.error(err)
        raise ValueError(err)

    if np.logical_and(overwrite, outfile is not None):
        err = (
            f"An output file was specified ({outfile}), "
            f"but the script was configured to overwrite the original file."
        )
        logger.error(err)
        raise ValueError(err)

    n_image = len(files)
    failures = []
    questionable = []
    multi_info = []

    for img_path in files:
        if len(files) > 1:
            logger.debug(f"Processing {img_path}")

        fit_info = run_autoastrometry_single(
            img_path=img_path,
            seeing=seeing,
            pixel_scale=pixel_scale,
            pa=pa,
            uncpa=uncpa,
            inv=inv,
            user_ra=user_ra,
            user_dec=user_dec,
            max_ellip=max_ellip,
            box_size=box_size,
            max_rad=max_rad,
            tolerance=tolerance,
            catalog=catalog,
            outfile=outfile,
            overwrite=overwrite,
            output_dir=output_dir,
            saturation=saturation,
            no_rot=no_rot,
            write_crosscheck_files=write_crosscheck_files,
        )

        # WTF?

        if isinstance(fit_info, int):
            fit_info = (0, 0, 0, 0, 0, 0)

        multi_info.append(fit_info)

        if fit_info[0] == 0:  # number of matches
            failures.append(img_path)
        if fit_info[5] > 2:  # stdev of offset
            questionable.append(img_path)

    if n_image > 1:
        if len(failures) == 0 and len(questionable) == 0:
            logger.debug("Successfully processed all images!")
        else:
            logger.warning("Finished processing all images, not all were successful.")

        if len(questionable) > 0:
            logger.warning(
                "The following images solved but have questionable astrometry: \n"
            )
            for warn_img in questionable:
                logger.warning(warn_img)
        if len(failures) > 0:
            logger.error("The following images failed to solve: \n")
            for fail_img in failures:
                logger.error(fail_img)

        logger.debug(f"{'Filename':25s} ")
        logger.debug(
            f"{'#match':6s} {'dPA ':8s} ({'stdev':6s})  "
            f"{'dRA':7s} {'dDec':7s} ({'stdev':6s})"
        )
        for i in range(len(files)):
            info = multi_info[i]
            logger.debug(f"{files:25s} ")
            if info[0] > 0:
                logger.debug(info)
            else:
                logger.debug("failed to solve")

    try:
        os.remove("temp.param")
    except FileNotFoundError:
        pass


######################################################################
# Running as executable
if __name__ == "__main__":
    run_autoastrometry_batch(*sys.argv)

######################################################################


# some possible future improvements:
# verification to assess low-confidence solutions
# full automatic retry mode (parity detection, etc.)
# dealing with unknown pixel scale
# run wcstools for distortion parameters
# merge catalog check with catalog search to save a query
# improve the CR rejection further... maybe think about recognizing elliptical "seeing"?
