#!/bin/env python

#  autoastrometry.py - a fast astrometric solver
#
#    author: Daniel Perley (dperley@astro.caltech.edu)
#    last significant modifications 2012-04-23
#  
#  Installation:
#     Save this file anywhere on disk, and call it from the command 
#       line: "python autoastrometry.py"
#     Required python packages:  numpy, pyfits, and optionally ephem.
#     You must also have sextractor installed: if the path is
#       nonstandard, edit the global variable below to specify. 
#     For help, type "python autoastrometry.py -help"

# 4/23: program can actually be overwhelmed by too many good matches (too high maxrad).
# need to fix this.

# Modified by Kishalay De (kde@astro.caltech.edu) for removing dependency on deprecated pyfits
# and making it compatible with astropy headers and python 3.6 (June 11, 2018)

# Modified/dissected by Robert Stein (rdstein@caltech.edu) in 2021/2022 for incorporation into winterdrp
# Also converted to python 3.10

import os
import sys
import urllib.error
import urllib.parse
import urllib.request
import math
import numpy as np
from astropy.io import fits
from winterdrp.paths import base_output_dir, ProcessingError
from winterdrp.processors.astromatic.sextractor.sourceextractor import run_sextractor_single, default_saturation
import logging
import ephem
from winterdrp.processors.astromatic.sextractor.settings import write_param_file, write_config_file, \
    default_config_path, default_conv_path, default_param_path, default_starnnw_path
from pathlib import Path

logger = logging.getLogger(__name__)

default_tolerance = 0.01  # these defaults should generally not be altered.
default_pa_tolerance = 1.4
default_min_fwhm = 1.5
default_max_fwhm = 40

fast_match = True
show_matches = False


class AstrometryError(ProcessingError):
    pass


class BaseSource:

    def __init__(
            self,
            ra_deg: float,
            dec_deg: float,
            in_mag: float
    ):
        self.ra_deg = float(ra_deg)
        self.dec_deg = dec_deg
        self.ra_rad = ra_deg * math.pi / 180
        self.dec_rad = dec_deg * math.pi / 180
        self.mag = in_mag

    def rotate(
            self,
            dpa_deg: float,
            ra0: float,
            dec0: float
    ):
        dpa_rad = dpa_deg * math.pi/180
        sin_dpa = math.sin(dpa_rad)
        cos_dpa = math.cos(dpa_rad)
        ra_scale = math.cos(dec0*math.pi/180)

        # this is only valid for small fields away from the pole.
        x = (self.ra_deg - ra0) * ra_scale
        y = (self.dec_deg - dec0)

        x_rot = cos_dpa * x - sin_dpa * y
        y_rot = sin_dpa * x + cos_dpa * y

        self.ra_deg = (x_rot / ra_scale) + ra0
        self.dec_deg = y_rot + dec0
        self.ra_rad = self.ra_deg * math.pi / 180
        self.dec_rad = self.dec_deg * math.pi / 180


class SextractorSource(BaseSource):

    def __init__(
            self,
            line: str
    ):
        inline_arg = [x.strip() for x in line.split(" ") if x not in [""]]

        if len(inline_arg) < 8:
            err = f"Expected 8 values in table, found {len(inline_arg)} ({inline_arg})"
            logger.error(err)
            raise ValueError(err)

        self.x = float(inline_arg[0])
        self.y = float(inline_arg[1])

        super().__init__(*[float(x) for x in inline_arg[2:5]])

        self.mag_err = float(inline_arg[5])
        self.ellip = float(inline_arg[6])
        self.fwhm = float(inline_arg[7])

        if len(inline_arg) >= 9:
            self.flag = int(inline_arg[8])
        else:
            self.flag = None


# Pixel distance
def pixel_distance(
        obj1: SextractorSource,
        obj2: SextractorSource
) -> float:
    return ((obj1.x - obj2.x)**2 + (obj1.y - obj2.y)**2)**0.5


# Great circle distance between two points.
def distance(
        obj1: BaseSource,
        obj2: BaseSource
) -> float:

    ddec = obj2.dec_rad - obj1.dec_rad
    dra = obj2.ra_rad - obj1.ra_rad
    dist_rad = 2 * math.asin(math.sqrt(
        (math.sin(ddec/2.))**2 +
        math.cos(obj1.dec_rad) * math.cos(obj2.dec_rad) * (math.sin(dra/2.))**2
    ))

    dist_deg = dist_rad * 180. / math.pi
    dist_arc_sec = dist_deg * 3600.
    return dist_arc_sec


# Non-great-circle distance is much faster
def quickdistance(
        obj1: BaseSource,
        obj2: BaseSource,
        cosdec: float
) -> float:
    ddec = obj2.dec_deg - obj1.dec_deg
    dra = obj2.ra_deg - obj1.ra_deg
    if dra > 180:
        dra = 360 - dra
    return 3600 * math.sqrt(ddec**2 + (cosdec*dra)**2)


# Calculate the (spherical) position angle between two objects.
def position_angle(
        obj1: BaseSource,
        obj2: BaseSource
) -> float:
    dra = obj2.ra_rad - obj1.ra_rad
    pa_rad = np.arctan2(
        math.cos(obj1.dec_rad) * math.tan(obj2.dec_rad) - math.sin(obj1.dec_rad) * math.cos(dra),
        math.sin(dra)
    )
    pa_deg = pa_rad * 180./math.pi
    pa_deg = 90. - pa_deg  # defined as degrees east of north
    while pa_deg > 200:
        pa_deg -= 360.   # make single-valued
    while pa_deg < -160:
        pa_deg += 360.  # note there is a crossing point at PA=200, images at this exact PA
    return pa_deg                        # will have the number of matches cut by half at each comparison level


# Compare objects using magnitude.
def compare_mag(
        source: SextractorSource
) -> float:
    """"useful for sorting; Altered by KD for compatibility with python 3"""
    return source.mag


def median(
        float_list: list[float]
) -> float:
    a = np.array(float_list)
    return float(np.median(a))


def stdev(
        float_list: list[float]
) -> float:
    a = np.array(float_list)
    return float(np.std(a))


def mode(
        float_list: list[float]
) -> float:

    if len(float_list) == 0:
        err = "Float list is empty, cannot calculate mode."
        logger.error(err)
        raise ValueError(err)

    s = np.array(sorted(float_list))
    d = s[1:] - s[:-1]
    nd = len(d)
    if nd >= 32:
        g = nd/16
    elif nd >= 6:
        g = 2
    else:
        g = 1

    min_mean = d.sum()
    i_mean = nd / 2

    for i in range(nd):
        r = [int(max(i-g, 0)), int(min(i+g, nd))]
        m = d[r[0]:r[1]].mean()
        if m < min_mean:
            min_mean = m
            i_mean = i

    list_mode = s[int(i_mean)]  # + s[i_mean+1])/2

    return list_mode


def ra_str_2_deg(
        ra_str: str
) -> float:
    ra_str = str(ra_str).strip()
    ra = ra_str.split(':')
    if len(ra) == 1:
        return float(ra_str)
    return 15*(float(ra[0])+float(ra[1])/60.0 + float(ra[2])/3600.0)


def dec_str_2_deg(
        dec_str: str
) -> float:
    dec_str = str(dec_str).strip()
    dec = dec_str.split(':')
    if len(dec) == 1:
        return float(dec_str)
    sign = 1
    if dec_str[0] == '-':
        sign = -1
    return sign*(abs(float(dec[0])) + float(dec[1])/60.0 + float(dec[2])/3600.0)


def unique(
        input_list: list
) -> list:
    lis = input_list[:]  # make a copy
    lis.sort()
    llen = len(lis)
    i = 0
    while i < llen-1:
        if lis[i+1] == lis[i]:
            del lis[i+1]
            llen = llen - 1
        else:
            i = i + 1
    return lis


def get_img_src_list(
        img_path: str,
        base_output_path: str,
        nx_pix: int,
        ny_pix: int,
        border: int = 3,
        corner: int = 12,
        min_fwhm: float = default_min_fwhm,
        max_fwhm: float = default_max_fwhm,
        max_ellip: float = 0.5,
        saturation: float = default_saturation,
        config_path: str = default_config_path,
        output_catalog: str = None,
        write_crosscheck_files: bool = False
) -> list[SextractorSource]:

    if output_catalog is None:
        output_catalog = Path(base_output_path).with_suffix(".cat")

    try:
        os.remove(output_catalog)
    except FileNotFoundError:
        pass

    # cmd = f"{sextractor_cmd} {sexfilename} -c {config_path} -SATUR_LEVEL {saturation} -CATALOG_NAME {output_catalog}"
    # execute_sextractor(cmd, output_dir=output_dir)

    run_sextractor_single(
        img=img_path,
        output_dir=os.path.dirname(output_catalog),
        config=config_path,
        saturation=saturation,
        catalog_name=output_catalog,
        parameters_name=default_param_path,
        filter_name=default_conv_path,
        starnnw_name=default_starnnw_path
    )

    # Read in the sextractor catalog
    with open(output_catalog, 'rb') as cat:
        catlines = [x.replace(b"\x00", b"").decode() for x in cat.readlines()][1:]

    # Delete the sextractor catalog again, if not requested
    if not write_crosscheck_files:
        os.remove(output_catalog)

    if len(catlines) == 0:
        logger.error('Sextractor catalog is empty: try a different catalog?')
        raise ValueError

    min_x = border
    min_y = border
    max_x = nx_pix - border    # This should be generalized
    max_y = ny_pix - border

    n_src_init = 0
    n_src_pass = 0
    src_list = []

    rejects = []

    for line in catlines:

        if line[0] == "#":
            continue

        src = SextractorSource(line)  # process the line into an object
        n_src_init += 1

        # Initial filtering
        if src.ellip > max_ellip:
            rejects.append("ellip")
            continue
        if src.fwhm < min_fwhm:
            rejects.append("min fwhm")
            continue
        if src.fwhm > max_fwhm:
            rejects.append("max fwhm")
            continue
        if src.x < min_x:
            rejects.append("min val")
            continue
        if src.y < min_y:
            rejects.append("min y")
            continue
        if src.x > max_x:
            rejects.append("max val")
            continue
        if src.y > max_y:
            rejects.append("max y")
            continue
        if src.x + src.y < corner:
            rejects.append("corner")
            continue
        if src.x + (ny_pix - src.y) < corner:
            rejects.append("corner")
            continue
        if (nx_pix - src.x) < corner:
            rejects.append("corner")
            continue
        if (nx_pix - src.x) + (ny_pix - src.y) < corner:
            rejects.append("corner")
            continue
        if saturation is not None:
            if src.flag > 0:
                rejects.append("saturation")
                continue  # this will likely overdo it for very deep fields.

        src_list.append(src)
        n_src_pass += 1

    # Remove detections along bad columns

    thresh_prob = 0.0001
    ct_bad_col = 0
    for i in range(5):

        for variable in ["x", "y"]:

            txp = 1.0
            val_thresh = 1
            while txp > thresh_prob:
                txp *= min((len(src_list) * 1.0 / nx_pix), 0.8)  # some strange way of estimating the threshold.
                val_thresh += 1                          # what I really want is a general analytic expression for

            remove_list = []                           # the 99.99% prob. threshold for value of n for >=n out
            val_list = [getattr(src, variable) for src in src_list]
            mode_val = mode(val_list)                       # of N total sources to land in the same bin (of NX total bins)
            for j, val in enumerate(val_list):
                if (val > mode_val-1) and (val < mode_val+1):
                    remove_list.append(j)

            remove_list.reverse()
            if len(remove_list) > val_thresh:
                for k in remove_list:
                    del src_list[k]
                    ct_bad_col += 1

    if ct_bad_col > 0:
        rejects += ["bad columns" for _ in range(ct_bad_col)]

    # Remove galaxies and cosmic rays

    fwhm_list = [src.fwhm for src in src_list]

    if len(fwhm_list) > 5:
        fwhm_list.sort()
        fwhm_20 = fwhm_list[int(len(fwhm_list)/5)]
        fwhm_mode = mode(fwhm_list)
    else:
        fwhm_mode = min_fwhm
        fwhm_20 = min_fwhm

    # formerly a max, but occasionally a preponderance of long CR's could cause fwhm_mode to be bigger than the stars
    refined_min_fwhm = median([0.75 * fwhm_mode, 0.9 * fwhm_20, min_fwhm])
    # if CR's are bigger and more common than stars, this is dangerous...
    logger.debug(f'Refined min FWHM: {refined_min_fwhm} pix')

    # Might also be good to screen for false detections created by bad columns/rows

    n_good = 0
    good_src_list = []
    for src in src_list:
        if src.fwhm > refined_min_fwhm:
            good_src_list.append(src)
            n_good += 1
        else:
            rejects.append("refined min fwhm")

    # Sort by magnitude
    good_src_list.sort(key=compare_mag)

    logger.debug(f'{n_good} objects detected in image {img_path} (a further {n_src_init - n_good} discarded)')

    reject_stats = [(x, rejects.count(x)) for x in list(set(rejects))]
    logger.debug(f"Reject reasons: {reject_stats}")

    return good_src_list


def get_catalog(
        catalog: str,
        ra: float,
        dec: float,
        box_size_arcsec: float,
        min_mag: float = 8.0,
        max_mag: float = None,
        max_pm: float = 60.
) -> list[BaseSource]:

    # Get catalog from USNO

    if max_mag is None:

        if catalog == 'ub2':
            max_mag = 21.0  # 19.5
        elif catalog == 'sdss':
            max_mag = 22.0
        elif catalog == 'tmc':
            max_mag = 20.0
        else:
            err = "Catalog not recognised. Please select 'ub2', 'sdss' or tmc'."
            logger.error(err)
            raise ValueError(err)

    ra_col = 1
    dec_col = 2

    if catalog == 'tmc':
        mag_col = 3
    else:
        mag_col = 6

    pm_ra_col = 10
    pm_dec_col = 11
    query_url = f"http://tdc-www.harvard.edu/cgi-bin/scat?catalog={catalog}" \
                f"&ra={ra}&dec={dec}&system=J2000&rad={-box_size_arcsec}" \
                f"&sort=mag&epoch=2000.00000&nstar=6400"

    with urllib.request.urlopen(query_url) as cat:
        cat_lines = cat.readlines()

    if len(cat_lines) > 6400-20:
        logger.warning('Reached maximum catalog query size. Gaps may be '
                       'present in the catalog, leading to a poor solution '
                       'or no solution. Decrease the search radius.')

    cat_list = list()

    for line in cat_lines:
        inline = line.strip()

        if len(inline) <= 2:
            continue

        if inline[0:2] == '#:':
            inline_arg = inline[2:].split(',')
            ra_col = int(inline_arg[0])-1
            dec_col = int(inline_arg[1])-1
            if len(inline_arg) > 2:
                mag_col = int(inline_arg[2])-1
            continue

        if (int(inline[0]) < ord('0') or int(inline[0]) > ord('9')) and str(inline[0]) != '.':
            continue  # this may be too overzealous about
        if (int(inline[1]) < ord('0') or int(inline[1]) > ord('9')) and str(inline[1]) != '.':
            continue  # removing comments...

        inline_arg_byte = inline.split()
        inline_arg = [str(bytes(a), 'utf-8') for a in inline_arg_byte]
        n_arg = len(inline_arg)

        if inline_arg[ra_col].find(':') == -1:
            ra = float(inline_arg[ra_col])
        else:
            ra = ra_str_2_deg(inline_arg[ra_col])

        if inline_arg[dec_col].find(':') == -1:
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


def distance_match(
        img_src_list: list[SextractorSource],
        ref_src_list: list[BaseSource],
        base_output_path: str,
        max_rad: float = 180.,
        min_rad: float = 10.,
        tolerance: float = 0.010,
        req_match: int = 3,
        pa_tolerance: float = 1.2,
        unc_pa: float = None,
        write_crosscheck_files: bool = False
) -> tuple[list[int], list[int], list[float]]:

    if tolerance <= 0:
        err = f'Tolerance cannot be negative! Using {abs(tolerance)} not {tolerance}.'
        logger.error(err)
        raise ValueError(err)
    elif pa_tolerance <= 0:
        err = f'PA tolerance cannot be negative! Using {abs(pa_tolerance)} not {pa_tolerance}.'
        logger.error(err)
        raise ValueError(err)

    if req_match < 2:
        logger.warning('Warning: reqmatch >=3 suggested')

    if unc_pa is None:
        unc_pa = 720.

    dec_list = list()
    for src in img_src_list:
        dec_list.append(src.dec_rad)
    median_dec_rad = median(dec_list)       # faster distance computation
    ra_scale = math.cos(median_dec_rad)          # will mess up meridian crossings, however

    # Calculate all the distances

    # In image catalog:
    img_src_dists = list()
    img_src_match_ids = list()
    for i, src in enumerate(img_src_list):
        d = []
        dj = []
        for j, src2 in enumerate(img_src_list):
            if i == j:
                continue

            if abs(src.dec_deg - src2.dec_deg) > max_rad:
                continue
            if ra_scale*abs(src.ra_deg - src2.ra_deg) > max_rad:
                continue

            dist = quickdistance(src, src2, ra_scale)

            if min_rad < dist < max_rad:
                d.append(dist)
                dj.append(j)

        img_src_dists.append(d)
        img_src_match_ids.append(dj)

    # In reference catalog:
    ref_src_dists = []
    ref_src_match_ids = []
    for i, src in enumerate(ref_src_list):
        d = []
        dj = []
        for j, src2 in enumerate(ref_src_list):

            if i == j:
                continue

            if abs(src.dec_deg - src2.dec_deg) > max_rad:
                continue
            if ra_scale*abs(src.ra_deg - src2.ra_deg) > max_rad:
                continue

            dist = quickdistance(src, src2, ra_scale)

            if min_rad < dist < max_rad:
                d.append(dist)
                dj.append(j)

        ref_src_dists.append(d)
        ref_src_match_ids.append(dj)

    # Now look for matches in the reference catalog to distances in the image catalog.

    n_great_matches = 0

    img_match = list()
    ref_match = list()
    mpa = list()
    match_ns = list()

    primary_match_img = []
    primary_match_ref = []

    for img_i, img_dist_array in enumerate(img_src_dists):

        if len(img_dist_array) < 2:
            continue

        for ref_i, ref_dist_array in enumerate(ref_src_dists):

            if len(ref_dist_array) < 2:
                continue

            match = 0
            img_match_in = []
            ref_match_in = []

            for img_j, img_src_dist in enumerate(img_dist_array):

                new_match = 1

                for ref_j, ref_src_dist in enumerate(ref_dist_array):

                    if abs((img_src_dist / ref_src_dist) - 1.0) < tolerance:
                        match += new_match

                        new_match = 0  # further matches before the next img_j loop indicate degeneracies

                        img_match_in.append(img_src_match_ids[img_i][img_j])
                        ref_match_in.append(ref_src_match_ids[ref_i][ref_j])

            if match >= req_match:

                dpa = []
                # Here, dpa[n] is the mean rotation of the PA from the primary star of this match
                #  to the stars in its match RELATIVE TO those same angles for those same stars
                #  in the catalog.  Therefore it is a robust measurement of the rotation.

                for i in range(len(img_match_in)):
                    ddpa = position_angle(img_src_list[img_i], img_src_list[img_match_in[i]]) \
                           - position_angle(ref_src_list[ref_i], ref_src_list[ref_match_in[i]])
                    while ddpa > 200.:
                        ddpa -= 360.
                    while ddpa < -160.:
                        ddpa += 360.
                    dpa.append(ddpa)

                # If user was confident the initial PA was right, remove bad PA'src right away
                for i in range(len(img_match_in)-1, -1, -1):
                    if abs(dpa[i]) > unc_pa:
                        del img_match_in[i]
                        del ref_match_in[i]
                        del dpa[i]

                if len(img_match_in) < 2:
                    continue

                mode_dpa = mode(dpa)

                # Remove deviant matches by PA
                for i in range(len(img_match_in)-1, -1, -1):
                    if abs(dpa[i] - mode_dpa) > pa_tolerance:
                        del img_match_in[i]
                        del ref_match_in[i]
                        del dpa[i]

                if len(img_match_in) < 2:
                    continue

                n_degeneracies = (
                        len(img_match_in)
                        - len(unique(img_match_in))
                        + len(ref_match_in)
                        - len(unique(ref_match_in))
                )
                # this isn't quite accurate (overestimates if degeneracies are mixed up)

                mpa.append(mode_dpa)
                primary_match_img.append(img_i)
                primary_match_ref.append(ref_i)
                img_match.append(img_match_in)
                ref_match.append(ref_match_in)
                match_ns.append(len(img_match_in) - n_degeneracies)

                if len(img_match_in) - n_degeneracies > 6:
                    n_great_matches += 1

        if n_great_matches > 16 and fast_match is True:
            break  # save processing time

    n_matches = len(img_match)
    if n_matches == 0:
        logger.error('Found no potential matches of any sort (including pairs).')
        logger.error('The algorithm is probably not finding enough real stars to solve the field.  Check seeing.')
        return [], [], []

    # Kill the bad matches
    rejects = 0

    # Get rid of matches that don't pass the reqmatch cut
    # if n_matches > 10 and max(match_ns) >= reqmatch:
    for i in range(len(primary_match_img)-1, -1, -1):
        if match_ns[i] < req_match:
            del mpa[i]
            del primary_match_img[i]
            del primary_match_ref[i]
            del img_match[i]
            del ref_match[i]
            del match_ns[i]

    if len(img_match) < 1:
        logger.error(f'Found no matching clusters of reqmatch = {req_match}')
        return [], [], []

    # If we still have lots of matches, get rid of those with the minimum number of submatches
    # (that is, increase reqmatch by 1)
    min_match = min(match_ns)
    count_not_min = 0

    for n in match_ns:
        if n > min_match:
            count_not_min += 1

    if len(match_ns) > 16 and count_not_min > 3:
        logger.debug(f'Too many matches: increasing reqmatch to {req_match + 1}')
        for i in range(len(primary_match_img)-1, -1, -1):
            if match_ns[i] == min_match:
                del mpa[i]
                del primary_match_img[i]
                del primary_match_ref[i]
                del img_match[i]
                del ref_match[i]
                del match_ns[i]

    n_matches = len(img_match)  # recalculate with the new reqmatch and with prunes supposedly removed
    logger.debug(f'Found {n_matches} candidate matches.')

    # Use only matches with a consistent PA

    offset_pa = mode(mpa)

    if len(img_match) > 2:
        # Coarse iteration for anything away from the mode
        for i in range(len(primary_match_img)-1, -1, -1):
            if abs(mpa[i] - offset_pa) > pa_tolerance:
                del mpa[i]
                del primary_match_img[i]
                del primary_match_ref[i]
                del img_match[i]
                del ref_match[i]
                del match_ns[i]
                rejects += 1

        st_dev_pa = stdev(mpa)
        refined_tolerance = (2.2 * st_dev_pa)

        # Fine iteration to flag outliers now that we know most are reliable
        for i in range(len(primary_match_img)-1, -1, -1):
            if abs(mpa[i] - offset_pa) > refined_tolerance:
                del mpa[i]
                del primary_match_img[i]
                del primary_match_ref[i]
                del img_match[i]
                del ref_match[i]
                del match_ns[i]
                rejects += 1  # these aren't necessarily bad, just making more manageable.

    # New verification step: calculate distances and PAs between central stars of matches
    n_dist_flags = [0] * len(primary_match_img)
    for v in range(2):  # two iterations
        # find bad pairs
        if len(primary_match_img) == 0:
            break

        for i, img_i in enumerate(primary_match_img):
            for j, img_j in enumerate(primary_match_img):

                if i == j:
                    continue

                ref_i = primary_match_ref[i]
                ref_j = primary_match_ref[j]

                img_dist_ij = distance(img_src_list[img_i], img_src_list[img_j])
                ref_dist_ij = distance(ref_src_list[ref_i], ref_src_list[ref_j])

                try:
                    if abs((img_dist_ij / ref_dist_ij) - 1.0) > tolerance:
                        n_dist_flags[i] += 1
                except ZeroDivisionError:  # (occasionally will get divide by zero)
                    pass

        # delete bad clusters
        n_test_matches = len(primary_match_img)
        for i in range(n_test_matches-1, -1, -1):
            if n_dist_flags[i] == n_test_matches-1:   # if every comparison is bad, this is a bad match
                del mpa[i]
                del primary_match_img[i]
                del primary_match_ref[i]
                del img_match[i]
                del ref_match[i]
                del match_ns[i]
                rejects += 1

    logger.debug(f'Rejected {rejects} bad matches.')
    n_matches = len(primary_match_img)
    logger.debug(f'Found {n_matches} good matches.')

    if n_matches == 0:
        return [], [], []

    # check the pixel scale while we're at it
    pix_scale_list = []
    if len(primary_match_img) >= 2:
        for i in range(len(primary_match_img)-1):
            for j in range(i+1, len(primary_match_img)):
                img_i = primary_match_img[i]
                ref_i = primary_match_ref[i]
                img_j = primary_match_img[j]
                ref_j = primary_match_ref[j]

                pix_scale_list.append(
                    distance(ref_src_list[ref_i], ref_src_list[ref_j])
                    / pixel_distance(img_src_list[img_i], img_src_list[img_j])
                )

        pix_scale = median(pix_scale_list)
        pix_scale_std = stdev(pix_scale_list)

        if len(primary_match_img) >= 3:
            logger.debug(f'Refined pixel scale measurement: {pix_scale:.4f}"/pix '
                         f'(+/- {pix_scale_std:.4f})')
        else:
            logger.debug(f'Refined pixel scale measurement: {pix_scale:.4f}"/pix')

    for i, img_i in enumerate(primary_match_img):

        ref_i = primary_match_ref[i]

        if show_matches:
            logger.debug(f'{img_i} matches {ref_i} (dPA ={mpa[i]:.3f})')
            if len(img_match[i]) < 16:
                logger.debug(f'  {img_i} --> {img_match[i]}')
                logger.debug(f'  {ref_i} --> {ref_match[i]}')
            else:
                logger.debug(f'  {img_i} --> {img_match[i][0:10]} {len(img_match[i])-10} more')
                logger.debug(f'  {ref_i} --> {ref_match[i][0:10]} {len(ref_match[i])-10} more')
            if i+1 >= 10 and len(primary_match_img)-10 > 0:
                logger.debug(f''
                             f'{(len(primary_match_img)-10)} additional matches not shown.')
                break
        else:
            logger.debug(
                f'{img_i} matches {ref_i} (dPA ={mpa[i]:.3f}):'
                f' {str(len(img_match[i])).strip()} rays'
            )

    if write_crosscheck_files:

        match_lines_im = os.path.splitext(base_output_path)[0] + '.matchlines.im.reg'

        logger.info(f"Writing match lines to {match_lines_im}")

        with open(match_lines_im, 'w') as out:
            color = 'red'
            out.write(
                f'# Region file format: DS9 version 4.0\nglobal color={color} '
                f'font="helvetica 10 normal" select=1 highlite=1 '
                f'edit=1 move=1 delete=1 include=1 fixed=0 source\n'
            )
            out.write('image\n')
            for i, img_i in enumerate(primary_match_img):
                for j, img_j in enumerate(img_match[i]):
                    out.write(
                        f"line({img_src_list[img_i].x:.3f},"
                        f"{img_src_list[img_i].y:.3f},"
                        f"{img_src_list[img_j].x:.3f},"
                        f"{img_src_list[img_j].y:.3f}) # line=0 0\n"
                    )

        match_lines_wcs = os.path.splitext(base_output_path)[0]+'.matchlines.wcs.reg'

        logger.info(f"Writing match lines to {match_lines_wcs}")

        with open(match_lines_wcs, 'w') as out:
            color = 'green'
            out.write(
                f'# Region file format: DS9 version 4.0\nglobal color={color} '
                f'font="helvetica 10 normal" select=1 highlite=1 '
                f'edit=1 move=1 delete=1 include=1 fixed=0 source\n'
            )
            out.write('fk5\n')
            for i, ref_i in enumerate(primary_match_ref):
                for j, ref_j in enumerate(ref_match[i]):
                    out.write(
                        f"line({ref_src_list[ref_i].ra_deg:.5f},"
                        f"{ref_src_list[ref_i].dec_deg:.5f},"
                        f"{ref_src_list[ref_j].ra_deg:.5f},"
                        f"{ref_src_list[ref_j].dec_deg:.5f}) # line=0 0\n"
                    )

    # future project: if not enough, go to the secondary offsets

    return primary_match_img, primary_match_ref, mpa

############################################


def write_text_file(
        file_path: str,
        src_list: list[BaseSource]
):

    logger.info(f"Saving text file to {file_path}")

    with open(file_path, 'w') as out:
        for ob in src_list:
            out.write(f"{ob.ra_deg:11.7f} {ob.dec_deg:11.7f} {ob.mag:5.2f}\n")


def write_region_file(
        file_path: str,
        src_list: list[BaseSource],
        color: str = "green",
        system: str = None
):
    if system is None:
        system = 'wcs'

    if system not in ["wcs", "img"]:
        err = f"Did not recognise system '{system}'. Valid values are 'wcs' and 'img'."
        logger.error(err)
        raise ValueError(err)

    logger.debug(f"Saving region file to {file_path}")

    with open(file_path, 'w') as out:
        out.write(
            f'# Region file format: DS9 version 4.0\n'
            f'global color={color} font="helvetica 10 normal" select=1 '
            f'highlite=1 edit=1 move=1 delete=1 include=1 fixed=0 source\n'
        )

        if system == 'wcs':
            out.write('fk5\n')
            for i, src in enumerate(src_list):
                out.write(f"point({src.ra_deg:.7f},{src.dec_deg:.7f}) # point=boxcircle text={{{i+1}}}\n")
        elif system == 'img':
            out.write('image\n')
            for i, src in enumerate(src_list):
                out.write(f"point({src.ra_deg:.3f},{src.dec_deg:.3f}) # point=boxcircle text={{{i+1}}}\n")


def parse_header(
        file_path: str,
        temp_path: str,
        pixel_scale: float = None,
        pa: float = None,
        inv: bool = False,
        user_ra_deg: float = None,
        user_dec_deg: float = None,
):
    sci_ext = 0

    # Get some basic info from the header
    with fits.open(file_path) as hdu:
        hdu.verify('silentfix')

        header = hdu[sci_ext].header

        if np.logical_and(pixel_scale is not None, pa is None):
            pa = 0

        # Check for old-style WCS header
        if pixel_scale is None:

            old_wcs_type = False

            for hkey in header.keys():
                if hkey in ['CDELT1', 'CDELT2']:
                    old_wcs_type = True

            if old_wcs_type:
                key = 'CDELT1'
                cdelt1 = header[key]
                key = 'CDELT2'
                cdelt2 = header[key]

                try:
                    c_rot = 0
                    key = 'CROTA1'
                    c_rot = header[key]
                    key = 'CROTA2'
                    c_rot = header[key]
                except KeyError:
                    pass

                if math.sqrt(cdelt1 ** 2 + cdelt2 ** 2) < 0.1:  # some images use CDELT to indicate nonstandard things
                    header['CD1_1'] = cdelt1 * math.cos(c_rot * math.pi / 180.)
                    header['CD1_2'] = -cdelt2 * math.sin(c_rot * math.pi / 180.)
                    header['CD2_1'] = cdelt1 * math.sin(c_rot * math.pi / 180.)
                    header['CD2_2'] = cdelt2 * math.cos(c_rot * math.pi / 180.)

        if np.logical_and(pixel_scale is not None, pa is not None):
            # Create WCS header information if pixel scale is specified
            pa_rad = pa * math.pi / 180.
            px_scale_deg = pixel_scale / 3600.

            if inv > 0:
                parity = -1
            else:
                parity = 1

            if user_ra_deg is not None:
                ra = user_ra_deg
            else:
                ra = ra_str_2_deg(header['CRVAL1'])

            if user_dec_deg is not None:
                dec = user_dec_deg
            else:
                dec = dec_str_2_deg(header['CRVAL2'])

            try:
                epoch = float(header.get('EPOCH', 2000))
            except KeyError:
                logger.warning("No EPOCH found in header. Assuming 2000")
                epoch = 2000.

            try:
                equinox = float(header.get('EQUINOX', epoch))  # If RA and DEC are not J2000 then convert
            except KeyError:
                logger.warning("No EQUINOX found in header. Assuming 2000")
                equinox = 2000.  # could be 'J2000'; try to strip off first character?

            if abs(equinox - 2000) > 0.5:
                logger.debug(f'Converting equinox from {equinox} to J2000')
                j2000 = ephem.Equatorial(ephem.Equatorial(
                    str(ra / 15), str(dec), epoch=str(equinox)), epoch=ephem.J2000
                )
                [ra, dec] = [ra_str_2_deg(j2000.ra), dec_str_2_deg(j2000.dec)]

            header["CD1_1"] = px_scale_deg * math.cos(pa_rad) * parity
            header["CD1_2"] = px_scale_deg * math.sin(pa_rad)
            header["CD2_1"] = -px_scale_deg * math.sin(pa_rad) * parity
            header["CD2_2"] = px_scale_deg * math.cos(pa_rad)
            header["CRPIX1"] = header['NAXIS1'] / 2
            header["CRPIX2"] = header['NAXIS2'] / 2
            header["CRVAL1"] = ra
            header["CRVAL2"] = dec
            header["CTYPE1"] = "RA---TAN"
            header["CTYPE2"] = "DEC--TAN"
            header["EQUINOX"] = 2000.0

            hdu[sci_ext].header = header
            hdu.writeto(temp_path, output_verify='silentfix', overwrite=True)  # ,clobber=True

        # Read the header info from the file.
        try:
            # no longer drawing RA and DEC from here.
            key = 'NAXIS1'
            nxpix = header[key]
            key = 'NAXIS2'
            nypix = header[key]
        except KeyError:
            err = f'Cannot find necessary WCS header keyword {key}'
            logger.debug(err)
            raise

        try:
            key = 'CRVAL1'
            cra = float(header[key])
            key = 'CRVAL2'
            cdec = float(header[key])

            key = 'CRPIX1'
            crpix1 = float(header[key])
            key = 'CRPIX2'
            crpix2 = float(header[key])

            key = 'CD1_1'
            cd11 = float(header[key])
            key = 'CD2_2'
            cd22 = float(header[key])
            key = 'CD1_2'
            cd12 = float(header[key])  # deg / pix
            key = 'CD2_1'
            cd21 = float(header[key])

            equinox = float(header.get('EQUINOX', 2000.))
            if abs(equinox - 2000.) > 0.2:
                logger.debug('Warning: EQUINOX is not 2000.0')

        except KeyError:
            if pixel_scale == -1:
                err = f"Cannot find necessary WCS header keyword '{key}' \n " \
                      f"Must specify pixel scale (-px VAL) or provide provisional basic WCS info via CD matrix."
                logger.error(err)
                raise
                # Some images might use CROT parameters, could try to be compatible with this too...?

        # Wipe nonstandard hdu info from the header (otherwise this will confuse verification)
        header_keys = list(header.keys())
        ctype_change = 0
        iraf_keys = []
        high_keys = []
        old_keys = []
        distortion_keys = []
        for hkey in header_keys:
            if hkey == 'RADECSYS' or \
                    hkey == 'WCSDIM' or \
                    hkey.find('WAT') == 0 or \
                    hkey.find('LTV') >= 0 or \
                    hkey.find('LTM') == 0:
                del header[hkey]
                iraf_keys.append(hkey)

            if hkey.find('CO1_') == 0 or \
                    hkey.find('CO2_') == 0 or \
                    hkey.find('PV1_') == 0 or \
                    hkey.find('PV2_') == 0 or \
                    hkey.find('PC00') == 0:
                del header[hkey]
                high_keys.append(hkey)

            if hkey.find('CDELT1') == 0 or \
                    hkey.find('CDELT2') == 0 or \
                    hkey.find('CROTA1') == 0 or \
                    hkey.find('CROTA2') == 0:
                del header[hkey]
                old_keys.append(hkey)

            if hkey.find('A_') == 0 or \
                    hkey.find('B_') == 0 or \
                    hkey.find('AP_') == 0 or \
                    hkey.find('BP_') == 0:
                del header[hkey]
                distortion_keys.append(hkey)

        if header['CTYPE1'] != 'RA---TAN':
            logger.info(f"Changing CTYPE1 from '{header['CTYPE1']}' to 'RA---TAN'")
            header["CTYPE1"] = "RA---TAN"
            ctype_change = 1

        if header['CTYPE2'] != 'DEC--TAN':
            if ctype_change:
                logger.debug(f"Changing CTYPE2 from '{header['CTYPE2']}' to 'DEC--TAN'")
            header["CTYPE2"] = "DEC--TAN"
            ctype_change = 1

        wcs_key_check = [
            'CRVAL1',
            'CRVAL2',
            'CRPIX1',
            'CRPIX2',
            'CD1_1',
            'CD1_2',
            'CD2_2',
            'CD2_1',
            'EQUINOX',
            'EPOCH'
        ]

        header_format_change = False

        for w in wcs_key_check:
            if isinstance(w, str):
                try:
                    header[w] = float(header[w])
                    header_format_change = True
                except KeyError:
                    pass

        if len(iraf_keys) > 0:
            logger.warning(f'Removed nonstandard WCS keywords: {iraf_keys}')
        if len(high_keys) > 0:
            logger.warning(f'Removed higher-order WCS keywords: {high_keys}')
        if len(old_keys) > 0:
            logger.warning(f'Removed old-style WCS keywords: {old_keys}')
        if len(distortion_keys) > 0:
            logger.warning(f'Removed distortion WCS keywords: {distortion_keys}')

        if len(high_keys) + len(distortion_keys) + ctype_change + header_format_change > 0:
            # Rewrite and reload the image if the header was modified in a significant way so
            # sextractor sees the same thing that we do.
            hdu[sci_ext].header = header
            hdu.writeto(temp_path, output_verify='silentfix', overwrite=True)  # ,clobber=True

    return nxpix, nypix, cd11, cd12, cd21, cd22, crpix1, crpix2, cra, cdec


def get_ref_sources_from_catalog(
        catalog: str,
        center_ra: float,
        center_dec: float,
        box_size_arcsec: float,
):
    # If no catalog specified, check availability of SDSS
    if catalog is None:
        try:
            trycats = ['sdss', 'ub2', 'tmc']
            for trycat in trycats:
                testqueryurl = f"http://tdc-www.harvard.edu/cgi-bin/scat?catalog={trycat}&ra={center_ra}" \
                               f"&dec={center_dec}&system=J2000&rad=-90"

                with urllib.request.urlopen(testqueryurl, timeout=30) as check:
                    checklines = check.readlines()

                if len(checklines) > 15:
                    catalog = trycat
                    logger.info(f'Using catalog {catalog}')
                    break
        except urllib.error.URLError:
            err = 'No catalog is available.  Check your internet connection.'
            logger.error(err)
            raise AstrometryError(err)

    # Load in reference star catalog

    ref_src_list = get_catalog(
        catalog=catalog,
        ra=center_ra,
        dec=center_dec,
        box_size_arcsec=box_size_arcsec
    )

    n_cat = len(ref_src_list)
    cat_density = n_cat / (2 * box_size_arcsec / 60.) ** 2
    logger.debug(f'{n_cat} good catalog objects.')
    logger.debug(f'Source density of {cat_density} /arcmin^2')

    if n_cat == 0:
        logger.error('No objects found in catalog.')
        logger.error('The web query failed, all stars were excluded by the FHWM clip, or the image')
        logger.error('is too small.  Check input parameters or your internet connection.')
        raise AstrometryError

    elif 0 < n_cat < 5:
        logger.error(f'Only {n_cat} catalog objects in the search zone.'
                     f'Increase the magnitude threshold or box size.')
        raise AstrometryError

    return ref_src_list, n_cat, cat_density


def crosscheck_source_lists(
        img_src_list: list[SextractorSource],
        n_img: int,
        img_density: float,
        ref_src_list: list[BaseSource],
        n_ref: int,
        ref_density: float,
        box_size_arcsec: float,
        area_sq_min: float
) -> tuple[list[SextractorSource], int, float, list[BaseSource], int, float]:

    # If this image is actually shallower than reference catalog, trim the reference catalog down
    if n_ref > 16 and ref_density > 3 * img_density:
        logger.debug('Image is shallow.  Trimming reference catalog...')
        while ref_density > 3 * img_density:
            ref_src_list = ref_src_list[0:int(len(ref_src_list)*4/5)]
            n_ref = len(ref_src_list)
            ref_density = n_ref / (2 * box_size_arcsec / 60.) ** 2

    # If the image is way deeper than USNO, trim the image catalog down
    if n_img > 16 and img_density > 4 * ref_density:
        logger.debug('Image is deep.  Trimming image catalog...')
        while img_density > 4 * ref_density and n_img > 8:
            img_src_list = img_src_list[0:int(len(img_src_list) * 4 / 5)]
            n_img = len(img_src_list)
            img_density = n_img / area_sq_min

    # If too many objects, do some more trimming
    if n_img*n_ref > 120*120*4:
        logger.debug('Image and/or catalog still too deep.  Trimming...')
        while n_img*n_ref > 120*120*4:
            if img_density > ref_density:
                img_src_list = img_src_list[0:int(len(img_src_list) * 4 / 5)]
                n_img = len(img_src_list)
                img_density = n_img / area_sq_min
            else:
                ref_src_list = ref_src_list[0:int(len(ref_src_list)*4/5)]
                n_ref = len(ref_src_list)
                ref_density = n_ref / (2 * box_size_arcsec / 60.) ** 2

    # Remove fainter object in close pairs for both lists
    min_sep = 3

    delete_list = []
    for i, src in enumerate(img_src_list):
        for j, src2 in enumerate(img_src_list[i + 1:]):
            dist = distance(src, src2)
            if dist < min_sep:
                if src.mag > src2.mag:
                    delete_list.append(i)
                else:
                    delete_list.append(i + j + 1)

    delete_list = unique(delete_list)
    delete_list.reverse()
    for d in delete_list:
        del img_src_list[d]

    for i, src in enumerate(ref_src_list):
        for j, src2 in enumerate(ref_src_list[i+1:]):

            dist = distance(src, src2)
            if dist < min_sep:
                if src.mag > src2.mag:
                    delete_list.append(i)
                else:
                    delete_list.append(i + j + 1)

    delete_list = unique(delete_list)
    delete_list.reverse()
    for d in delete_list:
        del ref_src_list[d]

    return img_src_list, n_img, img_density, ref_src_list, n_ref, ref_density


def export_src_lists(
        img_src_list: list[SextractorSource],
        ref_src_list: list[BaseSource],
        base_output_path: str
):
    write_text_file(
        file_path=os.path.splitext(base_output_path)[0] + '.det.init.txt',
        src_list=img_src_list
    )
    write_region_file(
        file_path=os.path.splitext(base_output_path)[0] + '.det.im.reg',
        src_list=img_src_list,
        color='red',
        system='img'
    )
    write_text_file(
        file_path=os.path.splitext(base_output_path)[0] + '.cat.txt',
        src_list=ref_src_list
    )
    write_region_file(
        file_path=os.path.splitext(base_output_path)[0] + '.cat.wcs.reg',
        src_list=ref_src_list,
        color='green',
        system='wcs'
    )


############################################
def autoastrometry(
        filename: str,
        pixel_scale: float = None,
        pa: float = None,
        inv: bool = False,
        unc_pa: float = None,
        user_ra_deg: float = None,
        user_dec_deg: float = None,
        max_ellip: float = 0.5,
        box_size_arcsec: float = None,
        max_rad: float = None,
        tolerance: float = default_tolerance,
        catalog: str = None,
        overwrite: bool = True,
        outfile: str = "",
        output_dir: str = base_output_dir,
        temp_file: str = None,
        saturation: float = default_saturation,
        no_rot: bool = False,
        min_fwhm: float = default_min_fwhm,
        max_fwhm: float = default_max_fwhm,
        write_crosscheck_files: bool = False
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
        output_dir,
        ".".join(os.path.basename(filename).split(".")[:-1])
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
        user_dec_deg=user_dec_deg
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
        write_crosscheck_files=write_crosscheck_files
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
    field_width = max(x_scale * nx_pix, y_scale * ny_pix) * 3600.
    area_sq_deg = x_scale * nx_pix * y_scale * ny_pix
    area_sq_min = area_sq_deg * 3600.
    center_x = nx_pix/2
    center_y = ny_pix/2
    center_dx = center_x - crpix1
    center_dy = center_y - crpix2
    center_ra = (
            cra
            - center_dx * x_scale * math.cos(init_pa*math.pi/180.)
            + center_dy * y_scale * math.sin(init_pa*math.pi/180.)
    )
    center_dec = (
            cdec
            + parity * center_dx * x_scale * math.sin(-init_pa * math.pi/180.)
            + center_dy * y_scale * math.cos(init_pa * math.pi/180.)
    )

    if box_size_arcsec is None:
        box_size_arcsec = field_width

    # this has only been checked for a PA of zero.

    logger.debug(
        f'Initial WCS info: \n'
        f'   pixel scale:     x={x_scale*3600:.4f}"/pix,   y={y_scale*3600:.4f}"/pix \n'
        f'   position angle: PA={init_pa:.2f}'
    )

    if parity == 1:
        logger.debug('   normal parity')
    if parity == -1:
        logger.debug('   inverse parity')

    logger.debug(f'   center:        RA={center_ra:.6f}, dec={center_dec:.6f}')

    n_img = len(img_src_list)

    if n_img < 4:

        err = f'Only {n_img} good stars were found in the image.  The image is too small or shallow, the ' \
              f'detection threshold is set too high, or stars and cosmic rays are being confused.'
        logger.error(err)
        write_text_file('det.init.txt', img_src_list)
        write_region_file('det.im.reg', img_src_list, 'red', 'img')
        raise AstrometryError(err)

    img_density = len(img_src_list) / area_sq_min
    logger.debug('Source img_density of %f4 /arcmin^2' % img_density)

    # Block D

    ref_src_list, n_ref, ref_density = get_ref_sources_from_catalog(
        catalog=catalog,
        center_ra=center_ra,
        center_dec=center_dec,
        box_size_arcsec=box_size_arcsec,
    )

    # Block E

    img_src_list, n_img, img_density, ref_src_list, n_ref, ref_density = crosscheck_source_lists(
        img_src_list=img_src_list,
        n_img=n_img,
        img_density=img_density,
        ref_src_list=ref_src_list,
        n_ref=n_ref,
        ref_density=ref_density,
        box_size_arcsec=box_size_arcsec,
        area_sq_min=area_sq_min
    )

    # Block F

    if write_crosscheck_files:
        export_src_lists(
            img_src_list=img_src_list,
            ref_src_list=ref_src_list,
            base_output_path=base_output_path
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
            min_rad = 10.0   # in theory could scale this up further to reduce #comparisons
        max_rad = min(max_rad, field_width * 3. / 4)

        # note that img_density is per arcmin^2, while the radii are in arcsec, hence the conversion factor.

    circ_density = img_density * min([
        area_sq_min,
        (math.pi * (max_rad / 60.) ** 2 - math.pi * (min_rad / 60) ** 2)
    ])
    circ_ref_density = ref_density * (
            math.pi * (max_rad / 60.) ** 2 - math.pi * (min_rad / 60) ** 2
    )
    n_ref_src_per_image = ref_density * area_sq_min

    logger.debug('After trimming: ')
    logger.debug(
        f'{len(img_src_list)} detected objects '
        f'({img_density:.2f}/arcmin^2, '
        f'{circ_density:.1f}/searchzone)'
    )
    logger.debug(
        f'{len(ref_src_list)} catalog objects '
        f'({ref_density:.2f}/arcmin^2, '
        f'{circ_ref_density:.1f}/searchzone)'
    )

    pa_tolerance = default_pa_tolerance

    # RS: WTF is x**1 for?
    expect_false_trios = (
            n_img * n_ref * circ_density**2
            * circ_ref_density**2 * tolerance**2 * (pa_tolerance/360.)**1
    )

    # fraction of stars in image that are also in catalog - a guess
    overlap_first_guess = 0.3 * min(1., ref_density/img_density)
    # but how many matches >3 and >4?  some annoying binomial thing
    true_matches_per_star = (circ_density * overlap_first_guess)

    req_matches = 3
    if expect_false_trios > 30 and true_matches_per_star >= 4:
        req_matches = 4
    # should check that this will actually work for the catalog, too.
    if n_ref_src_per_image <= 6 or n_img <= 6:
        req_matches = 2
    if n_ref_src_per_image <= 3 or n_img <= 3:
        req_matches = 1
    # for an extremely small or shallow image

    logger.debug('Pair comparison search radius: %.2f"' % max_rad)
    logger.debug(f'Using req_matches = {req_matches}')

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
        write_crosscheck_files=write_crosscheck_files
    )

    n_match = len(primary_match_img)
    if n_match == 0:
        err = (' No valid matches found!\n '
               'Possible issues: \n'
               '  - The specified pixel scale (or PA or parity) is incorrect.  Double-check the input value. \n'
               '  - The field is outside the catalog search region.  Check header RA/DEC or increase search radius. \n'
               ' - The routine is flooded by bad sources.  Specify or check the input seeing. \n'
               '  - The routine is flagging many real stars.  Check the input seeing. \n'
               ' You can display a list of detected/catalog sources using det.im.reg and cat.wcs.reg. \n'
               )
        logger.error(err)
        raise AstrometryError(err)

    if n_match <= 2:
        logger.warning(f'Warning: only {n_match} match(es).  Astrometry may be unreliable.')
        logger.warning('   Check the pixel scale and parity and consider re-running.')

    # We now have the PA and a list of stars that are almost certain matches.
    median_pa = median(mpa)  # get average PA from the excellent values
    stdev_pa = stdev(mpa)

    sky_offset_pa = -parity * median_pa
    # This appears to be necessary for the printed value to agree with our normal definition.

    logger.debug('PA offset:')
    logger.debug(f'  dPA = {sky_offset_pa:.3f}  (unc. {stdev_pa:.3f})')

    if no_rot <= 0:
        # Rotate the image to the new, correct PA
        #  NOTE: when CRPIX don't match CRVAL this shifts the center and screws things up.
        #  I don't understand why they don't always match.  [[I think this was an equinox issue.
        #  should be solved now, but be alert for further problems.]]

        # Rotate....
        rot = median_pa * math.pi/180
        # ...the image itself
        header["CD1_1"] = math.cos(rot)*cd11 - math.sin(rot)*cd21
        header["CD1_2"] = math.cos(rot)*cd12 - math.sin(rot)*cd22   # a parity issue may be involved here?
        header["CD2_1"] = math.sin(rot)*cd11 + math.cos(rot)*cd21
        header["CD2_2"] = math.sin(rot)*cd12 + math.cos(rot)*cd22
        # ...the coordinates (so we don't have to resex)
        for i in range(len(img_src_list)):  # do all of them, though this is not necessary
            img_src_list[i].rotate(median_pa, cra, cdec)

    else:
        if abs(sky_offset_pa) > 1.0:
            logger.warning(' (WARNING: image appears rotated, may produce bad shift)')
        logger.debug('  Skipping rotation correction ')

    im_ra_offset = []
    im_dec_offset = []
    for i, x in enumerate(primary_match_img):
        im_ra_offset.append(img_src_list[x].ra_deg - ref_src_list[primary_match_ref[i]].ra_deg)
        im_dec_offset.append(img_src_list[x].dec_deg - ref_src_list[primary_match_ref[i]].dec_deg)

    ra_offset = -median(im_ra_offset)
    dec_offset = -median(im_dec_offset)
    ra_std = stdev(im_ra_offset) * math.cos(cdec*math.pi/180)  # all of these are in degrees
    dec_std = stdev(im_dec_offset)
    std_offset = math.sqrt(ra_std**2 + dec_std**2)

    ra_offset_arcsec = ra_offset*3600*math.cos(cdec*math.pi/180)
    dec_offset_arcsec = dec_offset*3600
    tot_offset_arcsec = (ra_offset_arcsec**2 + dec_offset**2)**0.5
    std_offset_arcsec = std_offset*3600

    logger.debug('Spatial offset:')

    msg = f'  dra = {ra_offset_arcsec:.2f}",' \
          f'  ddec = {dec_offset_arcsec:.2f}"' \
          f'  (unc. {std_offset_arcsec:.3f}")'
    logger.debug(msg)

    if std_offset*3600 > 1.0:
        logger.debug('WARNING: poor solution - some matches may be bad.  Check pixel scale?')

    header["CRVAL1"] = cra + ra_offset
    header["CRVAL2"] = cdec + dec_offset

    logger.info(f'Updated header {header["CRVAL1"],header["CRVAL2"]}')
    try:
        oldcat = header['ASTR_CAT']
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
            file_path=os.path.splitext(base_output_path)[0] + '.det.wcs.txt',
            src_list=img_src_list
        )

        # Write out a match list to allow doing a formal fit with WCStools.

        match_list_path = os.path.splitext(base_output_path)[0] + '.match.list'

        logger.info(f"Writing match list to {match_list_path}")

        with open(match_list_path, 'w') as outmatch:
            for i, si in enumerate(primary_match_img):
                ci = primary_match_ref[i]
                outmatch.write(
                    f"{img_src_list[si].x} {img_src_list[si].y} "
                    f"{ref_src_list[ci].ra_deg} {ref_src_list[ci].dec_deg}\n"
                )

    logger.info(f"Finished deriving astrometry for {filename}")

    # Could repeat with scale adjustment
    # Could then go back to full good catalog and match all sources

    if overwrite:
        outfile = filename
    elif outfile == '':
        slash_pos = filename.rfind('/')
        dir_name = filename[0:slash_pos+1]
        fil = filename[slash_pos+1:]
        outfile = f"{dir_name}a{fil}"  # alternate behavior would always output to current directory

    if outfile is not None:
        with fits.open(temp_path) as hdu:
            hdu[sci_ext].header = header

            hdu.writeto(outfile, output_verify='silentfix', overwrite=True)
            logger.info(f'Written updated file to {outfile}')

        logger.info(f"Derived center coordinates of {header['CRVAL1']}, {header['CRVAL2']}.")


    return n_match, sky_offset_pa, stdev_pa, ra_offset_arcsec, dec_offset_arcsec, std_offset_arcsec


def run_autoastrometry_single(
        img_path: str,
        seeing: float = None,
        pixel_scale: float = None,
        pa: float = None,
        uncpa: float = None,
        inv: bool = False,
        user_ra: float = None,
        user_dec: float = None,
        max_ellip: float = 0.5,
        box_size: float = None,
        max_rad: float = None,
        tolerance: float = default_tolerance,
        catalog: str = None,
        overwrite: bool = False,
        outfile: str = None,
        output_dir: str = base_output_dir,
        saturation: float = default_saturation,
        no_rot: bool = False,
        write_crosscheck_files: bool = False
):

    if seeing is None:
        min_fwhm = default_min_fwhm  # 1.5
        max_fwhm = default_max_fwhm  # 40
    else:
        min_fwhm = 0.7 * seeing
        max_fwhm = 2. * seeing

        write_param_file()
    write_config_file(saturation=saturation)
    logger.info(f'Outfile is {outfile}')
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
        write_crosscheck_files=write_crosscheck_files
    )

    return fit_info


######################################################################

def run_autoastrometry_batch(
        files: str | list[str],
        seeing: float = None,
        pixel_scale: float = None,
        pa: float = None,
        uncpa: float = None,
        inv: bool = False,
        user_ra: float = None,
        user_dec: float = None,
        max_ellip: float = 0.5,
        box_size: float = None,
        max_rad: float = None,
        tolerance: float = default_tolerance,
        catalog: str = None,
        overwrite: bool = False,
        outfile: str = None,
        output_dir: str = base_output_dir,
        saturation: float = default_saturation,
        no_rot: bool = False,
        write_crosscheck_files: bool = False
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
    adding, e.g.#:1,2,6 to the first line.  In this case, this would indicate that the RA is in the
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
    no_rot: Some kind of bool
    write_crosscheck_files: Bool for whether to write region and other crosscheck files

    Returns
    -------

    """

    if isinstance(files, str):
        files = [files]

    if len(files) == 0:
        err = 'No files selected!'
        logger.error(err)
        raise ValueError(err)

    if np.logical_and(overwrite, outfile is not None):
        err = f"An output file was specified ({outfile}), but the script was configured to overwrite the original file."
        logger.error(err)
        raise ValueError(err)

    n_image = len(files)
    failures = []
    questionable = []
    multi_info = []

    for img_path in files:

        if len(files) > 1:
            logger.debug(f'Processing {img_path}')

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
            tolerance=tolerance,
            catalog=catalog,
            outfile=outfile,
            overwrite=overwrite,
            output_dir=output_dir,
            saturation=saturation,
            no_rot=no_rot,
            write_crosscheck_files=write_crosscheck_files
        )


        # WTF?

        if isinstance(fit_info, int):
            fit_info = (0, 0, 0, 0, 0, 0)

        multi_info.append(fit_info)

        if fit_info[0] == 0:   # number of matches
            failures.append(img_path)
        if fit_info[5] > 2:    # stdev of offset
            questionable.append(img_path)

    if n_image > 1:

        if len(failures) == 0 and len(questionable) == 0:
            logger.info('Successfully processed all images!')
        else:
            logger.warning(f'Finished processing all images, not all were successful.')

        if len(questionable) > 0:
            logger.warning('The following images solved but have questionable astrometry: \n')
            for f in questionable:
                logger.warning(f)
        if len(failures) > 0:
            logger.error('The following images failed to solve: \n')
            for f in failures:
                logger.error(f)

        logger.debug("%25s " % 'Filename')
        logger.debug("%6s %8s (%6s)  %7s %7s (%6s)" % ('#match', 'dPA ', 'stdev', 'dRA', 'dDec', 'stdev'))
        for i in range(len(files)):
            info = multi_info[i]
            logger.debug(f"{files:25s} ")
            if info[0] > 0:
                logger.debug(info)
            else:
                logger.debug("failed to solve")

    try:
        os.remove('temp.param')
    except FileNotFoundError:
        pass


######################################################################
# Running as executable
if __name__ == '__main__':
    run_autoastrometry_batch(*sys.argv)

######################################################################


# some possible future improvements:
# verification to assess low-confidence solutions
# full automatic retry mode (parity detection, etc.)
# dealing with unknown pixel scale
# run wcstools for distortion parameters
# merge catalog check with catalog search to save a query
# improve the CR rejection further... maybe think about recognizing elliptical "seeing"?
