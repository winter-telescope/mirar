"""
Component module for autoastrometry, dealing with crossmatching image sources with
reference sources
"""
import logging
import os
from typing import Optional

import numpy as np

from mirar.processors.astrometry.autoastrometry.sources import (
    BaseSource,
    SextractorSource,
    distance,
    pixel_distance,
    position_angle,
    quickdistance,
)
from mirar.processors.astrometry.autoastrometry.utils import median, mode, stdev, unique

logger = logging.getLogger(__name__)

FAST_MATCH = False
SHOW_MATCH = False


def distance_match(
    img_src_list: list[SextractorSource],
    ref_src_list: list[BaseSource],
    base_output_path: str,
    max_rad: float = 180.0,
    min_rad: float = 10.0,
    tolerance: float = 0.010,
    req_match: int = 3,
    pa_tolerance: float = 1.2,
    unc_pa: Optional[float] = None,
    write_crosscheck_files: bool = False,
) -> tuple[list[int], list[int], list[float]]:
    """
    Cross-match sources in image with reference sources

    :param img_src_list: sources in image
    :param ref_src_list: reference sources
    :param base_output_path: output path to save crossmatch
    :param max_rad: max radius (degrees)
    :param min_rad: min radius (degrees)
    :param tolerance:
    :param req_match: minimum number of matches needed
    :param pa_tolerance: pa tolerance needed
    :param unc_pa: uncertainty in pa
    :param write_crosscheck_files: boolean whether to write crosscheck files
    :return: image_matches, ref_matches, mpa
    """

    if tolerance <= 0:
        err = f"Tolerance cannot be negative! Using {abs(tolerance)} not {tolerance}."
        logger.error(err)
        raise ValueError(err)

    if pa_tolerance <= 0:
        err = (
            f"PA tolerance cannot be negative! "
            f"Using {abs(pa_tolerance)} not {pa_tolerance}."
        )
        logger.error(err)
        raise ValueError(err)

    if req_match < 2:
        logger.warning("Warning: reqmatch >=3 suggested")

    if unc_pa is None:
        unc_pa = 720.0

    dec_list = []
    for src in img_src_list:
        dec_list.append(src.dec_rad)
    median_dec_rad = median(dec_list)  # faster distance computation
    ra_scale = np.cos(median_dec_rad)  # will mess up meridian crossings, however

    # Calculate all the distances

    # In image catalog:
    img_src_dists = []
    img_src_match_ids = []
    for i, src in enumerate(img_src_list):
        distances = []
        distances_j_index = []
        for j, src2 in enumerate(img_src_list):
            if i == j:
                continue

            if abs(src.dec_deg - src2.dec_deg) > max_rad:
                continue
            if ra_scale * abs(src.ra_deg - src2.ra_deg) > max_rad:
                continue

            dist = quickdistance(src, src2, ra_scale)

            if min_rad < dist < max_rad:
                distances.append(dist)
                distances_j_index.append(j)

        img_src_dists.append(distances)
        img_src_match_ids.append(distances_j_index)

    # In reference catalog:
    ref_src_dists = []
    ref_src_match_ids = []
    for i, src in enumerate(ref_src_list):
        distances = []
        distances_j_index = []
        for j, src2 in enumerate(ref_src_list):
            if i == j:
                continue

            if abs(src.dec_deg - src2.dec_deg) > max_rad:
                continue
            if ra_scale * abs(src.ra_deg - src2.ra_deg) > max_rad:
                continue

            dist = quickdistance(src, src2, ra_scale)

            if min_rad < dist < max_rad:
                distances.append(dist)
                distances_j_index.append(j)

        ref_src_dists.append(distances)
        ref_src_match_ids.append(distances_j_index)

    # Now look for matches in the reference catalog to distances in the image catalog.

    n_great_matches = 0

    img_match = []
    ref_match = []
    mpa = []
    match_ns = []

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

                        # further matches before the next
                        # img_j loop indicate degeneracies
                        new_match = 0

                        img_match_in.append(img_src_match_ids[img_i][img_j])
                        ref_match_in.append(ref_src_match_ids[ref_i][ref_j])

            if match >= req_match:
                dpa = []
                # Here, dpa[n] is the mean rotation of the PA from
                # the primary star of this match to the stars in its match
                # RELATIVE TO those same angles for those same stars
                # in the catalog.  Therefore it is a robust measurement of the rotation.

                for i, match_index in enumerate(img_match_in):
                    ddpa = position_angle(
                        img_src_list[img_i], img_src_list[match_index]
                    ) - position_angle(
                        ref_src_list[ref_i], ref_src_list[ref_match_in[i]]
                    )
                    while ddpa > 200.0:
                        ddpa -= 360.0
                    while ddpa < -160.0:
                        ddpa += 360.0
                    dpa.append(ddpa)

                # If user was confident the initial PA was right, remove bad PA'src
                # right away
                for i in range(len(img_match_in) - 1, -1, -1):
                    if abs(dpa[i]) > unc_pa:
                        del img_match_in[i]
                        del ref_match_in[i]
                        del dpa[i]

                if len(img_match_in) < 2:
                    continue

                mode_dpa = mode(dpa)

                # Remove deviant matches by PA
                for i in range(len(img_match_in) - 1, -1, -1):
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

        if n_great_matches > 16 and FAST_MATCH is True:
            break  # save processing time

    n_matches = len(img_match)
    if n_matches == 0:
        logger.error("Found no potential matches of any sort (including pairs).")
        logger.error(
            "The algorithm is probably not finding enough real stars to solve "
            "the field.  Check seeing."
        )
        return [], [], []

    # Kill the bad matches
    rejects = 0

    # Get rid of matches that don't pass the reqmatch cut
    # if n_matches > 10 and max(match_ns) >= reqmatch:
    for i in range(len(primary_match_img) - 1, -1, -1):
        if match_ns[i] < req_match:
            del mpa[i]
            del primary_match_img[i]
            del primary_match_ref[i]
            del img_match[i]
            del ref_match[i]
            del match_ns[i]

    if len(img_match) < 1:
        logger.error(f"Found no matching clusters of reqmatch = {req_match}")
        return [], [], []

    # If we still have lots of matches,
    # get rid of those with the minimum number of submatches
    # (that is, increase reqmatch by 1)
    min_match = min(match_ns)
    count_not_min = 0

    for match_n in match_ns:
        if match_n > min_match:
            count_not_min += 1

    if len(match_ns) > 16 and count_not_min > 3:
        logger.debug(f"Too many matches: increasing reqmatch to {req_match + 1}")
        for i in range(len(primary_match_img) - 1, -1, -1):
            if match_ns[i] == min_match:
                del mpa[i]
                del primary_match_img[i]
                del primary_match_ref[i]
                del img_match[i]
                del ref_match[i]
                del match_ns[i]

    n_matches = len(
        img_match
    )  # recalculate with the new reqmatch and with prunes supposedly removed
    logger.debug(f"Found {n_matches} candidate matches.")

    # Use only matches with a consistent PA

    offset_pa = mode(mpa)

    if len(img_match) > 2:
        # Coarse iteration for anything away from the mode
        for i in range(len(primary_match_img) - 1, -1, -1):
            if abs(mpa[i] - offset_pa) > pa_tolerance:
                del mpa[i]
                del primary_match_img[i]
                del primary_match_ref[i]
                del img_match[i]
                del ref_match[i]
                del match_ns[i]
                rejects += 1

        st_dev_pa = stdev(mpa)
        refined_tolerance = 2.2 * st_dev_pa

        # Fine iteration to flag outliers now that we know most are reliable
        for i in range(len(primary_match_img) - 1, -1, -1):
            if abs(mpa[i] - offset_pa) > refined_tolerance:
                del mpa[i]
                del primary_match_img[i]
                del primary_match_ref[i]
                del img_match[i]
                del ref_match[i]
                del match_ns[i]
                rejects += (
                    1  # these aren't necessarily bad, just making more manageable.
                )

    # New verification step: calculate distances and PAs between central stars
    # of matches
    n_dist_flags = [0] * len(primary_match_img)
    for _ in range(2):  # two iterations
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
        for i in range(n_test_matches - 1, -1, -1):
            if (
                n_dist_flags[i] == n_test_matches - 1
            ):  # if every comparison is bad, this is a bad match
                del mpa[i]
                del primary_match_img[i]
                del primary_match_ref[i]
                del img_match[i]
                del ref_match[i]
                del match_ns[i]
                rejects += 1

    logger.debug(f"Rejected {rejects} bad matches.")
    n_matches = len(primary_match_img)
    logger.debug(f"Found {n_matches} good matches.")

    if n_matches == 0:
        return [], [], []

    # check the pixel scale while we're at it
    pix_scale_list = []
    if len(primary_match_img) >= 2:
        for i in range(len(primary_match_img) - 1):
            for j in range(i + 1, len(primary_match_img)):
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
            logger.debug(
                f'Refined pixel scale measurement: {pix_scale:.4f}"/pix '
                f"(+/- {pix_scale_std:.4f})"
            )
        else:
            logger.debug(f'Refined pixel scale measurement: {pix_scale:.4f}"/pix')

    for i, img_i in enumerate(primary_match_img):
        ref_i = primary_match_ref[i]

        if SHOW_MATCH:
            logger.debug(f"{img_i} matches {ref_i} (dPA ={mpa[i]:.3f})")
            if len(img_match[i]) < 16:
                logger.debug(f"  {img_i} --> {img_match[i]}")
                logger.debug(f"  {ref_i} --> {ref_match[i]}")
            else:
                logger.debug(
                    f"  {img_i} --> {img_match[i][0:10]} {len(img_match[i]) - 10} more"
                )
                logger.debug(
                    f"  {ref_i} --> {ref_match[i][0:10]} {len(ref_match[i]) - 10} more"
                )
            if i + 1 >= 10 and len(primary_match_img) - 10 > 0:
                logger.debug(
                    f"" f"{(len(primary_match_img) - 10)} additional matches not shown."
                )
                break
        else:
            logger.debug(
                f"{img_i} matches {ref_i} (dPA ={mpa[i]:.3f}):"
                f" {str(len(img_match[i])).strip()} rays"
            )

    if write_crosscheck_files:
        match_lines_im = os.path.splitext(base_output_path)[0] + ".matchlines.im.reg"

        logger.info(f"Writing match lines to {match_lines_im}")

        with open(match_lines_im, "w", encoding="utf8") as out:
            color = "red"
            out.write(
                f"# Region file format: DS9 version 4.0\nglobal color={color} "
                f'font="helvetica 10 normal" select=1 highlite=1 '
                f"edit=1 move=1 delete=1 include=1 fixed=0 source\n"
            )
            out.write("image\n")
            for i, img_i in enumerate(primary_match_img):
                for j, img_j in enumerate(img_match[i]):
                    out.write(
                        f"line({img_src_list[img_i].x:.3f},"
                        f"{img_src_list[img_i].y:.3f},"
                        f"{img_src_list[img_j].x:.3f},"
                        f"{img_src_list[img_j].y:.3f}) # line=0 0\n"
                    )

        match_lines_wcs = os.path.splitext(base_output_path)[0] + ".matchlines.wcs.reg"

        logger.info(f"Writing match lines to {match_lines_wcs}")

        with open(match_lines_wcs, "w", encoding="utf8") as out:
            color = "green"
            out.write(
                f"# Region file format: DS9 version 4.0\nglobal color={color} "
                f'font="helvetica 10 normal" select=1 highlite=1 '
                f"edit=1 move=1 delete=1 include=1 fixed=0 source\n"
            )
            out.write("fk5\n")
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


def crosscheck_source_lists(
    img_src_list: list[SextractorSource],
    n_img: int,
    img_density: float,
    ref_src_list: list[BaseSource],
    n_ref: int,
    ref_density: float,
    box_size_arcsec: float,
    area_sq_min: float,
) -> tuple[list[SextractorSource], int, float, list[BaseSource], int, float]:
    """
    Compares detected sources in image and reference, and trims so they are
    of comparable density

    :param img_src_list: list of image sources
    :param n_img: number of image sources
    :param img_density: img source density
    :param ref_src_list: list of reference images
    :param n_ref: number of reference sources
    :param ref_density: ref source density
    :param box_size_arcsec: radius of search box
    :param area_sq_min: area of image
    :return: trimmed list of sources, number, density,
        trimmed ref sources, number & density
    """
    # If this image is actually shallower than reference catalog, trim the
    # reference catalog down
    if n_ref > 16 and ref_density > 3 * img_density:
        logger.debug("Image is shallow.  Trimming reference catalog...")
        while ref_density > 3 * img_density:
            ref_src_list = ref_src_list[0 : int(len(ref_src_list) * 4 / 5)]
            n_ref = len(ref_src_list)
            ref_density = n_ref / (2 * box_size_arcsec / 60.0) ** 2

    # If the image is way deeper than USNO, trim the image catalog down
    if n_img > 16 and img_density > 4 * ref_density:
        logger.debug("Image is deep.  Trimming image catalog...")
        while img_density > 4 * ref_density and n_img > 8:
            img_src_list = img_src_list[0 : int(len(img_src_list) * 4 / 5)]
            n_img = len(img_src_list)
            img_density = n_img / area_sq_min

    # If too many objects, do some more trimming
    if n_img * n_ref > 120 * 120 * 4:
        logger.debug("Image and/or catalog still too deep.  Trimming...")
        while n_img * n_ref > 120 * 120 * 4:
            if img_density > ref_density:
                img_src_list = img_src_list[0 : int(len(img_src_list) * 4 / 5)]
                n_img = len(img_src_list)
                img_density = n_img / area_sq_min
            else:
                ref_src_list = ref_src_list[0 : int(len(ref_src_list) * 4 / 5)]
                n_ref = len(ref_src_list)
                ref_density = n_ref / (2 * box_size_arcsec / 60.0) ** 2

    # Remove fainter object in close pairs for both lists
    min_sep = 3

    delete_list = []
    for i, src in enumerate(img_src_list):
        for j, src2 in enumerate(img_src_list[i + 1 :]):
            dist = distance(src, src2)
            if dist < min_sep:
                if src.mag > src2.mag:
                    delete_list.append(i)
                else:
                    delete_list.append(i + j + 1)

    delete_list = unique(delete_list)
    delete_list.reverse()
    for del_index in delete_list:
        del img_src_list[del_index]

    for i, src in enumerate(ref_src_list):
        for j, src2 in enumerate(ref_src_list[i + 1 :]):
            dist = distance(src, src2)
            if dist < min_sep:
                if src.mag > src2.mag:
                    delete_list.append(i)
                else:
                    delete_list.append(i + j + 1)

    delete_list = unique(delete_list)
    delete_list.reverse()
    for del_index in delete_list:
        del ref_src_list[del_index]

    return img_src_list, n_img, img_density, ref_src_list, n_ref, ref_density
