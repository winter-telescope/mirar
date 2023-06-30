"""
Component of autoastrometry which detects sources in the science image
"""
import logging
import os
from pathlib import Path
from typing import Optional

from astropy.io import fits

from mirar.paths import SEXTRACTOR_HEADER_KEY
from mirar.processors.astromatic.sextractor.settings import (
    default_config_path,
    default_conv_path,
    default_param_path,
    default_starnnw_path,
)
from mirar.processors.astromatic.sextractor.sourceextractor import (
    DEFAULT_SATURATION,
    run_sextractor_single,
)
from mirar.processors.astrometry.autoastrometry.errors import AstrometrySourceError
from mirar.processors.astrometry.autoastrometry.sources import (
    SextractorSource,
    compare_mag,
)
from mirar.processors.astrometry.autoastrometry.utils import median, mode
from mirar.utils.ldac_tools import get_table_from_ldac

logger = logging.getLogger(__name__)

DEFAULT_MIN_FWHM = 1.5
DEFAULT_MAX_FWHM = 40


def get_img_src_list(
    img_path: str,
    base_output_path: str,
    nx_pix: int,
    ny_pix: int,
    border: int = 3,
    corner: int = 12,
    min_fwhm: float = DEFAULT_MIN_FWHM,
    max_fwhm: float = DEFAULT_MAX_FWHM,
    max_ellip: float = 0.5,
    saturation: float = DEFAULT_SATURATION,
    config_path: str = default_config_path,
    output_catalog: Optional[str | Path] = None,
    write_crosscheck_files: bool = False,
) -> list[SextractorSource]:
    """
    Run sextractor on an image, and then extract all the sources

    :param img_path: image to evaluate
    :param base_output_path:
    :param nx_pix: number of x pixels
    :param ny_pix: number of y pixels
    :param border: number of border pixels
    :param corner: number of corner pixels
    :param min_fwhm: min FWHM to allow
    :param max_fwhm: max FWHM to allow
    :param max_ellip: max ellipse
    :param saturation: saturation level
    :param config_path: sextractor config path
    :param output_catalog: output path
    :param write_crosscheck_files: boolean to write additional crosscheck files
    :return: list of sextractor sources
    """

    if output_catalog is None:
        output_catalog = Path(base_output_path).with_suffix(".cat")

    try:
        os.remove(output_catalog)
    except FileNotFoundError:
        pass

    header = fits.getheader(img_path)
    sextractor_catalog_path = None
    if SEXTRACTOR_HEADER_KEY in header.keys():
        sextractor_catalog_path = fits.getval(img_path, SEXTRACTOR_HEADER_KEY)

    if sextractor_catalog_path is not None:
        logger.info("Using existing sextractor catalog")
        output_catalog = sextractor_catalog_path
        output_catalog_table = get_table_from_ldac(output_catalog)
        output_catalog_table.write(
            output_catalog.replace(".cat", ".ascii.cat"),
            format="ascii.ecsv",
            overwrite=True,
        )
        output_catalog = output_catalog.replace(".cat", ".ascii.cat")
    else:
        run_sextractor_single(
            img=img_path,
            output_dir=os.path.dirname(output_catalog),
            config=config_path,
            saturation=saturation,
            catalog_name=output_catalog,
            parameters_name=default_param_path,
            filter_name=default_conv_path,
            starnnw_name=default_starnnw_path,
        )

    # Read in the sextractor catalog
    with open(output_catalog, "rb") as cat:
        catlines = [x.replace(b"\x00", b"").decode() for x in cat.readlines()][1:]

    # Delete the sextractor catalog again, if not requested
    if not write_crosscheck_files:
        os.remove(output_catalog)

    if len(catlines) == 0:
        logger.error("Sextractor catalog is empty: try a different catalog?")
        raise ValueError

    min_x = border
    min_y = border
    max_x = nx_pix - border  # This should be generalized
    max_y = ny_pix - border

    n_src_init = 0
    n_src_pass = 0
    src_list = []

    rejects = []

    for line in catlines:
        if line[0] == "#":
            continue

        # TODO : make more generic to skip column names
        if line[0] == "X":
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

    if n_src_pass == 0:
        reject_stats = [(x, rejects.count(x)) for x in list(set(rejects))]
        err = (
            f"Found no good sources in image to use for astrometry "
            f"({n_src_init} sources found initially). "
            f"Reject reasons: {reject_stats}"
        )
        logger.error(err)
        raise AstrometrySourceError(err)

    # Remove detections along bad columns

    thresh_prob = 0.0001
    ct_bad_col = 0
    for _ in range(5):
        for variable in ["x", "y"]:
            txp = 1.0
            val_thresh = 1
            while txp > thresh_prob:
                txp *= min(
                    (len(src_list) * 1.0 / nx_pix), 0.8
                )  # some strange way of estimating the threshold.
                val_thresh += (
                    1  # what I really want is a general analytic expression for
                )

            remove_list = []  # the 99.99% prob. threshold for value of n for >=n out
            val_list = [getattr(src, variable) for src in src_list]
            mode_val = mode(
                val_list
            )  # of N total sources to land in the same bin (of NX total bins)
            for j, val in enumerate(val_list):
                if (val > mode_val - 1) & (val < mode_val + 1):
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
        fwhm_20 = fwhm_list[int(len(fwhm_list) / 5)]
        fwhm_mode = mode(fwhm_list)
    else:
        fwhm_mode = min_fwhm
        fwhm_20 = min_fwhm

    # formerly a max, but occasionally a preponderance of long CR's could
    # cause fwhm_mode to be bigger than the stars
    refined_min_fwhm = median([0.75 * fwhm_mode, 0.9 * fwhm_20, min_fwhm])
    # if CR's are bigger and more common than stars, this is dangerous...
    logger.debug(f"Refined min FWHM: {refined_min_fwhm} pix")

    # Might also be good to screen for false detections created by bad columns/rows

    n_good = 0
    good_src_list = []
    for src in src_list:
        if src.fwhm > refined_min_fwhm:
            good_src_list.append(src)
            n_good += 1
        else:
            rejects.append("refined min fwhm")

    if n_good == 0:
        reject_stats = [(x, rejects.count(x)) for x in list(set(rejects))]
        err = (
            f"Found no good sources in image to use for astrometry "
            f"({n_src_init} sources found initially). "
            f"Reject reasons: {reject_stats}"
        )
        logger.error(err)
        raise AstrometrySourceError(err)

    # Sort by magnitude
    good_src_list.sort(key=compare_mag)

    logger.debug(
        f"{n_good} objects detected in image {img_path} "
        f"(a further {n_src_init - n_good} discarded)"
    )

    reject_stats = [(x, rejects.count(x)) for x in list(set(rejects))]
    logger.debug(f"Reject reasons: {reject_stats}")

    return good_src_list
