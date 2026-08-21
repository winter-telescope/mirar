"""
Module containing a processor to run autoastrometry.
"""

import logging
from pathlib import Path
from typing import Optional

from mirar.data import ImageBatch
from mirar.paths import BASE_NAME_KEY, get_output_dir
from mirar.processors.astromatic.sextractor.sourceextractor import DEFAULT_SATURATION
from mirar.processors.astrometry.autoastrometry.autoastrometry import (
    DEFAULT_TOLERANCE,
    run_autoastrometry_single,
)
from mirar.processors.astrometry.autoastrometry.errors import (
    AstrometryCrossmatchError,
    AstrometryReferenceError,
    AstrometrySourceError,
    AstrometryURLError,
)
from mirar.processors.astrometry.autoastrometry.utils import (
    dec_str_2_deg,
    ra_str_2_deg,
)
from mirar.processors.base_processor import BaseImageProcessor

logger = logging.getLogger(__name__)

ASTROMETRY_ERRORS = (
    AstrometryCrossmatchError,
    AstrometryReferenceError,
    AstrometrySourceError,
    AstrometryURLError,
)

# Header keywords searched, in order, for the telescope pointing
RA_KEYS = ["RA", "OBJRA", "TELRA", "CRVAL1"]
DEC_KEYS = ["DEC", "OBJDEC", "TELDEC", "CRVAL2"]


class AutoAstrometry(BaseImageProcessor):
    """Processor to run autoastrometry"""

    base_key = "autoastrometry"

    def __init__(
        self,
        temp_output_sub_dir: str = "autoastrometry",
        write_crosscheck_files: bool = False,
        catalog: Optional[str] = None,
        pixel_scale: Optional[float] = None,
        inv: bool = False,
        pa: Optional[float] = None,
        seeing: Optional[float] = None,
        box_size: Optional[float] = None,
        max_rad: Optional[float] = None,
        tolerance: float = DEFAULT_TOLERANCE,
        unc_pa: Optional[float] = None,
        max_ellip: float = 0.5,
        saturation: float = DEFAULT_SATURATION,
        no_rot: bool = False,
        sextractor_config_path: Optional[str | Path] = None,
        use_existing_catalog: bool = False,
        use_header_radec: bool = True,
        max_scatter_arcsec: Optional[float] = None,
        min_matches: Optional[int] = None,
        on_failure: str = "raise",
    ):
        """
        :param temp_output_sub_dir: subdirectory for temp files
        :param write_crosscheck_files: write .reg / .txt crosscheck files
        :param catalog: 'tmc', 'ub2' or 'sdss' (None -> try all)
        :param pixel_scale: arcsec/pix. Setting this DISCARDS the existing WCS
            and builds a provisional one from scratch. Must be good to ~1%.
        :param inv: reverse (positive) parity. Only has an effect when
            pixel_scale is also given.
        :param pa: position angle in deg. Defaults to 0 when pixel_scale is set.
        :param seeing: approximate seeing in PIXELS (not arcsec). Sets the
            FWHM acceptance window to [0.7*seeing, 2*seeing].
        :param box_size: half-width of the reference catalogue query, arcsec.
            None -> the image field width.
        :param max_rad: max separation for star pairs, arcsec. None -> auto.
        :param tolerance: fractional slack in the pair-distance match.
        :param unc_pa: PA uncertainty in deg; restricts the rotation search.
        :param max_ellip: max source ellipticity to keep.
        :param saturation: SATUR_LEVEL passed to SExtractor.
        :param no_rot: skip the rotation correction, shift only.
        :param sextractor_config_path: config for the internal SExtractor run.
            Must have CATALOG_TYPE ASCII_HEAD. None -> mirar's default.
        :param use_existing_catalog: reuse a catalogue from an upstream
            Sextractor step via the SRCCAT header key. Leave False whenever
            pixel_scale/pa/inv are set - the upstream RA/Dec were computed with
            the old WCS and will not match the provisional one.
        :param use_header_radec: read the pointing from RA/DEC-style keywords
            rather than relying on CRVAL1/CRVAL2 existing.
        :param max_scatter_arcsec: fail if ASTR_UNC exceeds this.
        :param min_matches: fail if ASTR_NUM falls below this.
        :param on_failure: 'raise' or 'flag'. 'flag' sets ASTRGOOD=False and
            lets the image continue, for filtering downstream.
        """
        super().__init__()

        if on_failure not in ["raise", "flag"]:
            raise ValueError(f"on_failure must be 'raise' or 'flag', got {on_failure}")

        self.temp_output_sub_dir = temp_output_sub_dir
        self.write_crosscheck_files = write_crosscheck_files
        self.catalog = catalog
        self.pixel_scale = pixel_scale
        self.inv = inv
        self.pa = pa
        self.seeing = seeing
        self.box_size = box_size
        self.max_rad = max_rad
        self.tolerance = tolerance
        self.unc_pa = unc_pa
        self.max_ellip = max_ellip
        self.saturation = saturation
        self.no_rot = no_rot
        self.sextractor_config_path = sextractor_config_path
        self.use_existing_catalog = use_existing_catalog
        self.use_header_radec = use_header_radec
        self.max_scatter_arcsec = max_scatter_arcsec
        self.min_matches = min_matches
        self.on_failure = on_failure

    def description(self) -> str:
        return "Processor to perform astrometric calibration."

    @staticmethod
    def get_pointing(image) -> tuple[Optional[float], Optional[float]]:
        """
        Pull the telescope pointing out of the header, tolerating both
        sexagesimal and decimal-degree conventions.

        :param image: image to inspect
        :return: (ra_deg, dec_deg), either may be None
        """
        ra_deg, dec_deg = None, None

        for key in RA_KEYS:
            if key in image.keys():
                try:
                    ra_deg = ra_str_2_deg(image[key])
                    break
                except (ValueError, TypeError):
                    continue

        for key in DEC_KEYS:
            if key in image.keys():
                try:
                    dec_deg = dec_str_2_deg(image[key])
                    break
                except (ValueError, TypeError):
                    continue

        return ra_deg, dec_deg

    def _validate(self, image, base_name: str) -> tuple[bool, str]:
        """
        Check the ASTR_* keywords written by autoastrometry.

        :param image: solved image
        :param base_name: name, for the log message
        :return: (passed, reason)
        """
        n_match = image.get("ASTR_NUM", 0)
        scatter = image.get("ASTR_UNC", 999.0)

        if self.min_matches is not None and n_match < self.min_matches:
            return False, (
                f"{base_name}: only {n_match} matches "
                f"(min_matches={self.min_matches})"
            )

        if self.max_scatter_arcsec is not None and scatter > self.max_scatter_arcsec:
            return False, (
                f'{base_name}: astrometric scatter {scatter:.3f}" '
                f'(max={self.max_scatter_arcsec}")'
            )

        return True, ""

    def _apply_to_images(self, batch: ImageBatch) -> ImageBatch:
        output_dir = get_output_dir(self.temp_output_sub_dir, self.night_sub_dir)
        output_dir.mkdir(parents=True, exist_ok=True)

        for i, image in enumerate(batch):
            base_name = image[BASE_NAME_KEY]
            temp_path = output_dir.joinpath(base_name)
            self.save_fits(image, temp_path, compress=False)

            user_ra, user_dec = None, None
            if self.use_header_radec:
                user_ra, user_dec = self.get_pointing(image)
                logger.debug(f"{base_name}: pointing {user_ra}, {user_dec}")

            failed_reason = ""

            try:
                fit_info = run_autoastrometry_single(
                    img_path=temp_path,
                    output_dir=output_dir,
                    write_crosscheck_files=self.write_crosscheck_files,
                    overwrite=True,
                    catalog=self.catalog,
                    pixel_scale=self.pixel_scale,
                    inv=self.inv,
                    pa=self.pa,
                    seeing=self.seeing,
                    box_size=self.box_size,
                    max_rad=self.max_rad,
                    tolerance=self.tolerance,
                    uncpa=self.unc_pa,
                    max_ellip=self.max_ellip,
                    saturation=self.saturation,
                    no_rot=self.no_rot,
                    user_ra=user_ra,
                    user_dec=user_dec,
                    sextractor_config_path=self.sextractor_config_path,
                    use_existing_catalog=self.use_existing_catalog,
                )
                logger.debug(f"{base_name}: fit_info {fit_info}")

            except ASTROMETRY_ERRORS as exc:
                failed_reason = f"{base_name}: autoastrometry failed - {exc}"
                logger.error(failed_reason)

                if self.on_failure == "raise":
                    temp_path.unlink(missing_ok=True)
                    raise

            # Load the updated header back, then clean up
            image = self.open_fits(temp_path)
            temp_path.unlink(missing_ok=True)

            if failed_reason == "":
                passed, failed_reason = self._validate(image, base_name)
                if not passed:
                    logger.error(failed_reason)
                    if self.on_failure == "raise":
                        raise AstrometryCrossmatchError(failed_reason)

            image["ASTRGOOD"] = failed_reason == ""
            if failed_reason != "":
                image["ASTRFAIL"] = failed_reason[:68]

            batch[i] = image

        return batch