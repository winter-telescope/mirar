"""
Module containing base Source classes used by autoastrometry
"""
import logging

import numpy as np

logger = logging.getLogger(__name__)


class BaseSource:
    """
    A standard source object
    """

    def __init__(self, ra_deg: float, dec_deg: float, in_mag: float):
        self.ra_deg = float(ra_deg)
        self.dec_deg = dec_deg
        self.ra_rad = ra_deg * np.pi / 180
        self.dec_rad = dec_deg * np.pi / 180
        self.mag = in_mag

    def rotate(self, dpa_deg: float, ra0: float, dec0: float):
        """
        Function to rotate a source by dpa around ra0/dec0?

        :param dpa_deg: delta-pa (deg)
        :param ra0: ra
        :param dec0: dec
        :return: None
        """
        dpa_rad = dpa_deg * np.pi / 180
        sin_dpa = np.sin(dpa_rad)
        cos_dpa = np.cos(dpa_rad)
        ra_scale = np.cos(dec0 * np.pi / 180)

        # this is only valid for small fields away from the pole.
        source_x = (self.ra_deg - ra0) * ra_scale
        source_y = self.dec_deg - dec0

        x_rot = cos_dpa * source_x - sin_dpa * source_y
        y_rot = sin_dpa * source_x + cos_dpa * source_y

        self.ra_deg = (x_rot / ra_scale) + ra0
        self.dec_deg = y_rot + dec0
        self.ra_rad = self.ra_deg * np.pi / 180
        self.dec_rad = self.dec_deg * np.pi / 180


class SextractorSource(BaseSource):
    """Source from sextractor"""

    def __init__(self, line: str):
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
def pixel_distance(obj1: SextractorSource, obj2: SextractorSource) -> float:
    """
    Calculate pixel distance between obj1 and obj2

    :param obj1: Object 1
    :param obj2: Object 2
    :return: pixel distance
    """
    return ((obj1.x - obj2.x) ** 2 + (obj1.y - obj2.y) ** 2) ** 0.5


def distance(obj1: BaseSource, obj2: BaseSource) -> float:
    """
    # Great circle distance between two points.

    :param obj1: object 1
    :param obj2: object 2
    :return: great circle distance
    """
    ddec = obj2.dec_rad - obj1.dec_rad
    dra = obj2.ra_rad - obj1.ra_rad
    dist_rad = 2 * np.arcsin(
        np.sqrt(
            (np.sin(ddec / 2.0)) ** 2
            + np.cos(obj1.dec_rad) * np.cos(obj2.dec_rad) * (np.sin(dra / 2.0)) ** 2
        )
    )

    dist_deg = dist_rad * 180.0 / np.pi
    dist_arc_sec = dist_deg * 3600.0
    return dist_arc_sec


def quickdistance(obj1: BaseSource, obj2: BaseSource, cosdec: float) -> float:
    """
    Cartestian-approximation distance between two objects
    (Non-great-circle distance is much faster, but beware poles...)

    :param obj1: object 1
    :param obj2: object 2
    :param cosdec: cos(declination)
    :return: approximate distance
    """
    ddec = obj2.dec_deg - obj1.dec_deg
    dra = obj2.ra_deg - obj1.ra_deg
    if dra > 180:
        dra = 360 - dra
    return 3600 * np.sqrt(ddec**2 + (cosdec * dra) ** 2)


def position_angle(obj1: BaseSource, obj2: BaseSource) -> float:
    """
    Calculate the (spherical) position angle between two objects.

    :param obj1: Object 1
    :param obj2: Object 2
    :return: angle (degrees)
    """
    dra = obj2.ra_rad - obj1.ra_rad
    pa_rad = np.arctan2(
        np.cos(obj1.dec_rad) * np.tan(obj2.dec_rad)
        - np.sin(obj1.dec_rad) * np.cos(dra),
        np.sin(dra),
    )
    pa_deg = pa_rad * 180.0 / np.pi
    pa_deg = 90.0 - pa_deg  # defined as degrees east of north
    while pa_deg > 200:
        pa_deg -= 360.0  # make single-valued
    while pa_deg < -160:
        pa_deg += (
            360.0  # note there is a crossing point at PA=200, images at this exact PA
        )
    return (
        pa_deg  # will have the number of matches cut by half at each comparison level
    )


def compare_mag(source: SextractorSource) -> float:
    """
    # Compare objects using magnitude.
    "useful for sorting; Altered by KD for compatibility with python 3

    :param source: Source to compare
    :return: magnitude of source
    """
    return source.mag
