"""
Module containing custom errors for autoastrometry
"""
from mirar.errors import ProcessorError


class AstrometryError(ProcessorError):
    """Parent Error relating to autoastrometry"""


class AstrometrySourceError(AstrometryError):
    """Error with astrometry source catalogue"""


class AstrometryURLError(AstrometryError):
    """Error related to astrometry URL query"""


class AstrometryReferenceError(AstrometryError):
    """Error related to astrometry reference catalogue"""


class AstrometryCrossmatchError(AstrometryError):
    """Error related to crossmatching astrometry source/reference catalogue"""
