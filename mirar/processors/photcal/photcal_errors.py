"""
Photometric calibration errors
"""

from mirar.errors import ProcessorError


class PhotometryError(ProcessorError):
    """Base error for photometric calibration"""


class PhotometryReferenceError(PhotometryError):
    """Error related to the photometric reference catalogue"""


class PhotometrySourceError(PhotometryError):
    """Error related to the photometric source catalogue"""


class PhotometryCrossMatchError(PhotometryError):
    """Error related to cross-matching photometric reference and source catalogues"""


class PhotometryCalculationError(PhotometryError):
    """Error related to the photometric calibration"""
