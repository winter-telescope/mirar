"""
Module to calculate zero point using outlier rejection
"""

from mirar.processors.photcal.zp_calculator.base_zp_calculator import (
    BaseZeroPointCalculator,
)
from mirar.processors.photcal.zp_calculator.color_term_zp_calculator import (
    ZPWithColorTermCalculator,
)
from mirar.processors.photcal.zp_calculator.outlier_rejection_zp_calculator import (
    OutlierRejectionZPCalculator,
)
