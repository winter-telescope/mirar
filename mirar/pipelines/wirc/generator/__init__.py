"""
Dynamically generation functions for WIRC pipeline
"""

from mirar.pipelines.wirc.generator.calibration import (
    wirc_astrometric_catalog_generator,
    wirc_photometric_catalog_generator,
)
from mirar.pipelines.wirc.generator.imsub import (
    wirc_source_table_filter_annotator,
    wirc_zogy_catalogs_purifier,
)
from mirar.pipelines.wirc.generator.references import (
    wirc_reference_generator,
    wirc_reference_image_resampler,
    wirc_reference_psfex,
    wirc_reference_sextractor,
)
from mirar.pipelines.wirc.generator.target import annotate_target_coordinates
