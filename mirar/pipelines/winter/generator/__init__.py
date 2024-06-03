"""
Module for generating custom processors and catalogues for WINTER
"""

from mirar.pipelines.winter.generator.astrometry import (
    winter_astrometric_ref_catalog_generator,
    winter_astrometry_sextractor_catalog_purifier,
)
from mirar.pipelines.winter.generator.candidates import (
    winter_candidate_annotator_filterer,
    winter_candidate_avro_fields_calculator,
    winter_candidate_quality_filterer,
    winter_new_source_updater,
    winter_skyportal_annotator,
    winter_source_entry_updater,
)
from mirar.pipelines.winter.generator.photometry import (
    winter_astrostat_catalog_purifier,
    winter_photcal_color_columns_generator,
    winter_photometric_catalog_generator,
    winter_photometric_catalogs_purifier,
    winter_ref_photometric_catalogs_purifier,
    winter_reference_phot_calibrator,
)
from mirar.pipelines.winter.generator.realbogus import apply_rb_to_table
from mirar.pipelines.winter.generator.reduce import (
    mask_stamps_around_bright_stars,
    select_winter_dome_flats_images,
    select_winter_sky_flat_images,
    winter_anet_sextractor_config_path_generator,
    winter_boardid_6_demasker,
    winter_fourier_filtered_image_generator,
    winter_history_deprecated_constraint,
    winter_imsub_catalog_purifier,
    winter_stackid_annotator,
)
from mirar.pipelines.winter.generator.references import (
    winter_astrometric_ref_catalog_namer,
    winter_photometric_ref_catalog_namer,
    winter_reference_generator,
    winter_reference_image_resampler_for_zogy,
    winter_reference_psf_phot_sextractor,
    winter_reference_psfex,
    winter_reference_sextractor,
    winter_reference_stack_annotator,
    winter_wfau_component_image_stacker,
)
