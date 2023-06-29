"""
Module for WINTER data reduction
"""
from mirar.paths import FITS_MASK_KEY
from mirar.pipelines.winter.config import (
    sextractor_anet_config,
    sextractor_autoastrometry_config,
    sextractor_photometry_config,
    swarp_config_path,
)
from mirar.pipelines.winter.generator import (
    scamp_config_path,
    winter_astrometric_catalog_generator,
    winter_mask_path,
    winter_photometric_catalog_generator,
    winter_reference_generator,
)
from mirar.pipelines.winter.load_winter_image import (
    load_proc_winter_image,
    load_raw_winter_image,
    load_stacked_winter_image,
)
from mirar.processors.astromatic import Scamp
from mirar.processors.astromatic.sextractor.sextractor import (
    Sextractor,
    sextractor_checkimg_map,
)
from mirar.processors.astromatic.swarp.swarp import Swarp
from mirar.processors.astrometry.anet.anet_processor import AstrometryNet
from mirar.processors.csvlog import CSVLog
from mirar.processors.dark import DarkCalibrator
from mirar.processors.flat import FlatCalibrator
from mirar.processors.mask import MaskPixelsFromPath, WriteMaskedCoordsToFile
from mirar.processors.photcal import PhotCalibrator
from mirar.processors.reference import GetReferenceImage
from mirar.processors.sky import NightSkyMedianCalibrator, SkyFlatCalibrator
from mirar.processors.split import SplitImage
from mirar.processors.utils import (
    ImageBatcher,
    ImageDebatcher,
    ImageLoader,
    ImageSaver,
    ImageSelector,
)
from mirar.processors.utils.multi_ext_parser import MultiExtParser

refbuild = [
    ImageDebatcher(),
    GetReferenceImage(
        ref_image_generator=winter_reference_generator,
    ),
    ImageSaver(output_dir_name="stacked_ref"),
]

board_id = 4
target_name = "EMGW_gal1"
split = [
    MultiExtParser(
        input_sub_dir="raw/",
        extension_num_header_key="BOARD_ID",
        only_extract_num=board_id,
        output_sub_dir=f"raw_split_{board_id}",
    )
]

split_all_boards = [
    MultiExtParser(
        input_sub_dir="raw/",
        extension_num_header_key="BOARD_ID",
        output_sub_dir="raw_split",
    )
]
load = [
    ImageLoader(
        input_sub_dir=f"raw_split_{board_id}", load_image=load_raw_winter_image
    ),
    ImageSelector(("OBSTYPE", ["FOCUS", "DARK", "FLAT", "SCIENCE"])),
]

load_all_boards = [
    ImageLoader(input_sub_dir="raw_split", load_image=load_raw_winter_image),
    ImageSelector(("OBSTYPE", ["FOCUS", "DARK", "FLAT", "SCIENCE"])),
]

load_proc = [
    ImageLoader(input_sub_dir=f"skysub_{board_id}", load_image=load_raw_winter_image),
    ImageSelector(("TARGNAME", f"{target_name}"), ("OBSTYPE", "SCIENCE")),
]

load_dark = [
    ImageLoader(input_sub_dir=f"darkcal_{board_id}", load_image=load_raw_winter_image),
    ImageSelector(("TARGNAME", f"{target_name}")),
]

load_anet = [
    ImageLoader(input_sub_dir=f"anet_{board_id}", load_image=load_proc_winter_image),
    ImageSelector(("TARGNAME", f"{target_name}"), ("OBSTYPE", "SCIENCE")),
]

load_stack = [
    ImageLoader(input_sub_dir=f"anet_{board_id}", load_image=load_proc_winter_image),
    ImageSelector(("TARGNAME", f"{target_name}"), ("OBSTYPE", "SCIENCE")),
    ImageBatcher("EXPTIME"),
]

load_multiboard_stack = [
    ImageLoader(
        input_sub_dir=f"stack_all_{target_name}", load_image=load_stacked_winter_image
    ),
    ImageSelector(
        ("TARGNAME", f"{target_name}"),
        ("OBSTYPE", "SCIENCE"),
    ),
]
# ImageBatcher("COADDS")]
log = (
    split
    + load
    + [
        CSVLog(
            export_keys=[
                "FILTER",
                "UTCTIME",
                "EXPTIME",
                "OBSTYPE",
                "UNIQTYPE",
                "BOARD_ID",
                "OBSCLASS",
                "TARGET",
                "FILTER",
                "BASENAME",
                "TARGNAME",
                "RADEG",
                "DECDEG",
                "MEDCOUNT",
                "STDDEV",
                "T_ROIC",
            ]
        )
    ]
)

log_all_boards = (
    split_all_boards
    + load_all_boards
    + [
        CSVLog(
            export_keys=[
                "FILTER",
                "UTCTIME",
                "EXPTIME",
                "OBSTYPE",
                "UNIQTYPE",
                "BOARD_ID",
                "OBSCLASS",
                "TARGET",
                "FILTER",
                "BASENAME",
                "TARGNAME",
                "RADEG",
                "DECDEG",
                "MEDCOUNT",
                "STDDEV",
                "T_ROIC",
            ]
        )
    ]
)

dark_cal = [
    ImageSelector(("BOARD_ID", f"{board_id}")),
    ImageBatcher(["BOARD_ID", "EXPTIME"]),
    WriteMaskedCoordsToFile(output_dir="mask_raw"),
    DarkCalibrator(cache_sub_dir=f"calibration_{board_id}"),
    ImageSaver(output_dir_name=f"darkcal_{board_id}"),
    ImageDebatcher(),
]

dark_cal_all_boards = [
    ImageBatcher(["BOARD_ID", "EXPTIME"]),
    WriteMaskedCoordsToFile(output_dir="mask_raw"),
    DarkCalibrator(cache_sub_dir="calibration"),
    ImageSaver(output_dir_name="darkcal"),
    ImageDebatcher(),
]

flat_cal = [
    # ImageSelector(("OBSTYPE", ["FOCUS", "SCIENCE", "FLAT"])),
    # ImageSelector(("TARGNAME", ["INTERESTING"])),
    # ImageBatcher(["BOARD_ID", "FILTER"]),
    # FlatCalibrator(flat_mask_key=FITS_MASK_KEY,
    #                cache_sub_dir=f"calibration_{board_id}"
    #                ),
    ImageSelector(("OBSTYPE", ["SCIENCE"]), ("TARGNAME", f"{target_name}")),
    ImageBatcher(["BOARD_ID", "FILTER", "TARGNAME", "EXPTIME"]),
    SkyFlatCalibrator(flat_mask_key=FITS_MASK_KEY, cache_sub_dir=f"skycals_{board_id}"),
    ImageSelector(("OBSTYPE", ["SCIENCE"])),
    ImageSaver(output_dir_name=f"skyflatcal_{board_id}"),
    # ImageSelector(("OBSTYPE", ["SCIENCE"])),
    # Sextractor(**sextractor_astrometry_config,
    #            write_regions_bool=True,
    #            output_sub_dir="sextractor",
    #            cache=True),
    # ImageSelector(("TARGNAME", [""])),
    # ImageBatcher(["BOARD_ID", "FILTER", "TARGNAME", "EXPTIME"]),
    NightSkyMedianCalibrator(flat_mask_key=FITS_MASK_KEY),
    ImageSaver(output_dir_name=f"skysub_{board_id}"),
    # Sextractor(**sextractor_astrometry_config,
    #            write_regions_bool=True,
    #            output_sub_dir="sextractor",
    #            cache=True),
    # AutoAstrometry(catalog="tmc", pixel_scale=1.0, pa=0, inv=True,
    #                write_crosscheck_files=True),
    # ImageSaver(output_dir_name="skysub"),
    # Sextractor(**sextractor_astrometry_config,
    #            write_regions_bool=True,
    #            output_sub_dir="sextractor"),
    # MultiExtParser(input_sub_dir="raw/mef/"),
    # SplitImage(),
    # MaskPixelsFromPath(mask_path=winter_mask_path),
    # DarkCalibrator(),
    # SkyFlatCalibrator(),
    # NightSkyMedianCalibrator(),
]

flat_cal_all_boards = [
    ImageSelector(("OBSTYPE", ["SCIENCE"]), ("TARGNAME", f"{target_name}")),
    ImageBatcher(["BOARD_ID", "FILTER", "TARGNAME", "EXPTIME"]),
    SkyFlatCalibrator(flat_mask_key=FITS_MASK_KEY, cache_sub_dir="skycals"),
    ImageSelector(("OBSTYPE", ["SCIENCE"])),
    ImageSaver(output_dir_name="skyflatcal"),
    NightSkyMedianCalibrator(flat_mask_key=FITS_MASK_KEY),
    ImageSaver(output_dir_name="skysub"),
]

process = dark_cal + flat_cal
process_proc = [
    ImageDebatcher(),
    AstrometryNet(
        output_sub_dir=f"anet_{board_id}",
        scale_bounds=[25, 40],
        scale_units="amw",
        use_sextractor=True,
        parity="neg",
        search_radius_deg=1.0,
        # sextractor_config_path=sextractor_autoastrometry_config[
        #     'config_path'],
        use_weight=False,
    ),
    ImageSaver(output_dir_name=f"anet_{board_id}", use_existing_weight=False),
    Sextractor(
        **sextractor_autoastrometry_config,
        write_regions_bool=True,
        output_sub_dir="scamp",
    ),
    Scamp(
        temp_output_sub_dir="scamp",
        ref_catalog_generator=winter_astrometric_catalog_generator,
        scamp_config_path=scamp_config_path,
        cache=False,
    ),
    ImageDebatcher(),
    Swarp(
        swarp_config_path=swarp_config_path,
        calculate_dims_in_swarp=True,
        include_scamp=True,
        subtract_bkg=False,
        cache=True,
        center_type="ALL",
        temp_output_sub_dir=f"stack_all_{target_name}",
    ),
]

process_proc_all_boards = [
    ImageDebatcher(),
    ImageBatcher(["UTCTIME", "BOARD_ID"]),
    AstrometryNet(
        output_sub_dir="anet",
        scale_bounds=[25, 40],
        scale_units="amw",
        use_sextractor=True,
        parity="neg",
        search_radius_deg=1.0,
        sextractor_config_path=sextractor_anet_config[
            'config_path'],
        use_weight=False,
    ),
    ImageSaver(output_dir_name=f"anet", use_existing_weight=False),
    Sextractor(
        **sextractor_autoastrometry_config,
        write_regions_bool=True,
        output_sub_dir="scamp",
    ),
    Scamp(
        temp_output_sub_dir="scamp",
        ref_catalog_generator=winter_astrometric_catalog_generator,
        scamp_config_path=scamp_config_path,
        cache=False,
    ),
    ImageDebatcher(),
    ImageBatcher(["BOARD_ID", "FILTER", "TARGNAME"]),
    # ImageSaver(output_dir_name="pre-swarp"),
    Swarp(
        swarp_config_path=swarp_config_path,
        calculate_dims_in_swarp=True,
        include_scamp=True,
        subtract_bkg=False,
        cache=False,
        center_type="ALL",
        temp_output_sub_dir=f"stack_all_{target_name}",
    ),
]

process_noise = [
    Sextractor(
        **sextractor_autoastrometry_config,
        write_regions_bool=True,
        output_sub_dir="mask_sextractor",
        checkimage_type="SEGMENTATION",
        cache=True,
        verbose_type="FULL",
    ),
    MaskPixelsFromPath(
        mask_path_key=sextractor_checkimg_map["SEGMENTATION"],
        write_masked_pixels_to_file=True,
        output_dir="mask1",
    ),
    ImageSaver(output_dir_name=f"noisemask_{board_id}"),
]
# process_proc = [ImageDebatcher(),
#                 ImageSelector(("OBSTYPE", ["SCIENCE"])),
#                 Sextractor(**sextractor_autoastrometry_config,
#                            write_regions_bool=True,
#                            output_sub_dir="sextractor",
#                            cache=False),
#                 AutoAstrometry(catalog="tmc", pixel_scale=1.07,
#                                write_crosscheck_files=True,
#                                inv=False,
#                                ), ]

stack_proc = [
    ImageBatcher(["TARGNAME", "FILTER"]),
    Swarp(
        swarp_config_path=swarp_config_path,
        calculate_dims_in_swarp=True,
        include_scamp=False,
        subtract_bkg=True,
        cache=True,
        center_type="MANUAL",
    ),
    ImageSaver(output_dir_name=f"stack_{target_name}"),
]

photcal = [
    # ImageSelector(("BOARD_ID", board_id)),
    ImageDebatcher(),
    ImageBatcher(["BOARD_ID"]),
    Sextractor(
        **sextractor_photometry_config,
        output_sub_dir=f"phot_{board_id}_{target_name}",
        checkimage_type="BACKGROUND_RMS",
    ),
    PhotCalibrator(
        ref_catalog_generator=winter_photometric_catalog_generator,
        temp_output_sub_dir=f"phot_{board_id}_{target_name}",
        write_regions=True,
        cache=True,
    ),
    # ImageSaver(output_dir_name=f"phot_{board_id}_{target_name}")
    ImageSaver(output_dir_name=f"phot_{target_name}"),
]

photcal_indiv = [
    ImageSelector(("BOARD_ID", board_id)),
    ImageDebatcher(),
    ImageBatcher(["UTCTIME"]),
    Sextractor(
        **sextractor_photometry_config,
        output_sub_dir=f"phot_{board_id}_{target_name}",
        checkimage_type="BACKGROUND_RMS",
    ),
    PhotCalibrator(
        ref_catalog_generator=winter_photometric_catalog_generator,
        temp_output_sub_dir=f"phot_{board_id}_{target_name}",
        write_regions=True,
        cache=True,
    ),
    ImageSaver(output_dir_name=f"phot_{board_id}_{target_name}"),
]

stack_multiboard = [
    Swarp(
        swarp_config_path=swarp_config_path,
        calculate_dims_in_swarp=True,
        include_scamp=False,
        subtract_bkg=False,
        cache=True,
        temp_output_sub_dir=f"multiboard_stack_{target_name}",
        center_type="MANUAL",
    )
]
# commissioning = \
#     log + process_proc

commissioning = log + process

commissioning_dark = log + dark_cal
commissioning_proc = load_proc + process_proc
commissioning_flat = load_dark + flat_cal
commissioning_reduce = log + dark_cal + flat_cal
commissioning_stack = load_stack + stack_proc
commissioning_multiboard_stack = load_multiboard_stack + stack_multiboard
commissioning_noise = load_anet + process_noise
commissioning_photcal = load_multiboard_stack + photcal
commissioning_photcal_indiv = load_anet + photcal_indiv
full_commissioning = log + process + process_proc  # + stack_proc
full_commissioning_all_boards = (
    log_all_boards + dark_cal_all_boards + flat_cal_all_boards + process_proc_all_boards
)
