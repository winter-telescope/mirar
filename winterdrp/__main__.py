#!/usr/bin/env python
import argparse
import sys
import logging
from winterdrp.pipelines import get_pipeline

logger = logging.getLogger(__name__)

parser = argparse.ArgumentParser(
    description="winterdrp: An automated image reduction pipeline, developed for WINTER"
)
# parser.add_argument(
#     '-d',
#     "--datadir",
#     help='Data directory : Specify the name of the directory with raw files'
# )
# parser.add_argument(
#     '--o',
#     default=None,
#     help='Object name : Give the name of the object to be reduced.'
# )
# parser.add_argument(
#     '--f',
#     default=None,
#     help='Filter name: Specify filter to be reduced (options are J, H and Ks).'
# )
# parser.add_argument(
#     '--g',
#     help='Extended: Look for frames with extended emission, and do not use them for sky background.',
#     action='store_true',
#     default=False
# )
# parser.add_argument(
#     '--dk',
#     default=None,
#     help='Dark file list : Give a list of files with darks for the images to be reduced.'
# )
# parser.add_argument(
#     '--so',
#     default=None,
#     help='Source file list : Give list of files with source images to be reduced.'
#          ' Must be combined with corresponding sky file list.'
# )
# parser.add_argument(
#     '--sk',
#     default=None,
#     help='Sky file list: Give list of files corresponding to the source images. '
#          'Must be combined with corresponding source file list.'
# )
# parser.add_argument(
#     '--m',
#     default=None,
#     help='Mode specific execution: Run only specific steps of the reduction pipeline. '
#          'Options are fp (FIRST-PASS), sp (SECOND-PASS), as (ASTROMETRY AND STACKING), '
#          'pc (PHOTOMETRIC CALIBRATION). Join multiple modes with with spaceless commas, '
#          'e.g "--m fp,sp". Note that running each step expects that the previous step '
#          'has already been run, and expects the associated data products in the '
#          'reduction directory (except first-pass, which is the first step anyways)'
# )
# parser.add_argument(
#     '--c',
#     help='Add this option along with -m if you want a specific step '
#          'and all subsequent steps to be executed',
#     action='store_true',
#     default=False
# )
parser.add_argument(
    "-n",
    "--night",
    default="",
    help="Sub-directory to use in the data directory"
)
parser.add_argument(
    "-p",
    "--pipeline",
    default="summer",
    help="Pipeline to be used"
)
parser.add_argument(
    "-c",
    "--config",
    default=None,
    help="Pipeline configuration to be used"
)
parser.add_argument(
    "--logfile",
    default=None,
    help="If a path is passed, all logs will be written to this file."
)
parser.add_argument(
    "--level",
    default="DEBUG",
    help="Python logging level"
)
parser.add_argument(
    "-b",
    "--batch",
    default=None,
    help="Only process a specific image batch"
)
parser.add_argument(
    '--skipcache',
    help='Skip building cache',
    action='store_true',
    default=False
)
# parser.add_argument(
#     '-skipfail',
#     help='If processing of one image set fails, proceed with other objects/filters',
#     action='store_true',
#     default=False
# )

args = parser.parse_args()

# Set up logging

log = logging.getLogger("winterdrp")

if args.logfile is None:
    handler = logging.StreamHandler(sys.stdout)
else:
    handler = logging.FileHandler(args.logfile)

formatter = logging.Formatter('%(name)s [l %(lineno)d] - %(levelname)s - %(message)s')
handler.setFormatter(formatter)
log.addHandler(handler)
log.setLevel(args.level)

# Parse running modes

# if args.datadir is None:
#     err = "No data directory specified. Please rerun with '--da /path/to/data'"
#     logger.error(err)
#     raise Exception(err)

# if args.m is None:
#     modes = wp.all_modes
# else:
#     modes = [x for x in wp.all_modes if x in args.m.split(",")]
#
# if args.c:
#
#     if len(modes) > 1:
#         warning = f"Running with '-c' option (continue), so running all steps after '{modes[0]}'. " \
#                   f"However, multiple modes were specified ({modes}), so the command is ambiguous."
#         logger.warning(warning)
#
#     modes = wp.all_modes[wp.all_modes.index(modes[0]):]

# logger.info(f"Running modes selected: {modes}")

pipe = get_pipeline(
    args.pipeline,
    configuration=args.config,
    night=args.night,
    skip_build_cache=args.skipcache
)
# pipe.make_calibration_files(sub_dir=args.subdir)
image_batches = pipe.split_raw_images_into_batches(
    select_batch=args.batch
)

for image_batch in image_batches:
    images, headers = pipe.open_image_batch(image_batch)
    pipe.reduce_images(images, headers)

logger.info('END OF WIRC-PIPE EXECUTION')
