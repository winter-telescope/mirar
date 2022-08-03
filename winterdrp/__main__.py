#!/usr/bin/env python
import argparse
import os
import sys
import logging
from winterdrp.pipelines import get_pipeline, Pipeline
from winterdrp.paths import raw_img_dir

logger = logging.getLogger(__name__)

parser = argparse.ArgumentParser(
    description="winterdrp: An automated image reduction pipeline, developed for WINTER"
)
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
    '--download',
    help='Download images from server',
    action='store_true',
    default=False
)
parser.add_argument(
    '--db',
    help='Set up database',
    action='store_true',
    default=False
)
parser.add_argument(
    "-e",
    "--errorpath",
    default="errorlog.txt",
    help="Path to output errors"
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

if args.download:

    Pipeline.pipelines[args.pipeline.lower()].download_raw_images_for_night(
        night=args.night
    )

    logger.info("Download complete")

pipe = get_pipeline(
    args.pipeline,
    selected_configurations=args.config,
    night=args.night,
)

pipe.reduce_images([[[], []]], output_error_path=args.errorpath)

logger.info('End of winterdrp execution')
