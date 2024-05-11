"""
Module containing WINTER CLI commands
"""

import argparse
import logging
import sys
import tempfile
from pathlib import Path

import numpy as np
from astropy.time import Time
from sqlalchemy import Select, text
from sqlalchemy.sql import func

from mirar.data import Dataset, ImageBatch, cache
from mirar.database.transactions.select import run_select
from mirar.io import open_raw_image
from mirar.paths import TEMP_DIR
from mirar.pipelines.winter.models import ExposuresTable, RawsTable, StacksTable
from mirar.pipelines.winter.winter_pipeline import WINTERPipeline
from mirar.processors.utils.image_loader import load_from_list

logger = logging.getLogger(__name__)


def run_winter(
    config: str,
    dataset: Dataset | None = None,
    log_level: str = "INFO",
    subdir: str = None,
):
    """
    Run a WINTER pipeline

    :param config:
    :param dataset:
    :param log_level:
    :param subdir:
    :return:
    """

    # Set up logging

    log = logging.getLogger("mirar")

    handler = logging.StreamHandler(sys.stdout)

    formatter = logging.Formatter(
        "%(name)s [l %(lineno)d] - %(levelname)s - %(message)s"
    )
    handler.setFormatter(formatter)
    log.addHandler(handler)
    log.setLevel(log_level)

    pipe = WINTERPipeline(
        selected_configurations=config,
        night=subdir,
    )

    _, errorstack = pipe.reduce_images(
        dataset=dataset,
        catch_all_errors=True,
    )

    errorstack.summarise_error_stack()


def run_stack_of_stacks():
    """
    CLI wrapper for stacking pre-existing stack images of a target

    :return:
    """
    parser = argparse.ArgumentParser(description="Stack of stacks")
    parser.add_argument("-t", "--target", help="Target to stack", type=str)
    parser.add_argument("-f", "--fieldid", help="Field ID", type=str)
    parser.add_argument("--level", help="Logging level", type=str, default="INFO")
    parser.add_argument(
        "--bindays", help="Window, in days, for bin", type=int, default=None
    )
    args = parser.parse_args()

    if args.target is not None:
        db_constraint = f"targname = '{args.target}'"
        subdir = args.target
        if args.fieldid is not None:
            err = "Cannot specify both target and fieldid"
            logger.error(err)
            raise ValueError(err)
    elif args.fieldid is not None:
        db_constraint = f"fieldid = '{args.fieldid}'"
        subdir = args.fieldid
    else:
        err = "Must specify either target or fieldid"
        logger.error(err)
        raise ValueError(err)

    sel = (
        Select(
            StacksTable.savepath,
            func.min(ExposuresTable.utctime).label("utctime"),
        )
        .join(StacksTable.raw)
        .join(RawsTable.exposure_ids)
        .group_by(StacksTable.savepath)
        .where(text(db_constraint))
    )

    df = run_select(sel, StacksTable)
    df.sort_values(by="utctime", inplace=True)
    df.reset_index(drop=True, inplace=True)

    with tempfile.TemporaryDirectory(dir=TEMP_DIR) as temp_dir_path:
        print(f"Using cache {temp_dir_path}")

        cache.set_cache_dir(temp_dir_path)

        savepaths = []

        for x in df["savepath"]:
            path = Path(x)
            if not path.exists():
                err = f"Stack file {path} does not exist"
                logger.warning(err)

            savepaths.append(path)

        print(f"Found {len(savepaths)} stack images")

        savepaths = [Path(x) for x in df["savepath"]]

        img_batch = load_from_list(savepaths, open_raw_image)

        if args.bindays is not None:

            times = np.array(
                [
                    Time(str(x["utctime"].to_datetime64()), format="isot").mjd
                    for _, x in df.iterrows()
                ]
            )

            bins = np.arange(min(times), max(times) + args.bindays, step=args.bindays)

            new_batch = ImageBatch()

            images = np.array(img_batch.get_batch())

            for i, t_start in enumerate(bins[:-1]):
                t_end = bins[i + 1]

                mask = (times >= t_start) & (times < t_end)

                for img in images[mask]:
                    targ_name = f"{subdir}_{i}"
                    img["TARGNAME"] = targ_name
                    img["RESTACKN"] = i
                    img["BINSTART"] = t_start
                    img["BINEND"] = t_end

                    new_batch.append(img)

            img_batch = new_batch

        dataset = Dataset(img_batch)

        run_winter(
            config="stack_stacks_db",
            dataset=dataset,
            log_level=args.level,
            subdir=f"restacks/{subdir}",
        )
