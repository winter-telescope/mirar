"""
Ingest reference (UKIRT/VIRTA) images into the database
"""

import argparse
import logging
import sys
from glob import glob
from pathlib import Path

from astropy.io import fits

from mirar.database.user import PostgresUser
from mirar.paths import LATEST_SAVE_KEY, base_output_dir, get_output_dir
from mirar.pipelines.winter.models import RefComponent

logger = logging.getLogger(__name__)


def get_logger(level="INFO"):
    """
    Get a logger for the module
    """
    log = logging.getLogger("mirar")
    handler = logging.StreamHandler(sys.stdout)
    formatter = logging.Formatter(
        "%(asctime)s %(name)s [l %(lineno)d] - %(levelname)s - %(message)s"
    )
    handler.setFormatter(formatter)
    log.addHandler(handler)
    log.setLevel(level)
    return log


def export_image_to_db(path: str, db_table=RefComponent, pg_user=PostgresUser()):
    """
    Export a fits image to the database
    """
    with fits.open(path, "update") as hdul:
        header = hdul[0].header  # pylint: disable=no-member
        header[LATEST_SAVE_KEY] = path
        sequence_key_names, sequence_values = pg_user.export_to_db(
            db_table=db_table, value_dict=header, duplicate_protocol="ignore"
        )
        for idx, key in enumerate(sequence_key_names):
            header[key] = sequence_values[idx]


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("-sub_dir", type=str, default="components")
    parser.add_argument("-base_output_dir", type=str, default=base_output_dir)

    parser.add_argument("-base_ref_dir", type=str, default=base_output_dir)
    parser.add_argument("-level", type=str, default="INFO")
    args = parser.parse_args()

    logger = get_logger(args.level)
    output_dir = get_output_dir(
        dir_root=args.sub_dir,
        sub_dir="winter/references",
        output_dir=Path(args.base_output_dir),
    )

    ref_imglist = glob(f"{output_dir}/*.fits")

    logger.info(f"Found {len(ref_imglist)} reference images in {output_dir}")
    for ind, ref_img in enumerate(ref_imglist):
        export_image_to_db(ref_img)
        logger.info(f"Exported {ref_img} to database : {ind+1}/{len(ref_imglist)}")
