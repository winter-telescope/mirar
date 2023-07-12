"""
Module to write an observation log for a directory of images
"""
import argparse
from glob import glob

from astropy.io import fits
from astropy.table import Table

from mirar.paths import get_output_dir


def write_observation_log(image_dir):
    """
    Write an observation log for a directory of images
    """
    filelist = glob(image_dir.as_posix() + "/*.fits")
    keywords = [
        "FILTER",
        "EXPTIME",
        "OBSTYPE",
        "TARGNAME",
        "RADEG",
        "DECDEG",
        "BOARD_ID",
        "FIELDID",
    ]

    all_entries = []
    for filename in filelist:
        header = fits.getheader(filename)
        entries = [header[key] for key in keywords if key in header]
        found_keys = [key for key in keywords if key in header]
        all_entries.append(entries)

    tab = Table(rows=all_entries, names=found_keys)
    tab["filename"] = [filename.split("/")[-1] for filename in filelist]
    tab.write(f"{directory}/obslog.csv", format="ascii.csv", overwrite=True)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Fix headers of images")
    parser.add_argument("night", type=str, help="Night of observation")
    parser.add_argument(
        "--sub_dir", type=str, default="raw", help="Subdirectory of images"
    )
    args = parser.parse_args()

    directory = get_output_dir(args.sub_dir, f"winter/{args.night}")
    write_observation_log(directory)
