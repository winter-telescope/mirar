"""
Script to fix headers of images based on a log file
"""
import argparse

import numpy as np
import pandas as pd
from astropy.io import fits

from mirar.paths import get_output_path


def fix_headers(logfile, night, keyword_list, sub_dir="raw"):
    """
    Function to fix headers of images based on a log file
    """
    obslog = pd.read_csv(logfile)
    for ind, image_basename in enumerate(obslog["BASENAME"]):
        for i in range(5):
            if f"_{i}.fits" in image_basename:
                image_basename = image_basename.replace(f"_{i}.fits", ".fits")
        image_filename = get_output_path(
            base_name=image_basename, dir_root=sub_dir, sub_dir="winter/" + night
        )
        img_hdulist = fits.open(image_filename, "update")
        img_header = img_hdulist[0].header
        for key in keyword_list:
            if obslog.loc[ind][key] is np.nan:
                img_header[key] = ""
            else:
                img_header[key] = obslog.loc[ind][key]

        img_hdulist.close()


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Fix headers of images")
    parser.add_argument("logfile", type=str, help="CSV file with image metadata")
    parser.add_argument("night", type=str, help="Night of observation")
    parser.add_argument(
        "--sub_dir", type=str, default="raw", help="Subdirectory of images"
    )
    parser.add_argument("-k", "--keywords", type=str, nargs="+")
    args = parser.parse_args()

    fix_headers(args.logfile, args.night, args.keywords, args.sub_dir)
