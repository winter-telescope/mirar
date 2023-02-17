"""
Script containing functions to run anet.net locally
"""
import logging
import os
from typing import Optional

from astropy.io import fits

from winterdrp.utils import ExecutionError, execute

logger = logging.getLogger(__name__)


class AstrometryNetError(ExecutionError):
    """Astrometry.net-related error"""


def run_astrometry_net(images: str | list, output_dir: str, *args, **kwargs):
    """
    function to execute `run_astrometry_net_single` on several images in batch
    """
    if not isinstance(images, list):
        images = [images]

    # make output directory if it doesn't exist

    try:
        os.makedirs(output_dir)
    except OSError:
        pass

    for img in images:
        run_astrometry_net_single(img, output_dir, *args, **kwargs)


def run_astrometry_net_single(
    img: str,
    output_dir: str,
    scale_bounds: Optional[tuple | list] = None,  # limits on scale (lower, upper)
    scale_units: Optional[str] = None,  # scale units ('degw', 'amw')
    downsample: Optional[float | int] = None,  # downsample by factor of __
):
    """
    function to run anet.net locally on one image, with options to adjust settings
    default: solve-field <img> -D <output_dir> -N <newname> -O
    """
    # name for new file if a-net solves (otherwise a-net writes to '<img>.new')
    newname = str(img).split("a-net/")[1]  # should generalize this

    # run a-net (solve-field)
    cmd = (
        f"solve-field {img} "
        f"-D {output_dir} "
        f"-N {newname} "  # instead of new-image
        f"-O "  # overwrite
        # f"-N {img} "  # instead of new-image
    )

    if scale_bounds is not None:
        cmd += f" --scale-high {max(scale_bounds)} "
        cmd += f" --scale-low {min(scale_bounds)} "

    if scale_units is not None:
        cmd += f"--scale-units {scale_units} "

    if downsample is not None:
        cmd += f"-z {downsample} "

    # cmd with a ra, dec first guess (speeds up solution)
    header = fits.open(img)[0].header  # pylint: disable=no-member
    ra_req, dec_req = header["RA"], header["DEC"]  # requested ra, dec
    cmd_loc = (
        cmd + f"--ra {ra_req}, --dec {dec_req} --radius 5"
    )  # radius takes on units of ra, dec

    try:
        execute(cmd_loc, output_dir)
        if os.path.isfile(img):
            print("solve file created,")
            print("ran with ra,dec guess \n command: ", cmd_loc)
        else:
            execute(cmd, output_dir)
            if os.path.isfile(img):
                print("solve file created,")
                print("ran without ra,dec guess \n command: ", cmd)
    except ExecutionError as err:
        raise AstrometryNetError(err) from err

    # remove 'HISTORY' keywords (breaks the pipeline)
    solved = fits.open(img)
    solved_hdr = solved[0].header  # pylint: disable=no-member
    del solved_hdr["HISTORY"]
    # save to same file, overwrite
    fits.writeto(  # pylint: disable=no-member
        img, solved[0].data, solved_hdr, overwrite=True  # pylint: disable=no-member
    )  # pylint: disable=no-member

    return output_dir
