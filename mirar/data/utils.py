"""
Utilities for working with data
"""
import logging
from pathlib import Path
from typing import List, Tuple

import matplotlib
import matplotlib.pyplot as plt
from astropy.stats import sigma_clipped_stats
from astropy.wcs import WCS

from mirar.data.image_data import Image
from mirar.paths import BASE_NAME_KEY

matplotlib.use("Agg")
logger = logging.getLogger(__name__)


def plot_fits_image(
    image: Image,
    savedir: str | Path,
    regions_wcs_coords: List[Tuple[float, float]] = None,
):
    """
    Plot the fits image with the regions
    Args:
        image: Image to plot
        savedir: Directory to save to
        regions_wcs_coords:If you want to mark specific coordinates on the image,
        provide a list of tuples of RA, Dec

    Returns:

    """
    if not isinstance(savedir, Path):
        savedir = Path(savedir)

    fig = plt.figure(figsize=(10, 10))
    ax = fig.add_subplot(111)
    data = image.get_data()
    _, median, std = sigma_clipped_stats(data)
    ax.imshow(
        data, origin="lower", cmap="gray", vmin=median - 5 * std, vmax=median + 5 * std
    )
    if regions_wcs_coords is not None:
        regions_image_coords = WCS(image.header).all_world2pix(regions_wcs_coords, 1)
        for _, (ra, dec) in enumerate(regions_image_coords):
            ax.scatter(ra, dec, marker="x", color="red", s=100)

    ax.set_xlabel("RA")
    ax.set_ylabel("Dec")
    plot_savepath = savedir / image.header[BASE_NAME_KEY].replace(".fits", ".png")
    logger.info(f"Saving plot to {plot_savepath}")
    plt.savefig(plot_savepath, bbox_inches="tight")
    plt.close()
