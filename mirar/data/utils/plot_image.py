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
    regions_wcs_coords: List[Tuple[float, float]] | None = None,
    plot_format: str = "png",
    title_fields: List[str] | None = None,
    fig: plt.Figure | None = None,
):
    """
    Plot the fits image with the specified regions
    Args:
        :param image: Image to plot
        :param savedir: Directory to save to.
        :param regions_wcs_coords:If you want to mark specific coordinates on the image,
        provide a list of tuples of RA, Dec
        :param plot_format: pdf or png
        :param title_fields: Image header fields to annotate the plot with
        :param fig: If you want to plot on an existing figure, provide it here
    """
    assert plot_format in ["pdf", "png"], (
        f"Only pdf and png formats are supported, " f"got {plot_format}."
    )
    if not isinstance(savedir, Path):
        savedir = Path(savedir)

    if fig is None:
        fig = plt.figure()
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

    if title_fields is not None:
        assert all(
            [field in image.header for field in title_fields]
        ), f"Title fields {title_fields} not present in header {image.header}"
        # Make a string for title in the format key:value
        hdr = image.get_header()
        title_str = ""
        for field in title_fields:
            title_str += field + ": "
            val = hdr[field]
            if isinstance(val, float):
                val = f"{val:.2f}"
            title_str += str(val) + ",   "

        # Wrap title if it's too long
        title = "\n".join([title_str[i : i + 80] for i in range(0, len(title_str), 80)])
        # Set title that has latex
        ax.set_title(rf"{title}", fontsize=10)

    ax.set_xlabel("RA")
    ax.set_ylabel("Dec")

    plot_savepath = savedir / image.header[BASE_NAME_KEY].replace(
        ".fits", f".{plot_format}"
    )
    logger.debug(f"Saving plot to {plot_savepath}")
    fig.savefig(plot_savepath, bbox_inches="tight")
    plt.close(fig)
