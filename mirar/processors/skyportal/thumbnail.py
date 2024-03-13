"""
Function for converting a FITS image to a skyportal "thumbnail".
"""

import base64
import io

import matplotlib
import matplotlib.pyplot as plt
import numpy as np
from astropy.visualization import (
    AsymmetricPercentileInterval,
    ImageNormalize,
    LinearStretch,
    LogStretch,
)

matplotlib.use("agg")


def make_thumbnail(
    image_data: np.ndarray,
    linear_stretch: bool = False,
) -> str:
    """
    Util function to convert a FITS image to a PNG image

    Skyportal requires, to quote the API:

        base64-encoded PNG image file contents.
        Image size must be between 16px and 500px on a side.

    :param image_data: Image data
    :param linear_stretch: boolean whether to use a linear stretch (default is log)
    :return: Skyportal-compliant PNG image string
    """
    with io.BytesIO() as buff:
        fig = plt.figure()
        fig.set_size_inches(4, 4, forward=False)
        ax_1 = plt.Axes(fig, [0.0, 0.0, 1.0, 1.0])
        ax_1.set_axis_off()
        fig.add_axes(ax_1)

        # replace nans with median:
        img = np.array(image_data)
        # replace dubiously large values
        xl_mask = np.greater(np.abs(img), 1e20, where=~np.isnan(img))
        if img[xl_mask].any():
            img[xl_mask] = np.nan
        if np.isnan(img).any():
            median = float(np.nanmean(img.flatten()))
            img = np.nan_to_num(img, nan=median)

        norm = ImageNormalize(
            img,
            stretch=LinearStretch() if linear_stretch else LogStretch(),
        )
        img_norm = norm(img)
        normalizer = AsymmetricPercentileInterval(
            lower_percentile=1, upper_percentile=100
        )
        vmin, vmax = normalizer.get_limits(img_norm)
        ax_1.imshow(img_norm, cmap="bone", origin="lower", vmin=vmin, vmax=vmax)
        plt.savefig(buff, dpi=42)

        buff.seek(0)
        plt.close()

        fritz_thumbnail = base64.b64encode(buff.read()).decode("utf-8")

    return fritz_thumbnail
