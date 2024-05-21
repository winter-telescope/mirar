"""
Module to plot 2D image-data as a pdf or png file.
"""

import logging

import matplotlib
import matplotlib.pyplot as plt

from mirar.data import ImageBatch
from mirar.data.utils.plot_image import plot_fits_image
from mirar.paths import get_output_dir
from mirar.processors.base_processor import BaseImageProcessor

matplotlib.use("Agg")

logger = logging.getLogger(__name__)


class ImagePlotter(BaseImageProcessor):
    """Processor to plot images.
    Attributes
    :param plot_format: pdf or png?
    """

    base_key = "plot"

    def __init__(
        self,
        output_sub_dir: str = "plots",
        plot_format: str = "png",
        annotate_fields: str | list[str] | None = None,
    ):
        super().__init__()
        self.output_sub_dir = output_sub_dir
        assert plot_format in ["pdf", "png"], (
            f"Only pdf and png formats are " f"supported, got {plot_format}."
        )
        self.plot_format = plot_format
        if isinstance(annotate_fields, str):
            annotate_fields = [annotate_fields]
        self.annotate_fields = annotate_fields

    def description(self):
        return (
            f"Processor to plot images as {self.plot_format} and save them "
            f"in the '{self.output_sub_dir}' sub-directory"
        )

    def _apply_to_images(self, batch: ImageBatch) -> ImageBatch:
        output_dir = get_output_dir(
            dir_root=self.output_sub_dir, sub_dir=self.night_sub_dir
        )
        output_dir.mkdir(parents=True, exist_ok=True)
        for image in batch:
            # We use multithreading to plot the images, so we need to make sure
            # the axes are not shared between threads
            fig = plt.figure()
            plot_fits_image(
                image=image,
                savedir=output_dir,
                title_fields=self.annotate_fields,
                plot_format=self.plot_format,
                fig=fig,
            )

        return batch
