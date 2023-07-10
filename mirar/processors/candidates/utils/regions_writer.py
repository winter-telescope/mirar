"""
Module with classes to write a regions file from a pandas dataframe
"""
import logging
import os
from pathlib import Path
from typing import Optional

from mirar.data import SourceBatch
from mirar.paths import base_output_dir, get_output_path
from mirar.processors.base_processor import BaseDataframeProcessor

logger = logging.getLogger(__name__)


def write_regions_file(
    regions_path, x_coords, y_coords, system="image", region_radius=5
):
    """
    Function to write a regions file
    Args:
        regions_path: str, path to regions file
        x_coords: list, x-coordinates or RA
        y_coords: list, y-coordinates or Dec
        system: str, image or wcs
        region_radius: float, radius of circle

    Returns:

    """
    logger.debug(f"Writing regions path to {regions_path}")
    with open(f"{regions_path}", "w", encoding="utf8") as regions_f:
        regions_f.write(f"{system}\n")
        for ind, x in enumerate(x_coords):
            regions_f.write(f"CIRCLE({x},{y_coords[ind]},{region_radius})\n")


class RegionsWriter(BaseDataframeProcessor):
    """
    Class to write a regions file from candidate table
    """

    base_key = "REGWRIT"

    def __init__(
        self,
        output_dir_name: Optional[str] = None,
        region_pix_radius: float = 8,
        output_dir: str | Path = base_output_dir,
        *args,
        **kwargs,
    ):
        super().__init__(*args, **kwargs)
        self.output_dir_name = output_dir_name
        self.region_pix_radius = region_pix_radius
        self.output_dir = Path(output_dir)

    def __str__(self) -> str:
        return (
            f"Processor to save 'region files' to "
            f"the {self.output_dir_name} directory. "
        )

    def _apply_to_candidates(
        self,
        batch: SourceBatch,
    ) -> SourceBatch:
        for source_list in batch:
            candidate_table = source_list.get_data()

            started_regions_paths = []
            for ind in range(len(candidate_table)):
                row = candidate_table.iloc[ind]
                regions_basepath = os.path.basename(row["diffimname"]).replace(
                    ".fits", ".reg"
                )
                regions_path = get_output_path(
                    regions_basepath,
                    dir_root=self.output_dir_name,
                    sub_dir=self.night_sub_dir,
                    output_dir=self.output_dir,
                )

                regions_path.parent.mkdir(parents=True, exist_ok=True)

                if regions_path not in started_regions_paths:
                    logger.debug(f"Writing regions path to {regions_path}")
                    with open(f"{regions_path}", "w", encoding="utf8") as regions_f:
                        regions_f.write("image\n")
                    started_regions_paths.append(regions_path)

                with open(f"{regions_path}", "w", encoding="utf8") as regions_f:
                    regions_f.write(
                        f"CIRCLE({row['X_IMAGE']},{row['Y_IMAGE']},"
                        f"{self.region_pix_radius})\n"
                    )

        return batch
