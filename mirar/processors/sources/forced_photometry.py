"""
Module to extract a candidates table from an image header
"""

import logging

import pandas as pd

from mirar.data import ImageBatch, SourceBatch, SourceTable
from mirar.data.utils import get_xy_from_wcs
from mirar.errors import ProcessorError
from mirar.paths import (
    CAND_DEC_KEY,
    CAND_RA_KEY,
    SOURCE_HISTORY_KEY,
    SOURCE_NAME_KEY,
    XPOS_KEY,
    YPOS_KEY,
)
from mirar.processors.base_processor import BaseSourceGenerator

logger = logging.getLogger(__name__)


class TableFromHeaderError(ProcessorError):
    """Error relating to writing a source table from a header"""


class HeaderKeyMissingError(ProcessorError):
    """Error relating to missing keys in headers"""


class ForcedPhotometryDetector(BaseSourceGenerator):
    """
    Class to create a candidates table for performing forced photometry
    """

    base_key = "forcedphot"

    def __init__(
        self,
        calculate_image_coordinates: bool = True,
        ra_header_key: str = CAND_RA_KEY,
        dec_header_key: str = CAND_DEC_KEY,
        name_header_key: str | None = None,
    ):
        super().__init__()
        self.calculate_image_coordinates = calculate_image_coordinates
        self.ra_header_key = ra_header_key
        self.dec_header_key = dec_header_key
        self.name_header_key = name_header_key

    def description(self) -> str:
        return (
            f"Perform forced photometry using metadata keys "
            f"'{self.ra_header_key}'&'{self.dec_header_key}'"
        )

    def _apply_to_images(self, batch: ImageBatch) -> SourceBatch:
        all_cands = SourceBatch()
        for image in batch:
            # Save the image here, to facilitate processing downstream
            header = image.get_header()

            metadata = self.get_metadata(image)

            name = None

            if self.name_header_key is not None:
                try:
                    name = header[self.name_header_key]
                except KeyError as exc:
                    err = (
                        f"Header key for name ({self.name_header_key}) not "
                        f"found in header. Available keys: {header.keys()}"
                    )
                    logger.error(err)
                    raise HeaderKeyMissingError(err) from exc

            new_dict = {
                CAND_RA_KEY: image[self.ra_header_key],
                CAND_DEC_KEY: image[self.dec_header_key],
                SOURCE_HISTORY_KEY: pd.DataFrame(),
                SOURCE_NAME_KEY: name,
            }

            if self.calculate_image_coordinates:
                x, y = get_xy_from_wcs(
                    ra_deg=image[self.ra_header_key],
                    dec_deg=image[self.dec_header_key],
                    header=header,
                )
                new_dict[XPOS_KEY] = x
                new_dict[YPOS_KEY] = y

            src_table = pd.DataFrame([new_dict])

            all_cands.append(SourceTable(src_table, metadata=metadata))

        return all_cands
