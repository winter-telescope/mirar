"""
Module to extract a candidates table from an image header
"""
from typing import Optional

import numpy as np
import pandas as pd
from astropy.io.fits import Header

from mirar.data import ImageBatch, SourceBatch, SourceTable
from mirar.errors import ProcessorError
from mirar.paths import (
    BASE_NAME_KEY,
    CAND_DEC_KEY,
    CAND_NAME_KEY,
    CAND_RA_KEY,
    LATEST_SAVE_KEY,
    NORM_PSFEX_KEY,
    UNC_IMG_KEY,
    XPOS_KEY,
    YPOS_KEY,
    ZP_KEY,
    base_output_dir,
    core_fields,
    get_output_path,
)
from mirar.processors.base_processor import BaseCandidateGenerator
from mirar.processors.candidates.utils import get_xy_from_wcs


class TableFromHeaderError(ProcessorError):
    """Error relating to writing a source table from a header"""


class HeaderKeyMissingError(ProcessorError):
    """Error relating to missing keys in headers"""


class SourceTablefromHeader(BaseCandidateGenerator):
    """
    Class to create a source table from an image header with user specified header keys
    """

    base_key = "header_source_exporter"

    def __init__(
        self,
        header_keys_list: list[str],
        calculate_image_coordinates: bool = False,
        ra_header_key: str = CAND_RA_KEY,
        dec_header_key: str = CAND_DEC_KEY,
        save: bool = True,
        output_sub_dir: Optional[str] = "sourcetable",
        output_dir: str = base_output_dir,
        column_names_list: list[str] = None,
    ):
        super().__init__()
        self.header_keys_list = np.array(header_keys_list)
        self.calculate_image_coordinates = calculate_image_coordinates
        self.ra_header_key = ra_header_key
        self.dec_header_key = dec_header_key
        self.save = save
        self.output_sub_dir = output_sub_dir
        self.output_dir = output_dir
        self.column_names_list = column_names_list
        if self.column_names_list is None:
            self.column_names_list = np.array(header_keys_list)

        assert len(self.column_names_list) == len(self.header_keys_list)
        if np.logical_and(
            self.calculate_image_coordinates,
            np.logical_or(
                CAND_RA_KEY not in self.column_names_list,
                CAND_DEC_KEY not in self.column_names_list,
            ),
        ):
            err = (
                f"The default keys for ra/dec : {CAND_RA_KEY}, {CAND_DEC_KEY} are "
                f"not present in the header_keys_list. Please provide ra_header_key"
                f" and dec_header_key, or disable x/y coordinate calculation."
            )
            raise TableFromHeaderError(err)

    def append_keys_to_header(self, header: Header):
        """
        Append additional keys to the header
        Args:
            header: header to append keys to

        Returns:
            header: Header with appended keys
        """
        return header

    def _apply_to_images(self, batch: ImageBatch) -> SourceBatch:
        all_cands = SourceBatch()
        for image in batch:
            # Save the image here, to facilitate processing downstream
            header = image.get_header()
            header = self.append_keys_to_header(header)
            image.set_header(header)

            if self.save:
                path = get_output_path(
                    base_name=image[BASE_NAME_KEY],
                    dir_root=self.output_sub_dir,
                    sub_dir=self.night_sub_dir,
                )
                if not path.parent.exists():
                    path.parent.mkdir(parents=True)
                self.save_fits(image, path)

            key_absent = [x not in header.keys() for x in self.header_keys_list]

            if np.sum(key_absent) > 0:
                err = f"Keys {self.header_keys_list[key_absent]} absent from header"
                raise HeaderKeyMissingError(err)

            cands_table = pd.DataFrame()

            for ind, key in enumerate(self.header_keys_list):
                val = header[key]
                if not isinstance(val, list):
                    val = [val]
                cands_table[self.column_names_list[ind]] = val

            if self.calculate_image_coordinates:
                x, y = get_xy_from_wcs(
                    ra_deg=list(cands_table[self.ra_header_key]),
                    dec_deg=list(cands_table[self.dec_header_key]),
                    header=header,
                )
                cands_table[XPOS_KEY] = x
                cands_table[YPOS_KEY] = y

            metadata = {}
            for key in core_fields:
                metadata[key] = image[key]

            all_cands.append(SourceTable(cands_table, metadata=metadata))

        return all_cands


class ForcedPhotometryCandidateTable(SourceTablefromHeader):
    """
    Class to create a candidates table for performing forced photometry
    """

    base_key = "forcedphot_table_exporter"

    def __init__(self, ra_header_key: str, dec_header_key: str, name_header_key: str):
        header_keys_list = [
            ra_header_key,
            dec_header_key,
            name_header_key,
            LATEST_SAVE_KEY,
            UNC_IMG_KEY,
            NORM_PSFEX_KEY,
            ZP_KEY,
        ]
        column_names_list = [
            CAND_RA_KEY,
            CAND_DEC_KEY,
            CAND_NAME_KEY,
            LATEST_SAVE_KEY,
            UNC_IMG_KEY,
            NORM_PSFEX_KEY,
            ZP_KEY,
        ]
        super().__init__(
            header_keys_list=header_keys_list,
            calculate_image_coordinates=True,
            column_names_list=column_names_list,
        )

    def append_keys_to_header(self, header: Header):
        if UNC_IMG_KEY not in header.keys():
            header[UNC_IMG_KEY] = None
        if NORM_PSFEX_KEY not in header.keys():
            header[NORM_PSFEX_KEY] = None
        return header


class SourceTablefromCoordinates(ForcedPhotometryCandidateTable):
    """
    Class to create a source table to run forced photometry from user-specified
    coordinates
    """

    base_key = "coords_table_exporter"

    def __init__(
        self,
        ra_deg: float | list[float],
        dec_deg: float | list[float],
        ra_header_key: str = "USERRA",
        dec_header_key: str = "USERDEC",
        name_header_key: str = "USERNAME",
    ):
        self.ra_deg = ra_deg
        self.dec_deg = dec_deg
        super().__init__(
            ra_header_key=ra_header_key,
            dec_header_key=dec_header_key,
            name_header_key=name_header_key,
        )

    def append_keys_to_header(self, header: Header):
        header["USERRA"] = self.ra_deg
        header["USERDEC"] = self.dec_deg
        if UNC_IMG_KEY not in header.keys():
            header[UNC_IMG_KEY] = None
        if NORM_PSFEX_KEY not in header.keys():
            header[NORM_PSFEX_KEY] = None
        return header
