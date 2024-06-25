"""
Functions to filter and purify the source table and photometric catalogs
"""

import numpy as np
import pandas as pd
from astropy.table import Table

from mirar.data import SourceBatch


def wirc_source_table_filter_annotator(source_table: SourceBatch) -> SourceBatch:
    """
    Function to remove bad candidates with None in sigmapsf, magpsf, magap, sigmagap,
    and update the source table with the keys required for the WIRC database
    :param source_table: source table
    :return: updated source table
    """

    new_batch = SourceBatch([])

    for source in source_table:
        src_df = source.get_data()

        none_mask = (
            src_df.loc[:, "sigmapsf"].isnull()
            | src_df.loc[:, "magpsf"].isnull()
            | src_df.loc[:, "magap"].isnull()
            | src_df.loc[:, "sigmagap"].isnull()
        )

        mask = none_mask.values

        # Needing to do this because the dataframe is big-endian
        mask_inds = np.where(~mask)[0]
        src_df = pd.DataFrame([src_df.loc[x] for x in mask_inds]).reset_index(drop=True)

        source.set_data(src_df)
        new_batch.append(source)

    return new_batch


def wirc_zogy_catalogs_purifier(
    sci_catalog: Table, ref_catalog: Table
) -> (Table, Table):
    """
    Function to purify the photometric catalog
    :param sci_catalog:
    :param ref_catalog:
    :return: sci_catalog, ref_catalog
    """
    good_sci_sources = (
        (sci_catalog["FLAGS"] == 0)
        & (sci_catalog["SNR_WIN"] > 5)
        & (sci_catalog["FWHM_WORLD"] < 4.0 / 3600)
        & (sci_catalog["FWHM_WORLD"] > 0.5 / 3600)
        & (sci_catalog["SNR_WIN"] < 1000)
    )

    good_ref_sources = (
        (ref_catalog["FLAGS"] == 0)
        & (ref_catalog["SNR_WIN"] > 5)
        & (ref_catalog["FWHM_WORLD"] < 5.0 / 3600)
        & (ref_catalog["FWHM_WORLD"] > 0.5 / 3600)
        & (ref_catalog["SNR_WIN"] < 1000)
    )
    return good_sci_sources, good_ref_sources
