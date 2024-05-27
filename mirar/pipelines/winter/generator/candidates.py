"""
Module for candidate-related functions for WINTER
"""

import numpy as np
import pandas as pd
from astropy.time import Time

from mirar.data import SourceBatch
from mirar.data.utils.compress import decode_img
from mirar.errors import ProcessorError
from mirar.paths import MAGLIM_KEY, SOURCE_HISTORY_KEY, TIME_KEY, ZP_KEY, ZP_STD_KEY
from mirar.pipelines.winter.constants import sncosmo_filters, winter_inv_filters_map
from mirar.processors.skyportal.skyportal_source import SNCOSMO_KEY


class NoGoodCandidatesError(ProcessorError):
    """Error raised when no candidates pass quality cuts"""


def winter_candidate_annotator_filterer(source_batch: SourceBatch) -> SourceBatch:
    """
    Function to perform basic filtering to weed out bad candidates with None
    magnitudes, to be added.
    :param source_batch: Source batch
    :return: updated batch
    """

    new_batch = []

    for source in source_batch:
        src_df = source.get_data()

        bad_sources_mask = (
            src_df["sigmapsf"].isnull()
            | src_df["magpsf"].isnull()
            | src_df["magap"].isnull()
            | src_df["sigmagap"].isnull()
            | (src_df["fwhm"] <= 0)
            | (src_df["scorr"] < 0)
        )

        mask = bad_sources_mask.values

        # Needing to do this because the dataframe is big-endian
        mask_inds = np.where(~mask)[0]
        filtered_df = pd.DataFrame([src_df.loc[x] for x in mask_inds]).reset_index(
            drop=True
        )

        if len(filtered_df) == 0:
            filtered_df = pd.DataFrame(columns=src_df.columns)

        # Pipeline (db) specific keywords
        source["magzpsci"] = source[ZP_KEY]
        source["magzpsciunc"] = source[ZP_STD_KEY]
        source["diffmaglim"] = source[MAGLIM_KEY]
        source["programpi"] = source["PROGPI"]
        source["programid"] = source["PROGID"]
        source["field"] = source["FIELDID"]
        source["jd"] = Time(source[TIME_KEY]).jd
        source["boardid"] = source["BOARD_ID"]

        source.set_data(filtered_df)
        if len(filtered_df) > 0:
            new_batch.append(source)

    if len(new_batch) == 0:
        raise NoGoodCandidatesError("No candidates passed quality filter")

    return SourceBatch(new_batch)


def winter_new_source_updater(source_table: SourceBatch) -> SourceBatch:
    """
    Function to add relevant fields for new sources

    :param source_table: Original source table
    :return: Updated source table
    """
    for source in source_table:
        src_df = source.get_data()

        src_df["ndet"] = 1
        src_df["average_ra"] = src_df["ra"]
        src_df["average_dec"] = src_df["dec"]
        src_df["first_det_utc"] = source[TIME_KEY]
        src_df["latest_det_utc"] = source[TIME_KEY]

        source.set_data(src_df)

    return source_table


def winter_source_entry_updater(source_table: SourceBatch) -> SourceBatch:
    """
    Function to update the source table with new source averages

    :param source_table: Original source table
    :return: Updated source table
    """
    for source in source_table:
        src_df = source.get_data()

        hist_dfs = [
            pd.DataFrame(src_df[SOURCE_HISTORY_KEY].loc[x]) for x in range(len(src_df))
        ]

        src_df["ndet"] = [len(x) for x in hist_dfs]

        new_fields = []

        # FIXME remove same detection by candid

        for i, hist_df in enumerate(hist_dfs):
            if len(hist_df) == 1:

                new_hist_df = pd.DataFrame(columns=hist_df.columns)

                new_fields.append(
                    {
                        "average_ra": src_df["ra"].iloc[i],
                        "average_dec": src_df["dec"].iloc[i],
                        "first_det_utc": source[TIME_KEY],
                        "latest_det_utc": source[TIME_KEY],
                        "jdstarthist": source["jd"],
                        "jdendhist": source["jd"],
                        "ndethist": 0,
                        SOURCE_HISTORY_KEY: new_hist_df,
                    }
                )

            else:

                ras = np.array(hist_df["ra"].tolist() + [src_df["ra"].iloc[i]])
                decs = np.array(hist_df["dec"].tolist() + [src_df["dec"].iloc[i]])
                weights = 1.0 / np.array(
                    hist_df["sigmapsf"].tolist() + [src_df["sigmapsf"].iloc[i]]
                )

                # Wrap around the RA if split at 0
                if np.max(ras) > 350.0:
                    ras[ras < 180.0] += 360.0

                av_ra = np.average(ras, weights=weights)
                av_dec = np.average(decs, weights=weights)

                # Unwrap the RA
                if av_ra > 360.0:
                    av_ra -= 360.0

                min_jd = min(hist_df["jd"].tolist() + [source["jd"]])
                max_jd = max(hist_df["jd"].tolist() + [source["jd"]])

                new_hist_df = hist_df[hist_df["candid"] != src_df["candid"].iloc[i]]

                new_fields.append(
                    {
                        "average_ra": av_ra,
                        "average_dec": av_dec,
                        "first_det_utc": Time(min_jd, format="jd").isot,
                        "latest_det_utc": Time(max_jd, format="jd").isot,
                        "jdstarthist": min_jd,
                        "jdendhist": max_jd,
                        "ndethist": len(new_hist_df),
                        SOURCE_HISTORY_KEY: new_hist_df,
                    }
                )

        new = pd.DataFrame(new_fields)

        for column in new.columns:
            src_df[column] = new[column]

        source.set_data(src_df)

    return source_table


def winter_candidate_avro_fields_calculator(source_table: SourceBatch) -> SourceBatch:
    """
    Function to calculate the AVRO fields for WINTER
    """

    new_batch = SourceBatch([])

    for source in source_table:
        src_df = source.get_data()

        src_df["magdiff"] = src_df["magpsf"] - src_df["magap"]
        src_df["magfromlim"] = source["diffmaglim"] - src_df["magpsf"]

        src_df["utctime"] = source[TIME_KEY]

        src_df["d_to_x"] = src_df["NAXIS1"] - src_df["xpos"]
        src_df["d_to_y"] = src_df["NAXIS2"] - src_df["ypos"]
        src_df["mindtoedge"] = src_df[["xpos", "ypos", "d_to_x", "d_to_y"]].min(axis=1)
        nnegs, nbads, sumrat, frac_masked = [], [], [], []
        for _, src in src_df.iterrows():
            diff_cutout_data = decode_img(src["cutout_difference"])
            # Get central 5x5 pixels
            nx, ny = diff_cutout_data.shape
            diff_stamp = diff_cutout_data[
                nx // 2 - 3 : nx // 2 + 2, ny // 2 - 3 : ny // 2 + 2
            ]
            nnegs.append(np.sum(diff_stamp < 0))
            nbads.append(np.sum(np.isnan(diff_stamp)))
            sumrat.append(np.sum(diff_stamp) / np.sum(np.abs(diff_stamp)))
            frac_masked.append(
                np.sum(np.isnan(diff_cutout_data))
                / (diff_cutout_data.shape[0] * diff_cutout_data.shape[1])
            )

        src_df["nneg"] = nnegs
        src_df["nbad"] = nbads
        src_df["sumrat"] = sumrat
        src_df["fracmasked"] = frac_masked

        for column in ["gaia_parallax_over_error1", "gaiabright_parallax_over_error1"]:
            mask = ~pd.isnull(src_df[column])
            src_df.loc[mask, [column]] = abs(src_df.loc[mask, [column]])

        src_df["distgaia"] = src_df["distgaianr1"]
        src_df["plxgaia"] = src_df["gaia_parallax_over_error1"]
        src_df["ruwegaia"] = src_df["gaia_ruwe1"]

        src_df["distgaiabright"] = src_df["distgaiabrightnr1"]
        src_df["plxgaiabright"] = src_df["gaiabright_parallax_over_error1"]
        src_df["ruwegaiabright"] = src_df["gaiabright_ruwe1"]

        source.set_data(src_df)
        new_batch.append(source)

    return new_batch


def winter_skyportal_annotator(source_batch: SourceBatch) -> SourceBatch:
    """
    Function to update the candidate table with the skyportal fields

    :param source_batch: Original source table
    :return: Updated source table
    """

    for source_table in source_batch:
        src_df = source_table.get_data()

        src_df["ndethist"] = [len(x) for x in src_df[SOURCE_HISTORY_KEY]]

        if "fid" not in src_df.columns:
            src_df["fid"] = source_table["FID"]

        if SNCOSMO_KEY not in src_df.columns:
            sncosmo_fs = [
                sncosmo_filters[winter_inv_filters_map[x].lower()[0]]
                for x in src_df["fid"]
            ]
            src_df[SNCOSMO_KEY] = sncosmo_fs

        for hist_df in src_df[SOURCE_HISTORY_KEY]:
            mjds = [Time(x, format="jd").mjd for x in hist_df["jd"]]
            hist_df["mjd"] = mjds
            sncosmo_fs = [
                sncosmo_filters[winter_inv_filters_map[x].lower()[0]]
                for x in hist_df["fid"]
            ]
            hist_df[SNCOSMO_KEY] = sncosmo_fs

        source_table.set_data(src_df)

    return source_batch


def winter_candidate_quality_filterer(source_table: SourceBatch) -> SourceBatch:
    """
    Function to perform quality filtering on WINTER candidates
    """
    new_batch = []

    min_dist_to_star = 7.0

    for source in source_table:
        src_df = source.get_data()

        mask = (
            ((src_df["rb"] > 0.1) | pd.isnull(src_df["rb"]))
            & (src_df["fwhm"] < 10.0)
            & (src_df["mindtoedge"] > 50.0)
            & (src_df["isdiffpos"])
            & (  # Cut on sgscore1
                (src_df["sgscore1"] < 0.5)
                | pd.isnull(src_df["sgscore1"])
                | (src_df["distpsnr1"] > min_dist_to_star)
            )
            & (  # Cut on PS1STRM Star Probability
                (src_df["ps1strmprobstar1"] < 0.5)
                | pd.isnull(src_df["ps1strmprobstar1"])
                | (src_df["distpsnr1"] > min_dist_to_star)
            )
            & (src_df["ndethist"] > 0)
        )
        filtered_df = src_df[mask].reset_index(drop=True)

        if len(filtered_df) > 0:
            source.set_data(filtered_df)
            new_batch.append(source)

    if len(new_batch) == 0:
        raise NoGoodCandidatesError("No candidates passed quality filter")

    return SourceBatch(new_batch)
