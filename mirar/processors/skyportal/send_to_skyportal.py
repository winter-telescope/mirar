"""
Module for sending candidates to Fritz.
"""
import base64
import gzip
import io
import logging
import time
from copy import deepcopy
from datetime import datetime
from typing import Mapping, Optional

import matplotlib
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import requests
from astropy.io import fits
from astropy.stats import sigma_clipped_stats
from astropy.time import Time
from astropy.visualization import (
    AsymmetricPercentileInterval,
    ImageNormalize,
    LinearStretch,
    LogStretch,
)

from mirar.data import SourceBatch, SourceTable
from mirar.paths import CAND_NAME_KEY, PACKAGE_NAME, SOURCE_HISTORY_KEY, __version__
from mirar.processors.base_processor import BaseSourceProcessor
from mirar.processors.skyportal.client import SkyportalClient

matplotlib.use("agg")

logger = logging.getLogger(__name__)


class SkyportalSender(BaseSourceProcessor):
    """
    Processor for sending candidates to Fritz.
    """

    base_key = "skyportalsender"

    def __init__(
        self,
        origin: str,
        group_ids: list[int],
        fritz_filter_id: int,
        instrument_id: int,
        stream_id: int,
        update_thumbnails: bool = True,
    ):
        super().__init__()
        self.group_ids = group_ids
        self.fritz_filter_id = fritz_filter_id
        self.instrument_id = instrument_id
        self.origin = origin  # used for sending updates to Fritz
        self.stream_id = stream_id
        self.update_thumbnails = update_thumbnails
        self.skyportal_client = SkyportalClient()

    def _apply_to_sources(
        self,
        batch: SourceBatch,
    ) -> SourceBatch:
        """
        Apply the processor to a batch of candidates.

        :param batch: SourceBatch to process
        :return: SourceBatch after processing
        """
        for source_table in batch:
            self.export_candidates_to_skyportal(source_table)
        return batch

    @staticmethod
    def read_input_df(candidate_df: pd.DataFrame):
        """Takes a DataFrame, which has multiple candidate
        and creates list of dictionaries, each dictionary
        representing a single candidate.

        Args:
            candidate_df (pandas.core.frame.DataFrame): dataframe of all candidates.

        Returns:
            (list[dict]): list of dictionaries, each a candidate.
        """
        all_candidates = []

        for i in range(0, len(candidate_df)):
            candidate = {}
            for key in candidate_df.keys():
                try:
                    if isinstance(candidate_df.iloc[i].get(key), (list, str)):
                        candidate[key] = candidate_df.iloc[i].get(key)
                    else:
                        # change to native python type
                        candidate[key] = candidate_df.iloc[i].get(key).item()
                except AttributeError:  # for IOBytes objs
                    candidate[key] = candidate_df.iloc[i].get(key)

            all_candidates.append(candidate)

        return all_candidates

    def skyportal_post_source(self, alert: dict, group_ids: Optional[list[int]] = None):
        """Add a new source to SkyPortal

        :param alert: dict of source info
        :param group_ids: list of group_ids to post source to, defaults to None
        :return: None
        """
        if group_ids is None:
            group_ids = self.group_ids

        data = {
            "ra": alert["ra"],
            "dec": alert["dec"],
            "id": alert[CAND_NAME_KEY],
            "group_ids": group_ids,
            "origin": self.origin,
        }

        logger.debug(
            f"Saving {alert[CAND_NAME_KEY]} {alert['candid']} as a Source on SkyPortal"
        )
        response = self.api("POST", "sources", data)

        if response.json()["status"] == "success":
            logger.debug(
                f"Saved {alert[CAND_NAME_KEY]} {alert['candid']} "
                f"as a Source on SkyPortal"
            )
        else:
            err = (
                f"Failed to save {alert[CAND_NAME_KEY]} {alert['candid']} "
                f"as a Source on SkyPortal"
            )
            logger.error(err)
            logger.error(response.json())

    def api(
        self, method: str, endpoint: str, data: Optional[Mapping] = None
    ) -> requests.Response:
        """Make an API call to a SkyPortal instance

        headers = {'Authorization': f'token {self.token}'}
        response = requests.request(method, endpoint, json_dict=data, headers=headers)

        :param method: HTTP method
        :param endpoint: API endpoint e.g sources
        :param data: JSON data to send
        :return: response from API call
        """
        return self.skyportal_client.api(method, endpoint, data)

    def skyportal_post_candidate(self, alert):
        """
        Post a candidate on SkyPortal. Creates new candidate(s) (one per filter)

        :param alert: dict of alert info
        :return: None
        """

        data = {
            "id": alert[CAND_NAME_KEY],
            "ra": alert["ra"],
            "dec": alert["dec"],
            "filter_ids": [self.fritz_filter_id],
            "passing_alert_id": self.fritz_filter_id,
            "passed_at": Time(datetime.utcnow()).isot,
            "origin": f"{PACKAGE_NAME}:{__version__}",
        }

        logger.debug(
            f"Posting metadata of {alert[CAND_NAME_KEY]} {alert['candid']} to SkyPortal"
        )
        response = self.api("POST", "candidates", data)

        if response.json()["status"] == "success":
            logger.debug(
                f"Posted {alert[CAND_NAME_KEY]} {alert['candid']} metadata to SkyPortal"
            )
        else:
            logger.error(
                f"Failed to post {alert[CAND_NAME_KEY]} {alert['candid']} "
                f"metadata to SkyPortal"
            )
            logger.error(response.json())

    def make_thumbnail(self, alert, skyportal_type: str, alert_packet_type: str):
        """
        Convert lossless FITS cutouts from ZTF-like alerts into PNGs.
        Make thumbnail for pushing to SkyPortal.

        :param alert: ZTF-like alert packet/dict
        :param skyportal_type: <new|ref|sub> thumbnail type expected by SkyPortal
        :param alert_packet_type: <Science|Template|Difference> survey naming
        :return:
        """
        alert = deepcopy(alert)
        cutout_data = alert[f"cutout_{alert_packet_type}"]

        with gzip.open(io.BytesIO(cutout_data), "rb") as cutout:
            with fits.open(
                io.BytesIO(cutout.read()), ignore_missing_simple=True
            ) as hdu:
                image_data = hdu[0].data  # pylint: disable=no-member

        buff = io.BytesIO()
        plt.close("all")
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
            stretch=LinearStretch()
            if alert_packet_type == "Difference"
            else LogStretch(),
        )
        img_norm = norm(img)
        normalizer = AsymmetricPercentileInterval(
            lower_percentile=1, upper_percentile=100
        )
        vmin, vmax = normalizer.get_limits(img_norm)
        ax_1.imshow(img_norm, cmap="bone", origin="lower", vmin=vmin, vmax=vmax)
        plt.savefig(buff, dpi=42)

        buff.seek(0)
        plt.close("all")

        thumbnail_dict = {
            "obj_id": alert[CAND_NAME_KEY],
            "data": base64.b64encode(buff.read()).decode("utf-8"),
            "ttype": skyportal_type,
        }

        return thumbnail_dict

    def skyportal_post_thumbnails(self, alert):
        """Post alert Science, Reference, and Subtraction thumbnails to SkyPortal

        :param alert: dict of source/candidate information
        :return:
        """
        for ttype, instrument_type in [
            ("new", "science"),
            ("ref", "template"),
            ("sub", "difference"),
        ]:
            logger.debug(
                f"Making {instrument_type} thumbnail for {alert[CAND_NAME_KEY]} "
                f"{alert['candid']}",
            )
            thumb = self.make_thumbnail(alert, ttype, instrument_type)

            logger.debug(
                f"Posting {instrument_type} thumbnail for {alert[CAND_NAME_KEY]} "
                f"{alert['candid']} to SkyPortal",
            )
            response = self.api("POST", "thumbnail", thumb)

            if response.json()["status"] == "success":
                logger.debug(
                    f"Posted {alert[CAND_NAME_KEY]} {alert['candid']} "
                    f"{instrument_type} cutout to SkyPortal"
                )
            else:
                logger.error(
                    f"Failed to post {alert[CAND_NAME_KEY]} {alert['candid']} "
                    f"{instrument_type} cutout to SkyPortal"
                )
                logger.error(response.json())

    def upload_thumbnail(self, alert):
        """Post new thumbnail to Fritz.

        NOTE: this is the original WINTER method for sending thumbnails,
        not full sized but higher contrast, similar to alert_make_thumbnail

        Format of thumbnail payload:
        { "obj_id": "string",  "data": "string",  "ttype": "string"}
        """
        fritz_to_cand = {"new": "SciBitIm", "ref": "RefBitIm", "sub": "DiffBitIm"}

        for fritz_key, cand_key in fritz_to_cand.items():
            cutout = alert[cand_key]

            buffer = io.BytesIO()
            plt.figure(figsize=(3, 3))
            mean, median, std = sigma_clipped_stats(cutout)
            plt.imshow(
                cutout,
                origin="lower",
                cmap="gray",
                vmin=mean - 1 * std,
                vmax=median + 3 * std,
            )
            plt.xticks([])
            plt.yticks([])

            plt.savefig(buffer, format="png")

            cutoutb64 = base64.b64encode(buffer.getvalue())
            cutoutb64_string = cutoutb64.decode("utf8")

            data_payload = {
                "obj_id": alert[CAND_NAME_KEY],
                "data": cutoutb64_string,
                "ttype": fritz_key,
            }

            response = self.api("POST", "thumbnail", data=data_payload)

            if response.json()["status"] == "success":
                logger.debug(
                    f"Posted {alert[CAND_NAME_KEY]} {alert['candid']} "
                    f"{cand_key} cutout to SkyPortal"
                )
            else:
                logger.error(
                    f"Failed to post {alert[CAND_NAME_KEY]} {alert['candid']} "
                    f"{cand_key} cutout to SkyPortal"
                )
                logger.error(response.json())

    def make_photometry(self, alert, jd_start: Optional[float] = None):
        """
        Make a de-duplicated pandas.DataFrame with photometry of alert[CAND_NAME_KEY]
        Modified from Kowalksi (https://github.com/dmitryduev/kowalski)

        :param alert: candidate dictionary
        :param jd_start: date from which to start photometry from
        """

        photometry_table = [
            {
                "mjd": alert["mjd"],
                "mag": alert["magpsf"],
                "magerr": alert["sigmapsf"],
                "filter": alert["sncosmofilter"],
                "ra": alert["ra"],
                "dec": alert["dec"],
            }
        ]

        if len(alert[SOURCE_HISTORY_KEY]) > 0:
            prv_detections = pd.DataFrame.from_records(alert[SOURCE_HISTORY_KEY])

            for _, row in prv_detections.iterrows():
                try:
                    photometry_table.append(
                        {
                            "mjd": row["mjd"],
                            "mag": row["magpsf"],
                            "magerr": row["sigmapsf"],
                            "filter": row["sncosmofilter"],
                            "ra": row["ra"],
                            "dec": row["dec"],
                        }
                    )
                except KeyError:
                    logger.warning(
                        f"Missing photometry information for previous "
                        f"detection of {alert[CAND_NAME_KEY]}"
                    )

        df_photometry = pd.DataFrame(photometry_table)

        df_photometry["magsys"] = "ab"
        df_photometry["limiting_mag"] = 23.9

        df_photometry["obj_id"] = alert[CAND_NAME_KEY]
        df_photometry["stream_ids"] = [int(self.stream_id)]
        df_photometry["instrument_id"] = self.instrument_id

        return df_photometry

    def skyportal_put_photometry(self, alert):
        """Send photometry to Fritz."""
        logger.debug(
            f"Making alert photometry of {alert[CAND_NAME_KEY]} {alert['candid']}"
        )
        df_photometry = self.make_photometry(alert)

        # post photometry
        photometry = df_photometry.to_dict("list")

        if len(photometry) > 0:
            logger.debug(
                f"Posting photometry of {alert[CAND_NAME_KEY]} {alert['candid']}, "
                f"stream_id={self.stream_id} to SkyPortal"
            )
            response = self.api("PUT", "photometry", photometry)
            if response.json()["status"] == "success":
                logger.debug(
                    f"Posted {alert[CAND_NAME_KEY]} photometry stream_id={self.stream_id} "
                    f"to SkyPortal"
                )
            else:
                logger.error(
                    f"Failed to post {alert[CAND_NAME_KEY]} photometry "
                    f"stream_id={self.stream_id} to SkyPortal"
                )
                logger.error(response.json())

    def skyportal_post_annotation(self, alert):
        """
        Post an annotation. Works for both candidates and sources.

        :param alert: alert data
        :return: None
        """
        data = {
            "chipsf": alert["chipsf"],
            "fwhm": alert["fwhm"],
            "scorr": alert["scorr"],
        }
        payload = {"origin": self.origin, "data": data, "group_ids": self.group_ids}

        path = f"sources/{str(alert[CAND_NAME_KEY])}/annotations"
        response = self.api("POST", path, payload)

        if response.json()["status"] == "success":
            logger.debug(f"Posted {alert[CAND_NAME_KEY]} annotation to SkyPortal")
        else:
            logger.error(
                f"Failed to post {alert[CAND_NAME_KEY]} annotation to SkyPortal"
            )
            logger.error(response.json())

    def skyportal_put_annotation(self, source):
        """
        Retrieve an annotation to check if it exists already.

        :param source: detection data
        :return: None
        """
        response = self.api(
            "GET",
            f"sources/{str(source[CAND_NAME_KEY])}/annotations",
        )

        if response.json()["status"] == "success":
            logger.debug(f"Got {source[CAND_NAME_KEY]} annotations from SkyPortal")
        else:
            logger.debug(
                f"Failed to get {source[CAND_NAME_KEY]} annotations from SkyPortal"
            )
            logger.debug(response.json())
            return

        existing_annotations = {
            annotation["origin"]: {
                "annotation_id": annotation["id"],
                "author_id": annotation["author_id"],
            }
            for annotation in response.json()["data"]
        }

        # no previous annotation exists on SkyPortal for this object? just post then
        if self.origin not in existing_annotations:
            self.skyportal_post_annotation(source)
        # annotation from this(WNTR) origin exists
        else:
            # annotation data
            data = {
                "fwhm": source["fwhm"],
                "scorr": source["scorr"],
                "chipsf": source["chipsf"],
            }
            new_annotation = {
                "author_id": existing_annotations[self.origin]["author_id"],
                "obj_id": source[CAND_NAME_KEY],
                "origin": self.origin,
                "data": data,
                "group_ids": self.group_ids,
            }

            logger.debug(
                f"Putting annotation for {source[CAND_NAME_KEY]} {source['candid']} "
                f"to SkyPortal",
            )
            response = self.api(
                "PUT",
                f"sources/{source[CAND_NAME_KEY]}"
                f"/annotations/{existing_annotations[self.origin]['annotation_id']}",
                new_annotation,
            )
            if response.json()["status"] == "success":
                logger.debug(
                    f"Posted updated {source[CAND_NAME_KEY]} annotation to SkyPortal"
                )
            else:
                logger.error(
                    f"Failed to post updated {source[CAND_NAME_KEY]} annotation "
                    f"to SkyPortal"
                )
                logger.error(response.json())

    def skyportal_candidate_exporter(self, alert):
        """
        Posts a candidate to SkyPortal.

        :param alert: _description_
        :type alert: _type_
        """
        # check if candidate exists in SkyPortal
        logger.debug(f"Checking if {alert[CAND_NAME_KEY]} is candidate in SkyPortal")
        response = self.api("HEAD", f"candidates/{alert[CAND_NAME_KEY]}")

        if response.status_code not in [200, 404]:
            response.raise_for_status()

        is_candidate = response.status_code == 200
        logger.debug(
            f"{alert[CAND_NAME_KEY]} {'is' if is_candidate else 'is not'} "
            f"candidate in SkyPortal"
        )

        # check if source exists in SkyPortal
        logger.debug(f"Checking if {alert[CAND_NAME_KEY]} is source in SkyPortal")
        response = self.api("HEAD", f"sources/{alert[CAND_NAME_KEY]}")

        if response.status_code not in [200, 404]:
            response.raise_for_status()

        is_source = response.status_code == 200
        logger.debug(
            f"{alert[CAND_NAME_KEY]} "
            f"{'is' if is_source else 'is not'} source in SkyPortal"
        )

        # object does not exist in SkyPortal: neither cand nor source
        if (not is_candidate) and (not is_source):
            # post candidate
            self.skyportal_post_candidate(alert)

            # post annotations
            self.skyportal_post_annotation(alert)

            # post full light curve
            self.skyportal_put_photometry(alert)

            # post thumbnails
            self.skyportal_post_thumbnails(alert)

        # obj already exists in SkyPortal
        else:
            # post candidate with new filter ids
            self.skyportal_post_candidate(alert)

            # put (*not* post) annotations
            self.skyportal_put_annotation(alert)

            # exists in SkyPortal & already saved as a source
            if is_source:
                # get info on the corresponding groups:
                logger.debug(
                    f"Getting source groups info on "
                    f"{alert[CAND_NAME_KEY]} from SkyPortal",
                )
                response = self.api(
                    "GET",
                    f"sources/{alert[CAND_NAME_KEY]}/groups",
                )
                if response.json()["status"] == "success":
                    existing_groups = response.json()["data"]
                    existing_group_ids = [g["id"] for g in existing_groups]

                    for existing_gid in existing_group_ids:
                        if existing_gid in self.group_ids:
                            self.skyportal_post_source(alert, [existing_gid])
                else:
                    logger.error(
                        f"Failed to get source groups info on {alert[CAND_NAME_KEY]}"
                    )
            else:  # exists in SkyPortal but NOT saved as a source
                self.skyportal_post_source(alert)

            # post alert photometry in single call to /api/photometry
            self.skyportal_put_photometry(alert)

            if self.update_thumbnails:
                self.skyportal_post_thumbnails(alert)

        logger.debug(f"SendToSkyportal Manager complete for {alert[CAND_NAME_KEY]}")

    def export_candidates_to_skyportal(self, source_table: SourceTable):
        """
        Function to export individual sources as candidates in SkyPortal

        :param source_table: Table containing the data to be processed
        :return: None
        """
        candidate_df = source_table.get_data()
        t_0 = time.time()
        # all_cands = self.read_input_df(candidate_df)
        # num_cands = len(all_cands)

        metadata = source_table.get_metadata()

        for _, src in candidate_df.iterrows():
            cand = self.generate_super_dict(metadata, src)
            self.skyportal_candidate_exporter(deepcopy(cand))

        t_1 = time.time()
        logger.debug(
            f"Took {(t_1 - t_0):.2f} seconds to Skyportal process "
            f"{len(candidate_df)} candidates."
        )
