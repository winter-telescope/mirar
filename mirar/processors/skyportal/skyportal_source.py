"""
Module for sending sources to Fritz.
"""

import logging
from copy import deepcopy
from typing import Mapping, Optional

import matplotlib
import numpy as np
import pandas as pd
import requests
from astropy.time import Time

from mirar.data import SourceBatch
from mirar.data.utils import decode_img
from mirar.paths import SOURCE_HISTORY_KEY, SOURCE_NAME_KEY, TIME_KEY
from mirar.processors.base_processor import BaseSourceProcessor
from mirar.processors.skyportal.client import SkyportalClient
from mirar.processors.skyportal.thumbnail import make_thumbnail

matplotlib.use("agg")

logger = logging.getLogger(__name__)

SNCOSMO_KEY = "sncosmof"


class SkyportalSourceUploader(BaseSourceProcessor):
    """
    Processor for sending source photometry to Skyportal.
    """

    base_key = "skyportalsender"

    def __init__(
        self,
        origin: str,
        group_ids: list[int],
        instrument_id: int,
        update_thumbnails: bool = False,
        skyportal_client: Optional[SkyportalClient] = SkyportalClient(),
    ):
        super().__init__()
        self.group_ids = group_ids
        self.instrument_id = instrument_id
        self.origin = origin  # used for sending updates to Fritz
        self.update_thumbnails = update_thumbnails
        self.skyportal_client = skyportal_client

    def _apply_to_sources(
        self,
        batch: SourceBatch,
    ) -> SourceBatch:
        """
        Apply the processor to a batch of sources/candidates.

        :param batch: SourceBatch to process
        :return: SourceBatch after processing
        """
        for source_table in batch:
            candidate_df = source_table.get_data()

            metadata = source_table.get_metadata()

            candidate_df["mjd"] = Time(metadata[TIME_KEY]).mjd
            for _, src in candidate_df.iterrows():
                super_dict = self.generate_super_dict(metadata, src)
                self.export_to_skyportal(deepcopy(super_dict))

        return batch

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
            "id": alert[SOURCE_NAME_KEY],
            "group_ids": group_ids,
            "origin": self.origin,
        }

        logger.debug(f"Saving {alert[SOURCE_NAME_KEY]} as a Source on SkyPortal")
        response = self.api("POST", "sources", data)

        if response.json()["status"] == "success":
            logger.debug(f"Saved {alert[SOURCE_NAME_KEY]} as a Source on SkyPortal")
        else:
            err = f"Failed to save {alert[SOURCE_NAME_KEY]} as a Source on SkyPortal"
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

    def make_thumbnail(
        self, source: pd.Series, skyportal_type: str, alert_packet_type: str
    ):
        """
        Convert lossless FITS cutouts from ZTF-like alerts into PNGs.
        Make thumbnail for pushing to SkyPortal.

        :param source: ZTF-like alert packet/dict
        :param skyportal_type: <new|ref|sub> thumbnail type expected by SkyPortal
        :param alert_packet_type: <Science|Template|Difference> survey naming
        :return:
        """
        source = deepcopy(source)
        cutout_data = source[f"cutout_{alert_packet_type}"]

        linear_stretch = alert_packet_type.lower() in ["difference"]

        skyportal_thumbnail = make_thumbnail(
            image_data=decode_img(cutout_data),
            linear_stretch=linear_stretch,
        )

        thumbnail_dict = {
            "obj_id": source[SOURCE_NAME_KEY],
            "data": skyportal_thumbnail,
            "ttype": skyportal_type,
        }

        return thumbnail_dict

    def skyportal_post_thumbnails(self, alert):
        """Post alert Science, Reference, and Subtraction thumbnails to SkyPortal

        :param alert: dict of source/candidate information
        :return: None
        """
        for ttype, instrument_type in [
            ("new", "science"),
            ("ref", "template"),
            ("sub", "difference"),
        ]:
            logger.debug(
                f"Making {instrument_type} thumbnail for {alert[SOURCE_NAME_KEY]} "
            )
            thumb = self.make_thumbnail(alert, ttype, instrument_type)

            logger.debug(
                f"Posting {instrument_type} thumbnail for {alert[SOURCE_NAME_KEY]} "
            )
            response = self.api("POST", "thumbnail", thumb)

            if response.json()["status"] == "success":
                logger.debug(
                    f"Posted {alert[SOURCE_NAME_KEY]} "
                    f"{instrument_type} cutout to SkyPortal"
                )
            else:
                logger.error(
                    f"Failed to post {alert[SOURCE_NAME_KEY]} "
                    f"{instrument_type} cutout to SkyPortal"
                )
                logger.error(response.json())

    def make_photometry(self, source: pd.Series) -> pd.DataFrame:
        """
        Make a de-duplicated pandas.DataFrame with photometry of alert[CAND_NAME_KEY]
        Modified from Kowalksi (https://github.com/dmitryduev/kowalski)

        :param source: source row
        :return: pandas.DataFrame with photometry
        """

        photometry_table = [
            {
                "mjd": source["mjd"],
                "mag": source["magpsf"],
                "magerr": source["sigmapsf"],
                "filter": source[SNCOSMO_KEY],
                "ra": source["ra"],
                "dec": source["dec"],
            }
        ]

        if SOURCE_HISTORY_KEY in source:
            if len(source[SOURCE_HISTORY_KEY]) > 0:
                prv_detections = pd.DataFrame.from_records(source[SOURCE_HISTORY_KEY])

                for _, row in prv_detections.iterrows():
                    try:
                        photometry_table.append(
                            {
                                "mjd": row["mjd"],
                                "mag": row["magpsf"],
                                "magerr": row["sigmapsf"],
                                "filter": row[SNCOSMO_KEY],
                                "ra": row["ra"],
                                "dec": row["dec"],
                            }
                        )
                    except KeyError:
                        logger.warning(
                            f"Missing photometry information for previous "
                            f"detection of {source[SOURCE_NAME_KEY]}"
                        )

        df_photometry = pd.DataFrame(photometry_table)

        df_photometry["magsys"] = "ab"

        # step 1: calculate the coefficient that determines whether the
        # flux should be negative or positive
        if "isdiffpos" in source:
            coeff = 1.0 if source["isdiffpos"] else -1.0
        else:
            coeff = 1.0

        # step 2: calculate the flux normalized to an arbitrary AB zeropoint of
        # 23.9 (results in flux in uJy)
        df_photometry["flux"] = coeff * 10 ** (-0.4 * (df_photometry["mag"] - 23.9))

        # step 4a: calculate fluxerr for detections using sigmapsf
        df_photometry["fluxerr"] = (
            df_photometry["magerr"] * df_photometry["flux"] * np.log(10) / 2.5
        )

        # step 5: set the zeropoint and magnitude system
        df_photometry["zp"] = 23.9

        df_photometry.drop(columns=["mag", "magerr"], inplace=True)

        return df_photometry

    def skyportal_put_photometry(self, alert):
        """Send photometry to Fritz."""
        logger.debug(f"Making alert photometry of {alert[SOURCE_NAME_KEY]}")
        df_photometry = self.make_photometry(alert)

        if len(df_photometry) > 0:
            photometry = df_photometry.to_dict("list")
            photometry["obj_id"] = alert[SOURCE_NAME_KEY]
            photometry["instrument_id"] = self.instrument_id
            if hasattr(self, "stream_id"):
                photometry["stream_ids"] = [int(self.stream_id)]
            logger.debug(f"Posting photometry of {alert[SOURCE_NAME_KEY]} to SkyPortal")
            response = self.api("PUT", "photometry", photometry)
            if response.json()["status"] == "success":
                logger.debug(f"Posted {alert[SOURCE_NAME_KEY]} photometry to SkyPortal")
            else:
                logger.error(
                    f"Failed to post {alert[SOURCE_NAME_KEY]} photometry to SkyPortal"
                )
                logger.error(response.json())

    def export_to_skyportal(self, alert):
        """
        Posts a source to SkyPortal.

        :param alert: _description_
        :type alert: _type_
        """
        # check if source exists in SkyPortal # pylint: disable=duplicate-code
        logger.debug(f"Checking if {alert[SOURCE_NAME_KEY]} is source in SkyPortal")
        response = self.api("HEAD", f"sources/{alert[SOURCE_NAME_KEY]}")

        if response.status_code not in [200, 404]:
            response.raise_for_status()

        is_source = response.status_code == 200
        logger.debug(
            f"{alert[SOURCE_NAME_KEY]} "
            f"{'is' if is_source else 'is not'} source in SkyPortal"
        )

        if not is_source:
            self.skyportal_post_source(alert, group_ids=self.group_ids)
            # post thumbnails
            self.skyportal_post_thumbnails(alert)

        # post full light curve
        self.skyportal_put_photometry(alert)

        if self.update_thumbnails:
            self.skyportal_post_thumbnails(alert)

        logger.debug(f"SendToSkyportal Manager complete for {alert[SOURCE_NAME_KEY]}")
