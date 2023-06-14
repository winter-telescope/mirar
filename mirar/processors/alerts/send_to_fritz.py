"""
Module for sending candidates to Fritz.
"""
import base64
import gzip
import io
import logging
import os
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
from requests.adapters import HTTPAdapter
from urllib3.util import Retry

from mirar.data import SourceBatch
from mirar.paths import PACKAGE_NAME, __version__
from mirar.processors.base_processor import BaseDataframeProcessor

matplotlib.use("agg")

logger = logging.getLogger(__name__)

DEFAULT_TIMEOUT = 5  # seconds


class TimeoutHTTPAdapter(HTTPAdapter):
    """
    HTTP adapter that sets a default timeout for all requests.
    """

    def __init__(self, *args, **kwargs):
        self.timeout = DEFAULT_TIMEOUT
        if "timeout" in kwargs:
            self.timeout = kwargs["timeout"]
            del kwargs["timeout"]
        super().__init__(*args, **kwargs)

    def send(self, request, *args, **kwargs):
        """
        Send a request with a default timeout.

        :param request: request to send
        :param args: args to pass to super().send
        :param kwargs: kwargs to pass to super().send
        :return: response from request
        """
        try:
            timeout = kwargs.get("timeout")
            if timeout is None:
                kwargs["timeout"] = self.timeout
            return super().send(request, *args, **kwargs)
        except AttributeError:
            kwargs["timeout"] = DEFAULT_TIMEOUT


class SendToFritz(BaseDataframeProcessor):
    """
    Processor for sending candidates to Fritz.
    """

    base_key = "fritzsender"

    def __init__(
        self,
        base_name: str,
        group_ids: list[int],
        fritz_filter_id: int,
        instrument_id: int,
        stream_id: int,
        update_thumbnails: bool = True,
        protocol: str = "http",
    ):
        super().__init__()
        self.token = None
        self.group_ids = group_ids
        self.fritz_filter_id = fritz_filter_id
        self.instrument_id = instrument_id
        self.origin = base_name  # used for sending updates to Fritz
        self.stream_id = stream_id
        self.protocol = protocol
        self.update_thumbnails = update_thumbnails

        self._session = None
        self.session_headers = None

    def set_up_session(self):
        """
        Set up a session for sending requests to Fritz.

        :return: None
        """
        # session to talk to SkyPortal/Fritz
        self._session = requests.Session()
        self.session_headers = {
            "Authorization": f"token {self._get_fritz_token()}",
            "User-Agent": "mirar",
        }

        retries = Retry(
            total=5,
            backoff_factor=2,
            status_forcelist=[405, 429, 500, 502, 503, 504],
            method_whitelist=["HEAD", "GET", "PUT", "POST", "PATCH"],
        )
        adapter = TimeoutHTTPAdapter(timeout=5, max_retries=retries)
        self._session.mount("https://", adapter)
        self._session.mount("http://", adapter)

    def get_session(self) -> requests.Session:
        """
        Wrapper for getting the session.
        If the session is not set up, it will be set up.

        :return: Session
        """
        if self._session is None:
            self.set_up_session()

        return self._session

    @staticmethod
    def _get_fritz_token():
        """
        Get Fritz token from environment variable.

        :return: Fritz token
        """
        token_fritz = os.getenv("FRITZ_TOKEN")

        if token_fritz is None:
            err = (
                "No Fritz token specified. Run 'export FRITZ_TOKEN=<token>' to "
                "set. The Fritz token will need to be specified manually "
                "for Fritz API queries."
            )
            logger.error(err)
            raise ValueError(err)

        return token_fritz

    def _apply_to_candidates(
        self,
        batch: SourceBatch,
    ) -> SourceBatch:
        """
        Apply the processor to a batch of candidates.

        :param batch: SourceBatch to process
        :return: SourceBatch after processing
        """
        for source_table in batch:
            candidate_table = source_table.get_data()
            self.make_alert(candidate_table)
        return batch

    @staticmethod
    def _get_author_id():
        """
        Fritz author id is used in update calls. Can be found

        :return: Fritz author id
        """
        authid_fritz = os.getenv("FRITZ_AUTHID")

        if authid_fritz is None:
            err = (
                "No Fritz author id specified. Run 'export FRITZ_AUTHID=<id>' to set. "
                "Author id needs to be specified for updates sent by Fritz API queries."
            )
            logger.error(err)
            raise ValueError(err)

        return authid_fritz

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
                    candidate[key] = candidate_df.iloc[i].get(key).getvalue()

            all_candidates.append(candidate)

        return all_candidates

    def api(
        self, method: str, endpoint: str, data: Optional[Mapping] = None
    ) -> requests.Response:
        """Make an API call to a SkyPortal instance

        headers = {'Authorization': f'token {self.token}'}
        response = requests.request(method, endpoint, json_dict=data, headers=headers)

        :param method: HTTP method
        :param endpoint: API endpoint
        :param data: JSON data to send
        :return: response from API call
        """
        method = method.lower()

        session = self.get_session()

        methods = {
            "head": session.head,
            "get": session.get,
            "post": session.post,
            "put": session.put,
            "patch": session.patch,
            "delete": session.delete,
        }

        if endpoint is None:
            raise ValueError("Endpoint not specified")
        if method not in ["head", "get", "post", "put", "patch", "delete"]:
            raise ValueError(f"Unsupported method: {method}")

        if method == "get":
            response = methods[method](
                f"{endpoint}",
                params=data,
                headers=self.session_headers,
            )
        else:
            response = methods[method](
                f"{endpoint}",
                json=data,
                headers=self.session_headers,
            )

        return response

    def alert_post_source(self, alert: dict, group_ids: Optional[list[int]] = None):
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
            "id": alert["objectId"],
            "group_ids": group_ids,
            "origin": self.origin,
        }

        logger.debug(
            f"Saving {alert['objectId']} {alert['candid']} as a Source on SkyPortal"
        )
        response = self.api("POST", "https://fritz.science/api/sources", data)

        if response.json()["status"] == "success":
            logger.debug(
                f"Saved {alert['objectId']} {alert['candid']} as a Source on SkyPortal"
            )
        else:
            err = (
                f"Failed to save {alert['objectId']} {alert['candid']} "
                f"as a Source on SkyPortal"
            )
            logger.error(err)
            logger.error(response.json())

    def alert_post_candidate(self, alert):
        """
        Post a candidate on SkyPortal. Creates new candidate(s) (one per filter)

        :param alert: dict of alert info
        :return: None
        """

        data = {
            "id": alert["objectId"],
            "ra": alert["ra"],
            "dec": alert["dec"],
            "filter_ids": [self.fritz_filter_id],
            "passing_alert_id": self.fritz_filter_id,
            "passed_at": Time(datetime.utcnow()).isot,
            "origin": f"{PACKAGE_NAME}:{__version__}",
        }

        logger.debug(
            f"Posting metadata of {alert['objectId']} {alert['candid']} to SkyPortal"
        )
        response = self.api("POST", "https://fritz.science/api/candidates", data)

        if response.json()["status"] == "success":
            logger.debug(
                f"Posted {alert['objectId']} {alert['candid']} metadata to SkyPortal"
            )
        else:
            logger.error(
                f"Failed to post {alert['objectId']} {alert['candid']} "
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
        cutout_data = alert[f"cutout{alert_packet_type}"]

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
            "obj_id": alert["objectId"],
            "data": base64.b64encode(buff.read()).decode("utf-8"),
            "ttype": skyportal_type,
        }

        return thumbnail_dict

    def alert_post_thumbnails(self, alert):
        """Post alert Science, Reference, and Subtraction thumbnails to SkyPortal

        :param alert: dict of source/candidate information
        :return:
        """
        for ttype, instrument_type in [
            ("new", "Science"),
            ("ref", "Template"),
            ("sub", "Difference"),
        ]:
            logger.debug(
                f"Making {instrument_type} thumbnail for {alert['objectId']} "
                f"{alert['candid']}",
            )
            thumb = self.make_thumbnail(alert, ttype, instrument_type)

            logger.debug(
                f"Posting {instrument_type} thumbnail for {alert['objectId']} "
                f"{alert['candid']} to SkyPortal",
            )
            response = self.api("POST", "https://fritz.science/api/thumbnail", thumb)

            if response.json()["status"] == "success":
                logger.debug(
                    f"Posted {alert['objectId']} {alert['candid']} "
                    f"{instrument_type} cutout to SkyPortal"
                )
            else:
                logger.error(
                    f"Failed to post {alert['objectId']} {alert['candid']} "
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
                "obj_id": alert["objectId"],
                "data": cutoutb64_string,
                "ttype": fritz_key,
            }

            response = self.api(
                "POST", "https://fritz.science/api/thumbnail", data=data_payload
            )

            if response.json()["status"] == "success":
                logger.debug(
                    f"Posted {alert['objectId']} {alert['candid']} "
                    f"{cand_key} cutout to SkyPortal"
                )
            else:
                logger.error(
                    f"Failed to post {alert['objectId']} {alert['candid']} "
                    f"{cand_key} cutout to SkyPortal"
                )
                logger.error(response.json())

    def make_photometry(self, alert, jd_start: Optional[float] = None):
        """
        Make a de-duplicated pandas.DataFrame with photometry of alert['objectId']
        Modified from Kowalksi (https://github.com/dmitryduev/kowalski)

        :param alert: candidate dictionary
        :param jd_start: date from which to start photometry from
        """
        alert = deepcopy(alert)
        top_level = [
            "schemavsn",
            "publisher",
            "objectId",
            "candid",
            "candidate",
            "prv_candidates",
            "cutoutScience",
            "cutoutTemplate",
            "cutoutDifference",
        ]
        alert["candidate"] = {}

        # (keys having value in 3.)
        delete = [key for key in alert.keys() if key not in top_level]

        # delete the key/s
        for key in delete:
            alert["candidate"][key] = alert[key]
            del alert[key]

        alert["candidate"] = [alert["candidate"]]
        df_candidate = pd.DataFrame(alert["candidate"], index=[0])

        df_prv_candidates = pd.DataFrame(alert["prv_candidates"])

        df_light_curve = pd.concat(
            [df_candidate, df_prv_candidates], ignore_index=True, sort=False
        )

        # note: WNTR (like PGIR) uses 2massj, which is not in sncosmo as of
        # 20210803, cspjs seems to be close/good enough as an approximation
        df_light_curve["filter"] = "cspjs"

        df_light_curve["magsys"] = "ab"
        df_light_curve["mjd"] = df_light_curve["jd"] - 2400000.5

        df_light_curve["mjd"] = df_light_curve["mjd"].astype(np.float64)
        df_light_curve["magpsf"] = df_light_curve["magpsf"].astype(np.float32)
        df_light_curve["sigmapsf"] = df_light_curve["sigmapsf"].astype(np.float32)

        df_light_curve = (
            df_light_curve.drop_duplicates(subset=["mjd", "magpsf"])
            .reset_index(drop=True)
            .sort_values(by=["mjd"])
        )

        # filter out bad data:
        mask_good_diffmaglim = df_light_curve["diffmaglim"] > 0
        df_light_curve = df_light_curve.loc[mask_good_diffmaglim]

        # convert from mag to flux

        # step 1: calculate the coefficient that determines whether the
        # flux should be negative or positive
        coeff = df_light_curve["isdiffpos"].apply(
            lambda x: 1.0 if x in [True, 1, "y", "Y", "t", "1"] else -1.0
        )

        # step 2: calculate the flux normalized to an arbitrary AB zeropoint of
        # 23.9 (results in flux in uJy)
        df_light_curve["flux"] = coeff * 10 ** (
            -0.4 * (df_light_curve["magpsf"] - 23.9)
        )

        # step 3: separate detections from non detections
        detected = np.isfinite(df_light_curve["magpsf"])
        undetected = ~detected

        # step 4: calculate the flux error
        df_light_curve["fluxerr"] = None  # initialize the column

        # step 4a: calculate fluxerr for detections using sigmapsf
        df_light_curve.loc[detected, "fluxerr"] = np.abs(
            df_light_curve.loc[detected, "sigmapsf"]
            * df_light_curve.loc[detected, "flux"]
            * np.log(10)
            / 2.5
        )

        # step 4b: calculate fluxerr for non detections using diffmaglim
        df_light_curve.loc[undetected, "fluxerr"] = (
            10 ** (-0.4 * (df_light_curve.loc[undetected, "diffmaglim"] - 23.9)) / 5.0
        )  # as diffmaglim is the 5-sigma depth

        # step 5: set the zeropoint and magnitude system
        df_light_curve["zp"] = 23.9
        df_light_curve["zpsys"] = "ab"

        # only "new" photometry requested?
        if jd_start is not None:
            w_after_jd = df_light_curve["jd"] > jd_start
            df_light_curve = df_light_curve.loc[w_after_jd]

        return df_light_curve

    def alert_put_photometry(self, alert):
        """Send photometry to Fritz."""
        logger.debug(
            f"Making alert photometry of {alert['objectId']} {alert['candid']}"
        )
        df_photometry = self.make_photometry(alert)

        # post photometry
        photometry = {
            "obj_id": alert["objectId"],
            "stream_ids": [int(self.stream_id)],
            "instrument_id": self.instrument_id,
            "mjd": df_photometry["mjd"].tolist(),
            "flux": df_photometry["flux"].tolist(),
            "fluxerr": df_photometry["fluxerr"].tolist(),
            "zp": df_photometry["zp"].tolist(),
            "magsys": df_photometry["zpsys"].tolist(),
            "filter": df_photometry["filter"].tolist(),
            "ra": df_photometry["ra"].tolist(),
            "dec": df_photometry["dec"].tolist(),
        }

        if (len(photometry.get("flux", ())) > 0) or (
            len(photometry.get("fluxerr", ())) > 0
        ):
            logger.debug(
                f"Posting photometry of {alert['objectId']} {alert['candid']}, "
                f"stream_id={self.stream_id} to SkyPortal"
            )
            response = self.api(
                "PUT", "https://fritz.science/api/photometry", photometry
            )
            if response.json()["status"] == "success":
                logger.debug(
                    f"Posted {alert['objectId']} photometry stream_id={self.stream_id} "
                    f"to SkyPortal"
                )
            else:
                logger.error(
                    f"Failed to post {alert['objectId']} photometry "
                    f"stream_id={self.stream_id} to SkyPortal"
                )
                logger.error(response.json())

    def alert_post_annotation(self, alert):
        """Post an annotation. Works for both candidates and sources."""
        data = {
            "chipsf": alert["chipsf"],
            "fwhm": alert["fwhm"],
            "scorr": alert["scorr"],
        }
        payload = {"origin": self.origin, "data": data, "group_ids": self.group_ids}

        path = f'https://fritz.science/api/sources/{str(alert["objectId"])}/annotations'
        response = self.api("POST", path, payload)

        if response.json()["status"] == "success":
            logger.debug(f"Posted {alert['objectId']} annotation to SkyPortal")
        else:
            logger.error(f"Failed to post {alert['objectId']} annotation to SkyPortal")
            logger.error(response.json())

    def alert_put_annotation(self, alert):
        """Retrieve an annotation to check if it exists already."""
        response = self.api(
            "GET",
            f'https://fritz.science/api/sources/{str(alert["objectId"])}/annotations',
        )

        if response.json()["status"] == "success":
            logger.debug(f"Got {alert['objectId']} annotations from SkyPortal")
        else:
            logger.debug(
                f"Failed to get {alert['objectId']} annotations from SkyPortal"
            )
            logger.debug(response.json())
            return False

        existing_annotations = {
            annotation["origin"]: {
                "annotation_id": annotation["id"],
                "author_id": annotation["author_id"],
            }
            for annotation in response.json()["data"]
        }

        # no previous annotation exists on SkyPortal for this object? just post then
        if self.origin not in existing_annotations:
            self.alert_post_annotation(alert)
        # annotation from this(WNTR) origin exists
        else:
            # annotation data
            data = {
                "fwhm": alert["fwhm"],
                "scorr": alert["scorr"],
                "chipsf": alert["chipsf"],
            }
            new_annotation = {
                "author_id": existing_annotations[self.origin]["author_id"],
                "obj_id": alert["objectId"],
                "origin": self.origin,
                "data": data,
                "group_ids": self.group_ids,
            }

            logger.debug(
                f"Putting annotation for {alert['objectId']} {alert['candid']} "
                f"to SkyPortal",
            )
            response = self.api(
                "PUT",
                f"https://fritz.science/api/sources/{alert['objectId']}"
                f"/annotations/{existing_annotations[self.origin]['annotation_id']}",
                new_annotation,
            )
            if response.json()["status"] == "success":
                logger.debug(
                    f"Posted updated {alert['objectId']} annotation to SkyPortal"
                )
            else:
                logger.error(
                    f"Failed to post updated {alert['objectId']} annotation "
                    f"to SkyPortal"
                )
                logger.error(response.json())

        return None

    def alert_skyportal_manager(self, alert):
        """Posts alerts to SkyPortal if criteria is met

        :param alert: _description_
        :type alert: _type_
        """
        # check if candidate exists in SkyPortal
        logger.debug(f"Checking if {alert['objectId']} is candidate in SkyPortal")
        response = self.api(
            "HEAD", f"https://fritz.science/api/candidates/{alert['objectId']}"
        )
        is_candidate = response.status_code == 200
        logger.debug(
            f"{alert['objectId']} {'is' if is_candidate else 'is not'} "
            f"candidate in SkyPortal"
        )

        # check if source exists in SkyPortal
        logger.debug(f"Checking if {alert['objectId']} is source in SkyPortal")
        response = self.api(
            "HEAD", f"https://fritz.science/api/sources/{alert['objectId']}"
        )
        is_source = response.status_code == 200
        logger.debug(
            f"{alert['objectId']} {'is' if is_source else 'is not'} source in SkyPortal"
        )

        # object does not exist in SkyPortal: neither cand nor source
        if (not is_candidate) and (not is_source):
            # post candidate
            self.alert_post_candidate(alert)

            # post annotations
            self.alert_post_annotation(alert)

            # post full light curve
            self.alert_put_photometry(alert)

            # post thumbnails
            self.alert_post_thumbnails(alert)

            # TODO autosave stuff, necessary?

        # obj already exists in SkyPortal
        else:
            # TODO passed_filters logic

            # post candidate with new filter ids
            self.alert_post_candidate(alert)

            # put (*not* post) annotations
            self.alert_put_annotation(alert)

            # exists in SkyPortal & already saved as a source
            if is_source:
                # get info on the corresponding groups:
                logger.debug(
                    f"Getting source groups info on {alert['objectId']} from SkyPortal",
                )
                response = self.api(
                    "GET",
                    f"https://fritz.science/api/sources/{alert['objectId']}/groups",
                )
                if response.json()["status"] == "success":
                    existing_groups = response.json()["data"]
                    existing_group_ids = [g["id"] for g in existing_groups]

                    for existing_gid in existing_group_ids:
                        if existing_gid in self.group_ids:
                            self.alert_post_source(alert, [existing_gid])
                else:
                    logger.error(
                        f"Failed to get source groups info on {alert['objectId']}"
                    )
            else:  # exists in SkyPortal but NOT saved as a source
                self.alert_post_source(alert)

            # post alert photometry in single call to /api/photometry
            self.alert_put_photometry(alert)

            if self.update_thumbnails:
                self.alert_post_thumbnails(alert)

        logger.debug(f'SendToFritz Manager complete for {alert["objectId"]}')

    def make_alert(self, candidate_df: pd.DataFrame):
        """
        Function to make an alert from a single row of a pandas DataFrame

        :param candidate_df: Candidate DataFrame
        :return: None
        """
        t_0 = time.time()
        all_cands = self.read_input_df(candidate_df)
        num_cands = len(all_cands)

        for cand in all_cands:
            self.alert_skyportal_manager(cand)

        t_1 = time.time()
        logger.info(
            f"Took {(t_1 - t_0):.2f} seconds to Fritz process {num_cands} candidates."
        )
