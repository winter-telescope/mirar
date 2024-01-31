"""
Module for sending candidates to Fritz.
"""

import logging
from datetime import datetime

from astropy.time import Time

from mirar.paths import PACKAGE_NAME, SOURCE_NAME_KEY, __version__
from mirar.processors.skyportal.skyportal_source import SkyportalSourceUploader

logger = logging.getLogger(__name__)


class SkyportalCandidateUploader(SkyportalSourceUploader):
    """
    Processor for sending candidates to Fritz.
    """

    def __init__(
        self,
        *args,
        stream_id: int,
        fritz_filter_id: int,
        **kwargs,
    ):
        super().__init__(*args, **kwargs)
        self.stream_id = stream_id
        self.fritz_filter_id = fritz_filter_id

    def skyportal_post_candidate(self, alert):
        """
        Post a candidate on SkyPortal. Creates new candidate(s) (one per filter)

        :param alert: dict of alert info
        :return: None
        """

        data = {
            "id": alert[SOURCE_NAME_KEY],
            "ra": alert["ra"],
            "dec": alert["dec"],
            "filter_ids": [self.fritz_filter_id],
            "passing_alert_id": self.fritz_filter_id,
            "passed_at": Time(datetime.utcnow()).isot,
            "origin": f"{PACKAGE_NAME}:{__version__}",
        }

        logger.debug(
            f"Posting metadata of {alert[SOURCE_NAME_KEY]} "
            f"{alert['candid']} to SkyPortal"
        )
        response = self.api("POST", "candidates", data)

        if response.json()["status"] == "success":
            logger.debug(
                f"Posted {alert[SOURCE_NAME_KEY]} {alert['candid']} "
                f"metadata to SkyPortal"
            )
        else:
            logger.error(
                f"Failed to post {alert[SOURCE_NAME_KEY]} {alert['candid']} "
                f"metadata to SkyPortal"
            )
            logger.error(response.json())

    def skyportal_post_annotation(self, alert):
        """
        Post an annotation. Works for both candidates and sources.

        :param alert: alert data
        :return: None
        """
        data = {}

        for key in [
            "chipsf",
            "fwhm",
            "scorr",
            "nneg",
            "mindtoedge",
            "diffmaglim",
            "distpsnr1",
            "sgmag1",
            "srmag1",
            "simag1",
        ]:
            if key in alert:
                data[key] = alert[key]

        payload = {"origin": self.origin, "data": data, "group_ids": self.group_ids}

        path = f"sources/{str(alert[SOURCE_NAME_KEY])}/annotations"
        response = self.api("POST", path, payload)

        if response.json()["status"] == "success":
            logger.debug(f"Posted {alert[SOURCE_NAME_KEY]} annotation to SkyPortal")
        else:
            logger.error(
                f"Failed to post {alert[SOURCE_NAME_KEY]} annotation to SkyPortal"
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
            f"sources/{str(source[SOURCE_NAME_KEY])}/annotations",
        )

        if response.json()["status"] == "success":
            logger.debug(f"Got {source[SOURCE_NAME_KEY]} annotations from SkyPortal")
        else:
            logger.debug(
                f"Failed to get {source[SOURCE_NAME_KEY]} annotations from SkyPortal"
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
                "obj_id": source[SOURCE_NAME_KEY],
                "origin": self.origin,
                "data": data,
                "group_ids": self.group_ids,
            }

            logger.debug(
                f"Putting annotation for {source[SOURCE_NAME_KEY]} {source['candid']} "
                f"to SkyPortal",
            )
            response = self.api(
                "PUT",
                f"sources/{source[SOURCE_NAME_KEY]}"
                f"/annotations/{existing_annotations[self.origin]['annotation_id']}",
                new_annotation,
            )
            if response.json()["status"] == "success":
                logger.debug(
                    f"Posted updated {source[SOURCE_NAME_KEY]} annotation to SkyPortal"
                )
            else:
                logger.error(
                    f"Failed to post updated {source[SOURCE_NAME_KEY]} annotation "
                    f"to SkyPortal"
                )
                logger.error(response.json())

    def export_to_skyportal(self, alert):
        """
        Posts a candidate to SkyPortal.

        :param alert: _description_
        :type alert: _type_
        """
        # check if candidate exists in SkyPortal
        logger.debug(f"Checking if {alert[SOURCE_NAME_KEY]} is candidate in SkyPortal")
        response = self.api("HEAD", f"candidates/{alert[SOURCE_NAME_KEY]}")

        if response.status_code not in [200, 404]:
            response.raise_for_status()

        is_candidate = response.status_code == 200
        logger.debug(
            f"{alert[SOURCE_NAME_KEY]} {'is' if is_candidate else 'is not'} "
            f"candidate in SkyPortal"
        )

        # check if source exists in SkyPortal
        logger.debug(f"Checking if {alert[SOURCE_NAME_KEY]} is source in SkyPortal")
        response = self.api("HEAD", f"sources/{alert[SOURCE_NAME_KEY]}")

        if response.status_code not in [200, 404]:
            response.raise_for_status()

        is_source = response.status_code == 200
        logger.debug(
            f"{alert[SOURCE_NAME_KEY]} "
            f"{'is' if is_source else 'is not'} source in SkyPortal"
        )

        # object does not exist in SkyPortal: neither cand nor source
        if (not is_candidate) and (not is_source):
            # post candidate
            self.skyportal_post_candidate(alert)

            # post annotations
            self.skyportal_post_annotation(alert)

            # post full light curve
            logger.debug(f"Using stream_id={self.stream_id}")
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
                    f"{alert[SOURCE_NAME_KEY]} from SkyPortal",
                )
                response = self.api(
                    "GET",
                    f"sources/{alert[SOURCE_NAME_KEY]}/groups",
                )
                if response.json()["status"] == "success":
                    existing_groups = response.json()["data"]
                    existing_group_ids = [g["id"] for g in existing_groups]

                    for existing_gid in existing_group_ids:
                        if existing_gid in self.group_ids:
                            self.skyportal_post_source(alert, [existing_gid])
                else:
                    logger.error(
                        f"Failed to get source groups info on {alert[SOURCE_NAME_KEY]}"
                    )
            else:  # exists in SkyPortal but NOT saved as a source
                self.skyportal_post_source(alert)

            # post alert photometry in single call to /api/photometry
            logger.debug(f"Using stream_id={self.stream_id}")
            self.skyportal_put_photometry(alert)

            if self.update_thumbnails:
                self.skyportal_post_thumbnails(alert)

        logger.debug(f"SendToSkyportal Manager complete for {alert[SOURCE_NAME_KEY]}")
