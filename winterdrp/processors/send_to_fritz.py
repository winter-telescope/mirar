import logging
import requests
from requests.packages.urllib3.util.retry import Retry
from requests.adapters import HTTPAdapter
import pandas as pd
import os, gzip, io

import base64
from contextlib import contextmanager
import astropy
import json, time
from datetime import datetime
from astropy.time import Time
from astropy.io import ascii
from astropy.io import fits
from astropy.stats import sigma_clipped_stats
from astropy.visualization import (
    AsymmetricPercentileInterval,
    LinearStretch,
    LogStretch,
    ImageNormalize,
)

from copy import deepcopy
import matplotlib.pyplot as plt
import numpy as np
from typing import Mapping, Optional, Sequence


from winterdrp.processors.base_processor import BaseDataframeProcessor
from winterdrp.paths import get_output_dir

logger = logging.getLogger(__name__)

DEFAULT_TIMEOUT = 5  # seconds

class TimeoutHTTPAdapter(HTTPAdapter):
    def __init__(self, *args, **kwargs):
        self.timeout = DEFAULT_TIMEOUT
        if "timeout" in kwargs:
            self.timeout = kwargs["timeout"]
            del kwargs["timeout"]
        super().__init__(*args, **kwargs)

    def send(self, request, **kwargs):
        try:
            timeout = kwargs.get("timeout")
            if timeout is None:
                kwargs["timeout"] = self.timeout
            return super().send(request, **kwargs)
        except AttributeError:
            kwargs["timeout"] = DEFAULT_TIMEOUT

class SendToFritz(BaseDataframeProcessor):
    def __init__(self, 
                token = None,
                group_ids = [1431],
                base_name = 'WIRC',
                filter_id = 1152,
                instrument_id = 5,
                stream_id = 1005,
                protocol = 'http',
                verbose = 2,
                *args,
                **kwargs):
        super(SendToFritz, self).__init__(**kwargs)
        self.token = None
        self.group_ids = group_ids
        self.base_name = base_name
        self.filter_id = filter_id
        self.instrument_id = instrument_id
        self.origin = base_name # used for sending updates to Fritz
        self.stream_id = stream_id
        self.protocol = protocol
        self.verbose = verbose

        def _get_fritz_token():
            token_fritz = os.getenv("FRITZ_TOKEN")

            if token_fritz is None:
                err = "No Fritz token specified. Run 'export FRITZ_TOKEN=<token>' to set. " \
                    "The Fritz token will need to be specified manually for Fritz API queries."
                logger.warning(err)
                raise ValueError
            
            return token_fritz
        
        # session to talk to SkyPortal/Fritz
        self.session = requests.Session()
        self.session_headers = {
            "Authorization": f"token {_get_fritz_token()}",
            "User-Agent": f"winterdrp",
        }

        retries = Retry(
            total=5,
            backoff_factor=2,
            status_forcelist=[405, 429, 500, 502, 503, 504],
            method_whitelist=["HEAD", "GET", "PUT", "POST", "PATCH"],
        )
        adapter = TimeoutHTTPAdapter(timeout=5, max_retries=retries)
        self.session.mount("https://", adapter)
        self.session.mount("http://", adapter)
        self.origin = base_name # used for sending updates to Fritz
        self.verbose = verbose

        def _get_fritz_token():
            token_fritz = os.getenv("FRITZ_TOKEN")

            if token_fritz is None:
                err = "No Fritz token specified. Run 'export FRITZ_TOKEN=<token>' to set. " \
                    "The Fritz token will need to be specified manually for Fritz API queries."
                logger.warning(err)
                raise ValueError
            
            return token_fritz
        
        # session to talk to SkyPortal/Fritz
        self.session = requests.Session()
        self.session_headers = {
            "Authorization": f"token {_get_fritz_token()}",
            "User-Agent": f"winterdrp",
        }

        retries = Retry(
            total=5,
            backoff_factor=2,
            status_forcelist=[405, 429, 500, 502, 503, 504],
            method_whitelist=["HEAD", "GET", "PUT", "POST", "PATCH"],
        )
        adapter = TimeoutHTTPAdapter(timeout=5, max_retries=retries)
        self.session.mount("https://", adapter)
        self.session.mount("http://", adapter)

    def _apply_to_candidates(
            self,
            candidate_table: pd.DataFrame,
    ) -> pd.DataFrame:
        logger.info("In SendToFritz")
        # self.token = self._get_fritz_token()
        self.make_alert(candidate_table) 
        return candidate_table

    def _get_author_id(self):
        """Fritz author id is used in update calls.
        Can be found """
        authid_fritz = os.getenv("FRITZ_AUTHID")

        if authid_fritz is None:
            err = "No Fritz author id specified. Run 'export FRITZ_AUTHID=<id>' to set. " \
                "Author id needs to be specified for updates sent by Fritz API queries."
            logger.warning(err)
            raise ValueError

        return authid_fritz

    def open_bytes_obj(self, bytes_obj):
        """Return numpy array of bytes_obj
        
        Args:
            bytes_obj (_io.BytesIO object in memory): BytesIO obj representing image.

        Returns:
            numpy.ndarrary: representation of the image        
        """
        bytes_io = io.BytesIO(gzip.open(io.BytesIO(bytes_obj.getvalue())).read()) # io.BytesIO obj, ready to be read by fits.open    
        cutout = fits.open(bytes_io)[0].data
        return cutout

    def read_input_df(self, df):
        """Takes a DataFrame, which has multiple candidate 
        and creates list of dictionaries, each dictionary 
        representing a single candidate.

        NOTE: saves images as np.arrays for thumbnail sending.

        Args:
            df (pandas.core.frame.DataFrame): dataframe of all candidates.
        
        Returns:
            (list[dict]): list of dictionaries, each a candidate.
        """
        all_candidates = []   
        
        for i in range(0, len(df)):
            candidate = {} 
            for key in df.keys():
                try: 
                    if type(df.iloc[i].get(key)) is str:
                        candidate[key] = df.iloc[i].get(key)   
                    else:
                        # change to native python type
                        candidate[key] = df.iloc[i].get(key).item()
                except AttributeError: # for IOBytes objs
                    candidate[key] = self.open_bytes_obj(df.iloc[i].get(key))
                                                 
            all_candidates.append(candidate)

        return all_candidates 

    def get_next_name(self, lastname, candjd, bwfile = 'badwords.txt', begcount = 'aaaaaaa'):
        """Creates candidate name following the naming format of 'WNTR22aaaaaaa' .
        Modified from https://github.com/dekishalay/pgirdps

        Args:
            lastname(str): last used candidate name.
            candjd (str): candidate's JD.
            bwfile (str): file name of .txt file of excluded words. 
            begcount (str): string to start naming convention.
        
        Returns:
            (str): next candidate name in sequence.
        """ 
        curyear = Time(candjd, format = 'jd').datetime.strftime('%Y')[2:4]

        if lastname is None:
            #If this is the first source being named
            newname = self.base_name + curyear + begcount
            return newname
        
        lastyear = lastname[4:6]

        if curyear != lastyear:
            #If this is the first candidate of the new year, start with aaaaa
            newname = self.base_name + curyear + begcount
            return newname
        else:
            lastcount = lastname[6:]
            charpos = len(lastcount) - 1
            # will iteratively try to increment characters starting from the last
            inctrue = False
            usestring = ''
            while charpos >= 0:
                cref = lastcount[charpos]
                if inctrue:
                    usestring = cref + usestring
                    charpos -= 1
                    continue
                creford = ord(cref)
                #increment each character, if at 'z', increment the next one
                if creford + 1 > 122:
                    usestring = 'a' + usestring
                else:
                    nextchar = chr(creford+1)
                    usestring = nextchar + usestring
                    inctrue = True
                charpos -= 1
                continue
            
            newname = self.base_name + curyear + usestring

            curdir = os.path.dirname(__file__) # /data/sulekha/winterdrp/winterdrp/processors/alert_packets
            file_path = os.path.join(curdir, bwfile)

            bwlist = ascii.read(file_path, format = 'no_header')
            isbw = False
            # check for bad word
            for i in range(len(bwlist)):
                if usestring.find(str(bwlist['col1'][i])) != -1:
                    #Shame shame
                    isbw = True
                    break
            if isbw:
                # increment the name with a recursive call
                return self.get_next_name(newname, candjd)
            else:
                return newname

    def api_old(self, method, endpoint, data=None):
        """Skyportal API Query"""
        headers = {'Authorization': f'token {self.token}'}
        response = requests.request(method, endpoint, json=data, headers=headers)
        return response

    def api(self, method: str, endpoint: str, data: Optional[Mapping] = None):
        """Make an API call to a SkyPortal instance

        :param method:
        :param endpoint:
        :param data:
        :return:
        """
        method = method.lower()
        methods = {
            "head": self.session.head,
            "get": self.session.get,
            "post": self.session.post,
            "put": self.session.put,
            "patch": self.session.patch,
            "delete": self.session.delete,
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

    # updated, formerly add_new_source
    def alert_post_source(self, cand, group_ids = None):
        if group_ids is None: 
            group_ids = self.group_ids


        data = {'ra': cand['ra'],
                        'dec': cand['dec'],
                        'id': cand['objectId'],
                        'group_ids': group_ids,
                        'origin': self.origin
        }   
        
        logger.info(f"Saving {cand['objectId']} {cand['candid']} as a Source on SkyPortal")
        response = self.api('POST', 'https://fritz.science/api/sources', data)
        
        if response.json()["status"] == "success":
            logger.info(f"Saved {cand['objectId']} {cand['candid']} as a Source on SkyPortal")
        else:
            logger.info(
                f"Failed to save {cand['objectId']} {cand['candid']} as a Source on SkyPortal"
            )
            logger.info(response.json())
    
    # updated 
    def make_thumbnail(
    self, alert, skyportal_type: str, alert_packet_type: str
    ):
        """
        Convert lossless FITS cutouts from ZTF-like alerts into PNGs

        :param alert: ZTF-like alert packet/dict
        :param skyportal_type: <new|ref|sub> thumbnail type expected by SkyPortal
        :param alert_packet_type: <Science|Template|Difference> survey naming
        :return:
        """
        alert = deepcopy(alert)

        cutout_data = alert[f"cutout{alert_packet_type}"]
    
        with gzip.open(io.BytesIO(cutout_data), "rb") as f:
            with fits.open(io.BytesIO(f.read()), ignore_missing_simple=True) as hdu:
                image_data = hdu[0].data

        # TODO  for WNTR
        # Survey-specific transformations to get North up and West on the right
        if self.instrument == "ZTF":
            image_data = np.flipud(image_data)
        elif self.instrument == "PGIR":
            image_data = np.rot90(np.fliplr(image_data), 3)

        buff = io.BytesIO()
        plt.close("all")
        fig = plt.figure()
        fig.set_size_inches(4, 4, forward=False)
        ax = plt.Axes(fig, [0.0, 0.0, 1.0, 1.0])
        ax.set_axis_off()
        fig.add_axes(ax)

        # replace nans with median:
        img = np.array(image_data)
        # replace dubiously large values
        xl = np.greater(np.abs(img), 1e20, where=~np.isnan(img))
        if img[xl].any():
            img[xl] = np.nan
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
        ax.imshow(img_norm, cmap="bone", origin="lower", vmin=vmin, vmax=vmax)
        plt.savefig(buff, dpi=42)

        buff.seek(0)
        plt.close("all")

        thumbnail_dict = {
            "obj_id": alert["objectId"],
            "data": base64.b64encode(buff.read()).decode("utf-8"),
            "ttype": skyportal_type,
        }

        return thumbnail_dict

    #updated thumbnail sending
    def alert_post_thumbnails(self, alert):
        """Post alert thumbnails to SkyPortal

        :param alert:
        :return:
        """
        for ttype, instrument_type in [
            ("new", "Science"),
            ("ref", "Template"),
            ("sub", "Difference"),
        ]:
            logger.info(
                f"Making {instrument_type} thumbnail for {alert['objectId']} {alert['candid']}",
                
            )
            thumb = self.make_thumbnail(alert, ttype, instrument_type)

            logger.info(
                f"Posting {instrument_type} thumbnail for {alert['objectId']} {alert['candid']} to SkyPortal", 
            )
            response = self.api("POST", "https://fritz.science/api/thumbnail", thumb)

            if response.json()["status"] == "success":
                logger.info(
                    f"Posted {alert['objectId']} {alert['candid']} {instrument_type} cutout to SkyPortal"
                )
            else:
                logger.info(
                    f"Failed to post {alert['objectId']} {alert['candid']} {instrument_type} cutout to SkyPortal"
                )
                logger.info(response.json())

    def upload_thumbnail(self, cand):
        """Post new thumbnail to Fritz.

        Format of thumbnail payload:
        { "obj_id": "string",  "data": "string",  "ttype": "string"}
        """
        fritz_to_cand = {"new": 'SciBitIm', "ref": 'RefBitIm', "sub": 'DiffBitIm'}

        for fritz_key in fritz_to_cand.keys():
            cand_key = fritz_to_cand[fritz_key]
            cutout = cand[cand_key]

            buffer = io.BytesIO()
            plt.figure(figsize=(3,3))
            mean, median, std = sigma_clipped_stats(cutout)
            plt.imshow(cutout, origin='lower', cmap='gray',vmin=mean-1*std,vmax=median+3*std)
            plt.xticks([])
            plt.yticks([])

            plt.savefig(buffer,format='png')

            cutoutb64 = base64.b64encode(buffer.getvalue())
            cutoutb64_string = cutoutb64.decode('utf8')

            data_payload = {'obj_id':cand['objectId'],
                            'data':cutoutb64_string,
                            'ttype':fritz_key
                        }

            response = self.api('POST', 'https://fritz.science/api/thumbnail', data=data_payload)
            # logger.info(f'candid {data_payload["obj_id"]}: {data_payload["ttype"]}, thumbnail response:{response}')
        return response

    # updated
    def make_photometry(self, cand, jd_start: float = None):
        """
        Make a de-duplicated pandas.DataFrame with photometry of cand['objectId']
        Modified from Kowalksi (https://github.com/dmitryduev/kowalski)

        :param cand: candidate dictionary
        :param jd_start: date from which to start photometry from
        """
        print(cand.keys())
        cand = deepcopy(cand)
        print(cand.keys())
        # df_candidate = pd.DataFrame(cand["candidate"], index=[0])
        df_candidate = pd.DataFrame(cand, index=[0])

        #TODO fix!!
        # df_prv_candidates = pd.DataFrame(cand["prv_candidates"])
        df_prv_candidates = pd.DataFrame(cand)

        df_light_curve = pd.concat(
            [df_candidate, df_prv_candidates], ignore_index=True, sort=False
        )
        
        # note: WNTR (like PGIR) uses 2massj, which is not in sncosmo as of 
        # 20210803, cspjs seems to be close/good enough as an approximation
        df_light_curve["filter"] = "cspjs"

        df_light_curve["magsys"] = "ab"
        df_light_curve["mjd"] = df_light_curve["jd"] - 2400000.5

        df_light_curve["mjd"] = df_light_curve["mjd"].apply(lambda x: np.float64(x))
        df_light_curve["magpsf"] = df_light_curve["magpsf"].apply(
            lambda x: np.float32(x)
        )
        df_light_curve["sigmapsf"] = df_light_curve["sigmapsf"].apply(
            lambda x: np.float32(x)
        )

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

    #updated to match K
    def alert_put_photometry(self, cand):
        """Send photometry to Fritz."""
        logger.info(f"Making alert photometry of {cand['objectId']} {cand['candid']}")
        df_photometry = self.make_photometry(cand)

        # post photometry
        photometry = {
            "obj_id": cand["objectId"],
            "stream_ids": [int(self.pgir_stream_id)],
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
            logger.info(f"Posting photometry of {cand['objectId']} {cand['candid']}, "
                    f"stream_id={self.stream_id} to SkyPortal")
            response = self.api_skyportal("PUT", "https://fritz.science/api/photometry", photometry)
            if response.json()["status"] == "success":
                logger.info(
                    f"Posted {cand['objectId']} photometry stream_id={self.stream_id} to SkyPortal"
                )
            else:
                logger.info(
                    f"Failed to post {cand['objectId']} photometry stream_id={self.stream_id} to SkyPortal"
                )
                logger.info(response.json())

    # updated: equivalent to create_new_cand:
    def alert_post_candidate(self, cand):
        """
        Post a candidate on SkyPortal. Creates new candidate(s) (one per filter)
        """
        data = { "id": cand["objectId"],
                "ra": cand["ra"],
                "dec": cand["dec"],
                "filter_ids": [self.filter_id],
                "passing_alert_id": self.filter_id,
                "passed_at": Time(datetime.utcnow()).isot,
                "origin": "WINTERdrp",
                }
        
        logger.info(f"Posting metadata of {cand['objectId']} {cand['candid']} to SkyPortal")
        response = self.api("POST", "https://fritz.science/api/candidates", data)

        if response.json()["status"] == "success":
            logger.info(f"Posted {cand['objectId']} {cand['candid']} metadata to SkyPortal")
        else:
            logger.info(f"Failed to post {cand['objectId']} {cand['candid']} metadata to SkyPortal")
            logger.info(response.json())

    # updated
    def alert_check_cand(self, cand):
        """Checks whether a candidate already exist."""
        response = self.api('HEAD', f'https://fritz.science/api/candidates/{str(cand["objectId"])}')
        return response        
    
    # updated to match K
    def alert_post_annotation(self, cand):
        """Post an annotation. Works for both candidates and sources.
        """
        data = {"chipsf": cand["chipsf"],
                "fwhm": cand["fwhm"],
                "scorr": cand["scorr"]}
        payload = {"origin": self.origin,
            "data": data,
            "group_ids": self.group_ids
            }

        path = f'https://fritz.science/api/sources/{str(cand["objectId"])}/annotations'
        response = self.api('POST', path, payload)

        if response.json()["status"] == "success":
                logger.info(f"Posted {cand['objectId']} annotation to SkyPortal")
        else:
            logger.info(f"Failed to post {cand['objectId']} annotation to SkyPortal")
            logger.info(response.json())

    #updated to match alert_put_annotations, depreciates update_annotation & retrieve_annotation_specified_source
    def alert_put_annotation(self, cand):
        """Retrieve an annotation to check if it exists already."""
        response = self.api('GET', f'https://fritz.science/api/sources/{str(cand["objectId"])}/annotations')

        if response.json()["status"] == "success":
            logger.info(f"Got {cand['objectId']} annotations from SkyPortal")
        else:
            logger.info(f"Failed to get {cand['objectId']} annotations from SkyPortal")
            logger.info(response.json())
            return False

        existing_annotations = {
            annotation["origin"]: {
                "annotation_id": annotation["id"],
                "author_id": annotation["author_id"],
            }
            for annotation in response.json()["data"]
        }

        # no annotation exists on SkyPortal for this object? just post then
        if self.origin not in existing_annotations:
            self.alert_post_annotation(cand)
        # annotation from this(WNTR) origin exists
        else:
            # annotation data
            data = {"fwhm":cand["fwhm"],
                "scorr": cand["scorr"],
                "chipsf": cand["chipsf"]
            }
            new_annotation = {
                "author_id": existing_annotations[self.origin]["author_id"],
                "obj_id": cand["objectId"],
                "origin": self.origin,
                "data": data,
                "group_ids": self.group_ids,
            }
            
            logger.info(
                f"Putting annotation for {cand['objectId']} {cand['candid']} to SkyPortal", 
            )
            response = self.api(
                    "PUT",
                    f"https://fritz.science/api/sources/{cand['objectId']}"
                    f"/annotations/{existing_annotations[self.origin]['annotation_id']}",
                    new_annotation,
                )
            if response.json()["status"] == "success":
                logger.info(f"Posted updated {cand['objectId']} annotation to SkyPortal")
            else:
                logger.info(f"Failed to post updated {cand['objectId']} annotation to SkyPortal")
                logger.info(response.json())

    def alert_skyportal_manager(self, alert):
        """Posts alerts to SkyPortal if criteria is met

        :param alert: _description_
        :type alert: _type_
        """
        # check if candidate exists in SkyPortal
        logger.info(f"Checking if {alert['objectId']} is candidate in SkyPortal")
        response = self.api("HEAD", f"https://fritz.science/api/candidates/{alert['objectId']}")
        is_candidate = response.status_code == 200
        logger.info(f"{alert['objectId']} {'is' if is_candidate else 'is not'} candidate in SkyPortal")
        
        # check if source exists in SkyPortal
        logger.info(f"Checking if {alert['objectId']} is source in SkyPortal")
        response = self.api("HEAD", f"https://fritz.science/api/sources/{alert['objectId']}")
        is_source = response.status_code == 200
        logger.info(f"{alert['objectId']} {'is' if is_source else 'is not'} source in SkyPortal")

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
                logger.info(
                    f"Getting source groups info on {alert['objectId']} from SkyPortal",
                )
                response = self.api(
                        "GET", f"https://fritz.science/api/sources/{alert['objectId']}/groups"
                    )
                if response.json()["status"] == "success":
                    existing_groups = response.json()["data"]
                    existing_group_ids = [g["id"] for g in existing_groups]

                    for existing_gid in existing_group_ids:
                        if existing_gid in self.group_ids:
                            self.alert_post_source(alert, str(existing_gid))    
                else:
                    logger.info(f"Failed to get source groups info on {alert['objectId']}")
            else: # exists in SkyPortal but NOT saved as a source
                self.alert_post_source(alert)

            # post alert photometry in single call to /api/photometry
            prv_cand_exists = "prv_candidates" in alert.keys()
            logger.info(f'!!!!!PREV CANDIDATES IN {alert["objectId"]}: {prv_cand_exists}!!!!!')
            # alert["prv_candidates"] = prv_candidates

            # TODO
            # self.alert_put_photometry(alert)

        logger.info(f'======== Manager complete for {alert["objectId"]} =======')

    def alert_skyportal_manager(self, alert):
        """Posts alerts to SkyPortal if criteria is met

        :param alert: _description_
        :type alert: _type_
        """
        # check if candidate exists in SkyPortal
        logger.info(f"Checking if {alert['objectId']} is candidate in SkyPortal")
        response = self.api("HEAD", f"https://fritz.science/api/candidates/{alert['objectId']}")
        is_candidate = response.status_code == 200
        logger.info(f"{alert['objectId']} {'is' if is_candidate else 'is not'} candidate in SkyPortal")
        
        # check if source exists in SkyPortal
        logger.info(f"Checking if {alert['objectId']} is source in SkyPortal")
        response = self.api("HEAD", f"https://fritz.science/api/sources/{alert['objectId']}")
        is_source = response.status_code == 200
        logger.info(f"{alert['objectId']} {'is' if is_source else 'is not'} source in SkyPortal")

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
                logger.info(
                    f"Getting source groups info on {alert['objectId']} from SkyPortal",
                )
                response = self.api(
                        "GET", f"https://fritz.science/api/sources/{alert['objectId']}/groups"
                    )
                if response.json()["status"] == "success":
                    existing_groups = response.json()["data"]
                    existing_group_ids = [g["id"] for g in existing_groups]

                    for existing_gid in existing_group_ids:
                        if existing_gid in self.group_ids:
                            self.alert_post_source(alert, str(existing_gid))    
                else:
                    logger.info(f"Failed to get source groups info on {alert['objectId']}")
            else: # exists in SkyPortal but NOT saved as a source
                self.alert_post_source(alert)

            # post alert photometry in single call to /api/photometry
            prv_cand_exists = "prv_candidates" in alert.keys()
            logger.info(f'!!!!!PREV CANDIDATES IN {alert["objectId"]}: {prv_cand_exists}!!!!!')
            # alert["prv_candidates"] = prv_candidates

            # TODO
            # self.alert_put_photometry(alert)

        logger.info(f'======== Manager complete for {alert["objectId"]} =======')

    def make_alert(self, cand_table):
        t0 = time.time()
        all_cands = self.read_input_df(cand_table)
        num_cands = len(all_cands)

        last_name = 'WIRC21aaaaaao'
        cand_id = 700
              
        for cand in all_cands:
            cand_jd = cand['jd'] # float
            cand_name = self.get_next_name(last_name, str(cand_jd))

            #TODO candid should be coming from naming database
            cand['candid'] = cand_id
            cand_id += 1
            
            # TODO check if cand is new? 
            # cand_name, new_status = self.check_and_insert_source(cand_name, cand)

            #TODO remove once new_status is up-to-date: dummy line
            new_status = True

            if new_status:
                last_name = cand_name    
            cand['objectId'] = cand_name

            # Old: 
            # source_response = self.add_new_source(cand)
            # source_response = self.create_new_cand(cand, id)
            # logger.info(f'Add source {cand["objectId"]}: {source_response}')
            # thumbnail_response = self.upload_thumbnail(cand)
            # logger.info(f'Upload thumbnail {cand["objectId"]}: {thumbnail_response}')
            # photometry_response = self.update_photometry(cand)
            # logger.info(f'Photometry {cand["objectId"]}: {photometry_response}')
            # self.retrieve_annotation_specified_source(cand)
            # annotation_response = self.post_annotation(cand)

            self.alert_skyportal_manager(cand)
            

        t1 = time.time()
        logger.info('###########################################')
        logger.info(f"Took {(t1 - t0):.2f} seconds to Fritz process {num_cands} candidates.")
        logger.info('###########################################')    