import os
import numpy as np
import astropy
import logging
from astropy.coordinates import SkyCoord
from astropy.io import fits
from astropy.time import Time
from astropy import units as u

from winterdrp.paths import LATEST_SAVE_KEY, RAW_IMG_KEY, BASE_NAME_KEY, PROC_HISTORY_KEY, PROC_FAIL_KEY, __version__

logger = logging.getLogger(__name__)

def load_raw_sedmv2_image(path: str) -> tuple[np.array, astropy.io.fits.Header]:

    with fits.open(path) as data:
        header = data[0].header


        # science / flat / bias / etc...
        if header['IMGTYPE'] == 'object':
            header['OBSTYPE'] = 'SCIENCE'
        else:
            header['OBSTYPE'] = header['IMGTYPE'].upper()
        header['TARGET'] = header['OBSTYPE'].lower()

        header['OBSCLASS'] = ['calibration', 'science'][header['OBSTYPE'] == 'SCIENCE']

        # coordinates
        ## change to deg units
        crd = SkyCoord(ra=header['RA'], dec=header['DEC'], unit=(u.hourangle, u.deg)) # HOURANGLE!!!!!
        # TODO: access RAD instead
        header['RA'] = crd.ra.deg
        header['DEC'] = crd.dec.deg

        tel_crd = SkyCoord(ra=header["TELRA"], dec=header["TELDEC"], unit=(u.hourangle, u.deg)) # HOURANGLE!
        header["TELRA"] = tel_crd.ra.deg
        header["TELDEC"] = tel_crd.dec.deg

        # ASSUME THAT CRVAL1 and CRVAL2 are correct based on a-net solution!!!!!!!!
        #header['CRVAL1'] = header['RA']
        #header['CRVAL2'] = header['DEC']


        data[0].data = data[0].data * 1.0
        data_save = data[0].data
        #except:
        #    data[0].data = data[1].data * 1.0
        #    data_save = data[0].data

        # filter
        header['FILTERID'] = header['FILTER'].split(' ')[1][0]  # overwrites numerical filterid
        header['FILTER'] = header['FILTERID']

        # keys...
        header[LATEST_SAVE_KEY] = path
        header[RAW_IMG_KEY] = path
        header[PROC_HISTORY_KEY] = ""
        header["PROCFLAG"] = 0
        header[PROC_FAIL_KEY] = ""

        base_name = os.path.basename(path)
        header[BASE_NAME_KEY] = base_name
        header["EXPID"] = int("".join(base_name.split("_")[1:3])[2:])

        pipeline_version = __version__
        pipeline_version_padded_str = "".join(
            [x.rjust(2, "0") for x in pipeline_version.split(".")]
        )
        header["PROCID"] = int(str(header["EXPID"]) + str(pipeline_version_padded_str))

        header.append(('GAIN', 1.0, 'Gain in electrons / ADU'), end=True)

        # times
        header['UTCTIME'] = header['UTC']
        header['TIMEUTC'] = header['UTCTIME']
        header['OBSDATE'] = int(header['UTC'].split('_')[0])

        obstime = Time(header['DATE'], format='isot')
        t0 = Time('2018-01-01', format='iso')  # TODO: change this to sedmv2 start date? otherwise NIGHT starts in 1700s
        header['NIGHT'] = int(obstime.jd) - int(t0.jd)  # integer value, night 1, night 2...
        header["EXPMJD"] = header["OBSDATE"] # TODO: verify
        header["SHUTOPEN"] = obstime.jd
        header["SHUTCLSD"] = obstime.jd

        # IDs # TODO: clean this up...
        default_id = 0

        for key in ["PROGID", "OBSID"]:
            if key not in header.keys():
                # logger.warning(f"No {key} found in header of {path}")
                header[key] = default_id
            else:
                try:
                    header[key] = int(header[key])
                except ValueError:
                    header[key] = default_id

        if "SUBPROG" not in header.keys():
            # logger.warning(f"No SUBPROG found in header of {path}")
            header["SUBPROG"] = "none"

        if 'OBSID' not in header.keys():
            logger.warning(f"No {key} found in header of {path}")
            header['OBSID'] = default_id
        else:
            try:
                header['OBSID'] = int(header['OBSID'])
            except ValueError:
                header['OBSID'] = default_id

        header['PROGID'] = 'SEDMv2'
        if header['PROGID'] == 'SEDMv2':
            header['PROGID'] = 3
        try:
            header['PROGID'] = int(header['PROGID'])
        except ValueError:
            try:
                progpi = header['PROGID']
                header['PROGID'] = int(header['PROGPI'])
                header['PROGPI'] = progpi
            except KeyError:
                header['PROGID'] = 0

        itid_dict = {
            'SCIENCE': 1,
            'BIAS': 2,
            'FLAT': 2,
            'DARK': 2,
            'FOCUS': 3,
            'POINTING': 4,
            'OTHER': 5
        }  # may be unnecessary

        if not header['OBSTYPE'] in itid_dict.keys():
            header['ITID'] = 5
        else:
            header['ITID'] = itid_dict[header['OBSTYPE']]

        header['FIELDID'] = 999999999 # TODO: get rid of FIELDID dependence in other files - is that okay?

        # others
        if 'COADDS' not in header.keys():
            header['COADDS'] = 1
        header["BZERO"] = 0

        sunmoon_keywords = [
            "MOONRA",
            "MOONDEC",
            "MOONILLF",
            "MOONPHAS",
            "MOONALT",
            "SUNAZ",
            "SUNALT",
        ]

        for key in sunmoon_keywords:
            val = 0
            if key in header.keys():
                if header[key] not in [""]:
                    val = header[key]
            header[key] = val


    return data_save, data[0].header


def load_proc_sedmv2_image(
        path: str
) -> tuple[np.array, astropy.io.fits.Header]:
    img = fits.open(path)
    data = img[0].data
    header = img[0].header
    if 'ZP' not in header.keys():
        header['ZP'] = header['ZP_AUTO']
        header['ZP_std'] = header['ZP_AUTO_std']
    header['CENTRA'] = header['CRVAL1']
    header['CENTDEC'] = header['CRVAL2']
    pipeline_version = __version__
    pipeline_version_padded_str = "".join([x.rjust(2, "0") for x in pipeline_version.split(".")])
    header['DIFFID'] = int(str(header["EXPID"])+str(pipeline_version_padded_str))
    data[data == 0] = np.nan
    # logger.info(header['CRVAL2'])
    return data, header
