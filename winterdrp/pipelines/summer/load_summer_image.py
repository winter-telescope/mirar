import os
import numpy as np
import pkg_resources
import astropy
import logging
from astropy.coordinates import SkyCoord
from astropy.io import fits
from astropy.time import Time
from astropy import units as u

from winterdrp.paths import latest_save_key, raw_img_key, base_name_key, proc_history_key, proc_fail_key

logger = logging.getLogger(__name__)


def load_raw_summer_image(
        path: str
) -> tuple[np.array, astropy.io.fits.Header]:
    with fits.open(path) as data:
        header = data[0].header
        header["OBSCLASS"] = ["calibration", "science"][header["OBSTYPE"] == "SCIENCE"]
        header['UTCTIME'] = header['UTCSHUT']
        header['TARGET'] = header['OBSTYPE'].lower()

        crd = SkyCoord(ra=data[0].header['RA'], dec=data[0].header['DEC'], unit=(u.deg, u.deg))
        header['RA'] = crd.ra.deg
        header['DEC'] = crd.dec.deg

        header['CRVAL1'] = header['RA']
        header['CRVAL2'] = header['DEC']

        tel_crd = SkyCoord(ra=data[0].header['TELRA'], dec=data[0].header['TELDEC'], unit=(u.deg, u.deg))
        header['TELRA'] = tel_crd.ra.deg
        header['TELDEC'] = tel_crd.dec.deg
        header['BZERO'] = 0

        header[latest_save_key] = path
        header[raw_img_key] = path

        data[0].data = data[0].data * 1.0

        if 'other' in header['FILTERID']:
            header['FILTERID'] = 'r'

        header[proc_history_key] = ""

        base_name = os.path.basename(path)
        header[base_name_key] = base_name
        header["EXPID"] = int("".join(base_name.split("_")[1:3])[2:])
        # header["EXPID"] = str(header["NIGHT"]) + str(header["OBSHISTID"])

        pipeline_version = pkg_resources.require("winterdrp")[0].version
        pipeline_version_padded_str = "".join([x.rjust(2, "0") for x in pipeline_version.split(".")])
        header["PROCID"] = int(str(header["EXPID"]) + str(pipeline_version_padded_str))
        # header.append(('GAIN', summer_gain, 'Gain in electrons / ADU'), end=True)

        header['OBSDATE'] = int(header['UTC'].split('_')[0])

        obstime = Time(header['UTCISO'], format='iso')
        t0 = Time('2018-01-01', format='iso')
        header['NIGHT'] = int(obstime.jd) - int(t0.jd)
        header['EXPMJD'] = header['OBSMJD']

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
            header['SUBPROG'] = 'none'

        header['FILTER'] = header['FILTERID']
        try:
            header['SHUTOPEN'] = Time(header['SHUTOPEN'], format='iso').jd
        except ValueError:
            logger.warning(f"Error parsing 'SHUTOPEN' of {path}: ({header['SHUTOPEN']})")

        try:
            header['SHUTCLSD'] = Time(header['SHUTCLSD'], format='iso').jd
        except ValueError:
            logger.warning(f"Error parsing 'SHUTCLSD' of {path}: ({header['SHUTCLSD']})")

        header['PROCFLAG'] = 0

        header[proc_fail_key] = ""

        sunmoon_keywords = [
            'MOONRA', 'MOONDEC', 'MOONILLF', 'MOONPHAS', 'MOONALT', 'SUNAZ', 'SUNALT',
        ]
        for key in sunmoon_keywords:
            val = 0
            if key in header.keys():
                if header[key] not in ['']:
                    val = header[key]
            header[key] = val

        itid_dict = {
            'SCIENCE': 1,
            'BIAS': 2,
            'FLAT': 2,
            'DARK': 2,
            'FOCUS': 3,
            'POINTING': 4,
            'OTHER': 5
        }

        if not header['OBSTYPE'] in itid_dict.keys():
            header['ITID'] = 5
        else:
            header['ITID'] = itid_dict[header['OBSTYPE']]

        if header['FIELDID'] == 'radec':
            header['FIELDID'] = 999999999

        if header['ITID'] != 1:
            header['FIELDID'] = -99

        if 'COADDS' not in header.keys():
            header['COADDS'] = 1

        if header['PROGID'] in ['', 'WINTER']:
            header['PROGID'] = 2

        try:
            header['PROGID'] = int(header['PROGID'])
        except ValueError:
            try:
                progpi = header['PROGID']
                header['PROGID'] = int(header['PROGPI'])
                header['PROGPI'] = progpi
            except KeyError:
                header['PROGID'] = 0

        crds = SkyCoord(ra=header['RA'], dec=header['DEC'], unit=(u.deg, u.deg))
        header['RA'] = crds.ra.deg
        header['DEC'] = crds.dec.deg

        data[0].header = header
    return data[0].data, data[0].header


def load_proc_summer_image(
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
    pipeline_version = pkg_resources.require("winterdrp")[0].version
    pipeline_version_padded_str = "".join([x.rjust(2, "0") for x in pipeline_version.split(".")])
    header['DIFFID'] = int(str(header["EXPID"])+str(pipeline_version_padded_str))
    data[data == 0] = np.nan
    # logger.info(header['CRVAL2'])
    return data, header
