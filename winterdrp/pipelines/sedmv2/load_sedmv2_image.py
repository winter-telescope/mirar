import os
import numpy as np
import astropy
import logging
from astropy.coordinates import SkyCoord
from astropy.io import fits
from astropy.time import Time
from astropy import units as u

from winterdrp.paths import latest_save_key, raw_img_key, base_name_key, proc_history_key, proc_fail_key, __version__

logger = logging.getLogger(__name__)

def load_raw_sedmv2_image(
        path: str
) -> tuple[np.array, astropy.io.fits.Header]:
    with fits.open(path) as data:
        header, header_1 = data[0].header, data[1].header

        # science / flat / bias / etc...
        if header['IMGTYPE'] == 'object':
            header['OBSTYPE'] = 'SCIENCE'
        else:
            header['OBSTYPE'] = header['IMGTYPE'].upper()
        header['TARGET'] = header['OBSTYPE'].lower()

        header['OBSCLASS'] = ['calibration', 'science'][header['OBSTYPE'] == 'SCIENCE']

        # coordinates
        ## change to deg units
        crd = SkyCoord(ra=header['RA'], dec=header['DEC'], unit=(u.deg, u.deg))
        header['RA'] = crd.ra.deg
        header['DEC'] = crd.dec.deg

        header['CRVAL1'] = header['RA']  # TODO: use header_1['CRVAL1'] instead??
        header['CRVAL2'] = header['DEC']

        data[1].data = data[1].data * 1.0

        # filter
        header['FILTERID'] = header['FILTER'].split(' ')[1][0]  # overwrites numerical filterid
        header['FILTER'] = header['FILTERID']

        # keys...
        header[latest_save_key] = path
        header[raw_img_key] = path
        header[proc_history_key] = ""
        header[proc_fail_key] = ""

        base_name = os.path.basename(path)
        header[base_name_key] = base_name

        # header.append(('GAIN', sedmv2_gain, 'Gain in electrons / ADU'), end=True)  # TODO: is this gain true?

        # times
        header['OBSDATE'] = int(header_1['UTC'].split('_')[0])

        obstime = Time(header['DATE'], format='isot')
        t0 = Time('2018-01-01', format='iso')  # TODO: change this to sedmv2 start date? otherwise NIGHT starts in 1700s
        header['NIGHT'] = int(obstime.jd) - int(t0.jd)  # integer value, night 1, night 2...

        # IDs # TODO: clean this up...
        default_id = 0
        if 'OBSID' not in header.keys():
            # logger.warning(f"No {key} found in header of {path}") # TODO: uncomment
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
        }  # may be unecessary

        if not header['OBSTYPE'] in itid_dict.keys():
            header['ITID'] = 5
        else:
            header['ITID'] = itid_dict[header['OBSTYPE']]

        # if header['FIELDID'] == 'radec':
        #    header['FIELDID'] = 999999999 # TODO: get rid of FIELDID dependence in other files - is that okay?

        if 'COADDS' not in header.keys():
            header['COADDS'] = 1

    return data[1].data, data[0].header


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
