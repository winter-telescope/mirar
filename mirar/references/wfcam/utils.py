import numpy as np
from astropy.io import fits
from astropy.wcs import WCS
from astrosurveyutils import get_known_ukirt_surveys
from astrosurveyutils.surveys import MOCSurvey

from mirar.data import Image
from mirar.data.utils import get_corners_ra_dec_from_header, get_image_center_wcs_coords
from mirar.paths import (
    BASE_NAME_KEY,
    COADD_KEY,
    EXPTIME_KEY,
    GAIN_KEY,
    OBSCLASS_KEY,
    TIME_KEY,
    ZP_KEY,
    ZP_STD_KEY,
    core_fields,
)

MULTIFRAME_ID_KEY = "MFID"
EXTENSION_ID_KEY = "XTNSNID"
LX_KEY = "LX"
LY_KEY = "LY"
HX_KEY = "HX"
HY_KEY = "HY"
COMPID_KEY = "COMPID"

QUERY_RA_KEY = "QRY_RA"
QUERY_DEC_KEY = "QRY_DEC"
QUERY_FILT_KEY = "QRY_FILT"


def get_query_coordinates_from_header(
    header: fits.Header, numpoints: int = 1
) -> (list[float], list[float]):
    """
    Function to get break an image into numpoints sections and get the
    relevsnt coordinates
    Args:
        header:
        numpoints:

    Returns:

    """
    nx, ny = header["NAXIS1"], header["NAXIS2"]

    wcs = WCS(header)

    if numpoints == 1:
        xcrd_list, ycrd_list = [nx / 2], [ny / 2]
    else:
        numpoints = int(np.sqrt(numpoints))
        xcrd_list, ycrd_list = np.linspace(0, nx, numpoints), np.linspace(
            0, ny, numpoints
        )

        crd_list = []
        for i in xcrd_list:
            for j in ycrd_list:
                crd_list.append((i, j))

        xcrd_list = [x[0] for x in crd_list]
        ycrd_list = [x[1] for x in crd_list]

    ra_list, dec_list = wcs.all_pix2world(xcrd_list, ycrd_list, 1)
    ra_list[ra_list < 0] = ra_list[ra_list < 0] + 360
    return ra_list, dec_list


def find_ukirt_surveys(ra: float, dec: float, band: str) -> list[MOCSurvey]:
    """
    Find which UKIRT survey does the given RA/Dec belong to
    Args:
        ra:
        dec:
        band:

    Returns:

    """
    surveys = get_known_ukirt_surveys()
    band_surveys = np.array([x for x in surveys if x.filter_name == band])
    in_survey_footprint = [x.contains(ra, dec)[0] for x in band_surveys]
    return band_surveys[in_survey_footprint]


def combine_headers(primary_header, header_to_append):
    """
    Function to append a header to another
    Args:
        primary_header:
        header_to_append:

    Returns:

    """
    if "SIMPLE" not in primary_header.keys():
        primary_header.insert(0, ("SIMPLE", True))
    if "XTENSION" in primary_header.keys():
        del primary_header["XTENSION"]
    for k in header_to_append.keys():
        if k not in primary_header.keys():
            try:
                primary_header[k] = header_to_append[k]
            except ValueError:
                continue

    return primary_header


def make_wfcam_image_from_hdulist(
    ukirt_hdulist: [fits.hdu.image.PrimaryHDU, fits.hdu.image.ImageHDU],
    url: str,
) -> Image:
    """
    Function to convert a ukirt image with two headers to a single header image
    Args:
        ukirt_hdulist:
        url:
    Returns:

    """
    assert len(ukirt_hdulist) == 2
    # combined_header = ukirt_hdulist[1].header.copy()

    (
        ukirt_filename,
        multiframeid,
        extension_id,
        frame_lx,
        frame_hx,
        frame_ly,
        frame_hy,
    ) = get_wfcam_file_identifiers_from_url(url)
    basename = (
        f"{multiframeid}_{extension_id}_{frame_lx}_"
        f"{frame_hx}_{frame_ly}_{frame_hy}.fits"
    )

    combined_header = combine_headers(
        primary_header=ukirt_hdulist[1].header,
        header_to_append=ukirt_hdulist[0].header,
    )

    combined_header[EXPTIME_KEY] = combined_header["EXP_TIME"]
    combined_header[BASE_NAME_KEY] = basename
    combined_header[GAIN_KEY] = combined_header["GAIN"]
    combined_header[TIME_KEY] = combined_header["DATE-OBS"]
    combined_header[OBSCLASS_KEY] = "ref"
    combined_header[COADD_KEY] = 1
    for key in core_fields:
        if key not in combined_header.keys():
            combined_header[key] = ""
    data = ukirt_hdulist[1].data
    image = Image(header=combined_header, data=data)

    comp_ra_cent, comp_dec_cent = get_image_center_wcs_coords(image, origin=1)
    (
        (ra0_0, dec0_0),
        (ra0_1, dec0_1),
        (ra1_0, dec1_0),
        (ra1_1, dec1_1),
    ) = get_corners_ra_dec_from_header(image.header)
    image.header["RA_CENT"] = comp_ra_cent
    image.header["DEC_CENT"] = comp_dec_cent
    image.header["RA0_0"] = ra0_0
    image.header["DEC0_0"] = dec0_0
    image.header["RA0_1"] = ra0_1
    image.header["DEC0_1"] = dec0_1
    image.header["RA1_0"] = ra1_0
    image.header["DEC1_0"] = dec1_0
    image.header["RA1_1"] = ra1_1
    image.header["DEC1_1"] = dec1_1
    image.header[MULTIFRAME_ID_KEY] = multiframeid
    image.header[EXTENSION_ID_KEY] = extension_id
    image.header[LX_KEY] = frame_lx
    image.header[HX_KEY] = frame_hx
    image.header[LY_KEY] = frame_ly
    image.header[HY_KEY] = frame_hy

    image.header["UKIRPATH"] = ukirt_filename
    image.header[ZP_KEY] = image.header["MAGZPT"]
    image.header[ZP_STD_KEY] = image.header["MAGZRR"]
    image.header[COMPID_KEY] = f"{multiframeid}_{extension_id}"

    if "SEEING" not in image.header.keys():
        image.header["SEEING"] = -99

    return image


def get_wfcam_file_identifiers_from_url(url: str) -> list:
    """
    Function to get the UKIRT file identifiers from the URL
    Args:
        url:
    Returns:

    """
    ukirt_filename = url.split("?")[1].split("&")[0].split("=")[1]
    multiframeid = url.split("&")[1].split("=")[1]
    extension_id = url.split("&")[2].split("=")[1]
    frame_lx = url.split("&")[3].split("=")[1]
    frame_hx = url.split("&")[4].split("=")[1]
    frame_ly = url.split("&")[5].split("=")[1]
    frame_hy = url.split("&")[6].split("=")[1]
    return [
        ukirt_filename,
        int(multiframeid),
        int(extension_id),
        int(frame_lx),
        int(frame_hx),
        int(frame_ly),
        int(frame_hy),
    ]
