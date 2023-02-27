"""
Module for querying reference images from the UKIRT survey
"""
import logging
from collections.abc import Callable

import astropy.units as u
import numpy as np
from astropy.coordinates import SkyCoord
from astropy.io import fits
from astropy.units import Quantity
from astropy.wcs import WCS
from astroquery.ukidss import UkidssClass
from astroquery.utils.commons import FileContainer
from survey_coverage import get_known_ukirt_surveys
from survey_coverage.surveys import MOCSurvey

from winterdrp.data import Image
from winterdrp.processors.astromatic.swarp.swarp import Swarp
from winterdrp.references.base_reference_generator import BaseReferenceGenerator

logger = logging.getLogger(__name__)

ukirt_image_height = 90 * u.arcmin
ukirt_image_width = 15 * u.arcmin


def get_centre_ra_dec_from_header(header: fits.Header) -> (float, float):
    nx, ny = header["NAXIS1"], header["NAXIS2"]

    wcs = WCS(header)

    ra_cent, dec_cent = wcs.all_pix2world(nx / 2, ny / 2, 0)

    return ra_cent, dec_cent


def get_image_dims_from_header(header: fits.Header) -> (Quantity, Quantity):
    nx, ny = header["NAXIS1"], header["NAXIS2"]
    dx, dy = header["CD1_1"], header["CD2_2"]
    img_dims = (dx * nx * u.deg, dy * ny * u.deg)
    return img_dims


def find_ukirt_surveys(ra: float, dec: float, band: str) -> list[MOCSurvey]:
    surveys = get_known_ukirt_surveys()
    band_surveys = np.array([x for x in surveys if x.filter_name == band])
    in_survey_footprint = [x.contains(ra, dec)[0] for x in band_surveys]
    return band_surveys[in_survey_footprint]


def combine_headers(primary_header, header_to_append):
    for k in header_to_append.keys():
        if k not in primary_header.keys():
            primary_header[k] = header_to_append[k]

    return primary_header


class UKIRTRef(BaseReferenceGenerator):
    """
    UKIRT Ref generator
    """

    def __init__(
        self,
        filter_name: str,
        stack_multiple_images: bool = True,
        swarp_resampler: Callable[..., Swarp] = None,
    ):
        super(UKIRTRef, self).__init__(filter_name)
        self.stack_multiple_images = stack_multiple_images
        self.swarp_resampler = swarp_resampler
        if np.logical_and(self.stack_multiple_images, self.swarp_resampler is None):
            err = (
                "Stacking mutiple reference images requested, but no swarp resampler "
                "suppied"
            )
            raise UKIRTRefError(err)

    def make_image_from_hdulist(
        self, ukirt_hdulist: [fits.hdu.image.PrimaryHDU, fits.hdu.image.ImageHDU]
    ) -> Image:
        assert len(ukirt_hdulist) == 2
        combined_header = combine_headers(
            primary_header=ukirt_hdulist[0].header,
            header_to_append=ukirt_hdulist[1].header,
        )
        data = ukirt_hdulist[1].data
        image = Image(header=combined_header, data=data)
        return image

    def get_reference(self, image: Image) -> fits.PrimaryHDU:
        header = image.get_header()
        ra_cent, dec_cent = get_centre_ra_dec_from_header(header)
        image_x_size, image_y_size = get_image_dims_from_header(header)

        center_crds = SkyCoord(ra=ra_cent, dec=dec_cent, unit=(u.deg, u.deg))
        ukirt_query = UkidssClass()

        # image_width = np.min([15 * u.arcmin, image_x_size.to(u.arcmin)])
        # image_height = np.min([90 * u.arcmin, image_y_size.to(u.arcmin)])
        search_radius = np.max(
            [image_x_size.to(u.arcmin), image_y_size.to(u.arcmin), ukirt_image_height]
        )

        # Sort surveys in descending order of limiting mags
        lim_mags = [x.lim_mag for x in ukirt_surveys]
        ukirt_surveys = ukirt_surveys[np.argsort(lim_mags)[::-1]]

        ukirt_surveys = find_ukirt_surveys(ra_cent, dec_cent, self.filter_name)
        if len(ukirt_surveys) == 0:
            err = "Coordinates not in UKIRT surveys"
            raise NotinUKIRTError(err)

        ukirt_image_urls = []
        for survey in ukirt_surveys:
            ukirt_query.database = survey["database"]
            # url_list = ukirt_query.get_image_list(center_crds,
            #                                       image_width=image_width,
            #                                       image_height=image_height
            #                                       )
            url_list = ukirt_query.get_image_list(center_crds, radius=search_radius)

            ukirt_image_urls += url_list

        if len(ukirt_image_urls) == 0:
            err = "No image found at the given coordinates in the UKIRT database"
            raise UKIRTRefNotFoundError(err)

        if not self.stack_multiple_images:
            # Choose the image from the most sensitive survey,
            # if you do not wish to stack
            ukirt_image_urls = [ukirt_image_urls[0]]

        ukirt_image_objects = [
            FileContainer(
                url,
                encoding="binary",
                remote_timeout=ukirt_query.TIMEOUT,
                show_progress=True,
            )
            for url in ukirt_image_urls
        ]

        ukirt_image_hdulists = [obj.get_fits() for obj in ukirt_image_objects]
        # UKIRT ref images are stored as multiHDU files, need to combine the hdus so
        # no info from the headers is lost.
        ukirt_images = [self.make_image_from_hdulist(x) for x in ukirt_image_hdulists]

        # TODO : Figure out how to stack references (they need to be stacked,
        #  then re-phot-caled)
        if np.logical_and(self.stack_multiple_images, len(ukirt_images) > 1):
            # stack if multiple images are returned
            raise NotImplementedError

        if (len(ukirt_image_urls) == 1) | self.stack_multiple_images:
            pass
