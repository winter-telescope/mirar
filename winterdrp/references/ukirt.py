"""
Module for querying reference images from the UKIRT survey
"""
import logging
from collections.abc import Callable
from typing import Type

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

from winterdrp.data import Image, ImageBatch
from winterdrp.errors import ProcessorError
from winterdrp.paths import BASE_NAME_KEY, COADD_KEY, core_fields, get_output_path
from winterdrp.processors.astromatic.swarp.swarp import Swarp
from winterdrp.processors.base_processor import ImageHandler
from winterdrp.processors.photcal import PhotCalibrator
from winterdrp.processors.sqldatabase.basemodel import BaseDB
from winterdrp.references.base_reference_generator import BaseReferenceGenerator

logger = logging.getLogger(__name__)

ukirt_image_height = 90 * u.arcmin
ukirt_image_width = 15 * u.arcmin


class UKIRTRefError(ProcessorError):
    """
    Base UKIRTRef error
    """


class NotinUKIRTError(ProcessorError):
    """
    Error when the coordinates are not in UKIRT footprint
    """


class UKIRTRefNotFoundError(ProcessorError):
    """
    Error when UKIRT ref is not found for some reason
    """


def get_query_coordinates_from_header(
    header: fits.Header, numpoints=1
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

    ra_list, dec_list = wcs.all_pix2world(xcrd_list, ycrd_list, 0)

    return ra_list, dec_list


def get_centre_ra_dec_from_header(header: fits.Header) -> (float, float):
    """
    Function to get center RA/Dec of the image from the header
    Args:
        header:

    Returns:

    """
    nx, ny = header["NAXIS1"], header["NAXIS2"]

    wcs = WCS(header)

    ra_cent, dec_cent = wcs.all_pix2world(nx / 2, ny / 2, 0)

    return ra_cent, dec_cent


def get_image_dims_from_header(header: fits.Header) -> (Quantity, Quantity):
    """
    Get image dimensions from the header
    Args:
        header:

    Returns:

    """
    nx, ny = header["NAXIS1"], header["NAXIS2"]
    dx, dy = header["CD1_1"], header["CD2_2"]
    img_dims = (dx * nx * u.deg, dy * ny * u.deg)
    return img_dims


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


def make_image_from_hdulist(
    ukirt_hdulist: [fits.hdu.image.PrimaryHDU, fits.hdu.image.ImageHDU], basename
) -> Image:
    """
    Function to convert a ukirt image with two headers to a single header image
    Args:
        ukirt_hdulist:
        basename:

    Returns:

    """
    assert len(ukirt_hdulist) == 2
    # combined_header = ukirt_hdulist[1].header.copy()

    combined_header = combine_headers(
        primary_header=ukirt_hdulist[1].header,
        header_to_append=ukirt_hdulist[0].header,
    )
    combined_header["EXPTIME"] = combined_header["EXP_TIME"]
    combined_header[BASE_NAME_KEY] = basename
    combined_header[COADD_KEY] = 1
    for key in core_fields:
        if key not in combined_header.keys():
            combined_header[key] = ""
    data = ukirt_hdulist[1].data
    image = Image(header=combined_header, data=data)
    return image


class UKIRTRef(BaseReferenceGenerator, ImageHandler):
    """
    UKIRT Ref generator
    """

    abbreviation = "ukirt_ref_download"

    def __init__(
        self,
        filter_name: str,
        num_query_points: int = 4,
        stack_multiple_images: bool = True,
        swarp_resampler: Callable[..., Swarp] = None,
        phot_calibrator: Callable[..., PhotCalibrator] = None,
        check_local_database: bool = True,
        db_table: Type[BaseDB] = None,
    ):
        super().__init__(filter_name)
        self.stack_multiple_images = stack_multiple_images
        self.swarp_resampler = swarp_resampler
        self.num_query_points = num_query_points
        self.phot_calibrator = phot_calibrator
        if np.logical_and(self.stack_multiple_images, self.swarp_resampler is None):
            err = (
                "Stacking mutiple reference images requested, but no swarp resampler "
                "suppied"
            )
            raise UKIRTRefError(err)

        self.check_local_database = check_local_database
        self.db_table = db_table
        if np.logical_and(self.check_local_database, self.db_table is None):
            err = (
                "You have requested checking locally, but no database "
                "has been specified"
            )
            raise UKIRTRefError(err)

    def check_ref_images_exist_locally(self, query_ra, query_dec, db_table):
        """
        Function to check if component images exist locally
        Args:
            query_ra:
            query_dec:
            db_table:

        Returns:

        """
        # table = db_table(query_ra=query_ra,
        #                  query_dec=query_dec,
        #                  savepath=' ',
        #                  query_url=' ')
        # exists = table.exists()
        results = db_table.sql_model().select_query(
            select_keys="savepath",
            compare_values=[query_ra, query_dec],
            compare_keys=["query_ra", "query_dec"],
            comparators=["__eq__", "__eq__"],
        )
        logger.debug(results)
        if len(results) == 0:
            savepaths = []
        else:
            savepaths = [x[0] for x in results]
        return savepaths

    def get_reference(self, image: Image) -> fits.PrimaryHDU:
        header = image.get_header()
        # ra_cent, dec_cent = get_centre_ra_dec_from_header(header)
        # image_x_size, image_y_size = get_image_dims_from_header(header)

        query_ra_list, query_dec_list = get_query_coordinates_from_header(
            header, numpoints=self.num_query_points
        )

        query_crds = SkyCoord(ra=query_ra_list, dec=query_dec_list, unit=(u.deg, u.deg))
        logger.debug(f"Querying around {query_crds}")

        ukirt_query = UkidssClass()

        # image_width = np.min([15 * u.arcmin, image_x_size.to(u.arcmin)])
        # image_height = np.min([90 * u.arcmin, image_y_size.to(u.arcmin)])
        # search_radius = np.max(
        #     [image_x_size.to(u.arcmin), image_y_size.to(u.arcmin), ukirt_image_height]
        # )

        ra_cent, dec_cent = get_centre_ra_dec_from_header(header)
        logger.debug(f"Center RA: {ra_cent} Dec: {dec_cent}")
        ukirt_surveys = find_ukirt_surveys(ra_cent, dec_cent, self.filter_name)
        if len(ukirt_surveys) == 0:
            err = "Coordinates not in UKIRT surveys"
            raise NotinUKIRTError(err)

        # Sort surveys in descending order of limiting mags
        lim_mags = [x.lim_mag for x in ukirt_surveys]
        ukirt_surveys = ukirt_surveys[np.argsort(lim_mags)[::-1]]
        logger.debug(f"Surveys are {[x.survey_name for x in ukirt_surveys]}")
        ukirt_image_urls, ukirt_query_ras, ukirt_query_decs = [], [], []
        for survey in ukirt_surveys:
            ukirt_query.database = survey.wfau_dbname
            # url_list = ukirt_query.get_image_list(query_crds,
            #                                       image_width=image_width,
            #                                       image_height=image_height
            #                                       )
            # This seems to pull multiextension images, which we don't really want.
            # url_list = ukirt_query.get_image_list(center_crds, radius=search_radius)
            url_list, ukirt_qra_list, ukirt_qdec_list = [], [], []
            for crd in query_crds:
                # Need to add a cache and check there.
                url = []
                if self.check_local_database:
                    url = self.check_ref_images_exist_locally(
                        query_ra=crd.ra.deg,
                        query_dec=crd.dec.deg,
                        db_table=self.db_table,
                    )

                if len(url) == 0:
                    url = ukirt_query.get_image_list(
                        crd,
                        image_width=ukirt_image_width,
                        image_height=ukirt_image_height,
                    )
                ra_list, dec_list = [crd.ra.deg] * len(url), [crd.dec.deg] * len(url)
                url_list += url
                ukirt_qra_list += ra_list
                ukirt_qdec_list += dec_list

            ukirt_image_urls += url_list
            ukirt_query_ras += ukirt_qra_list
            ukirt_query_decs += ukirt_qdec_list

        logger.debug(f"UKIRT image url list {ukirt_image_urls}")

        if len(ukirt_image_urls) == 0:
            err = "No image found at the given coordinates in the UKIRT database"
            raise UKIRTRefNotFoundError(err)

        if not self.stack_multiple_images:
            # Choose the image from the most sensitive survey,
            # if you do not wish to stack
            ukirt_image_urls = [ukirt_image_urls[0]]

        ukirt_images = []
        for url, ra, dec in zip(ukirt_image_urls, ukirt_query_ras, ukirt_query_decs):
            if "http" in url:
                obj = FileContainer(
                    url,
                    encoding="binary",
                    remote_timeout=ukirt_query.TIMEOUT,
                    show_progress=True,
                )
                ukirt_img_hdulist = obj.get_fits()
                # UKIRT ref images are stored as multiHDU files, need to combine the
                # hdus so no info from the headers is lost.
                ukirt_image = make_image_from_hdulist(
                    ukirt_img_hdulist,
                    basename=f"{url.split('/')[-1].split('&')[0].split('.')[0]}"
                    f".fits",
                )
                savepath = get_output_path(
                    ukirt_image[BASE_NAME_KEY],
                    dir_root="mock/ref",
                    sub_dir="ir_refbuild",
                )

                if self.check_local_database:
                    new = self.db_table(
                        query_ra=ra,
                        query_dec=dec,
                        savepath=savepath.as_posix(),
                        query_url=url,
                    )
                    ret = new.insert_entry()
                ukirt_image[ret[0][0]] = ret[1][0]
                self.save_fits(ukirt_image, savepath)

            else:
                ukirt_hdulist = fits.open(url, ignore_missing_simple=True)
                ukirt_image = Image(
                    header=ukirt_hdulist[0].header, data=ukirt_hdulist[0].data
                )

            ukirt_images.append(ukirt_image)

        ukirt_images = np.array(ukirt_images)
        # TODO : Figure out how to stack references (they need to be stacked,
        mag_zps = np.array(
            [
                x["MAGZPT"]
                + 2.5 * np.log10(x["EXPTIME"])
                - x["EXTINCT"] * ((x["AMSTART"] + x["AMEND"]) / 2)
                for x in ukirt_images
            ]
        )
        magerr_zps = np.array([x["MAGZRR"] for x in ukirt_images])
        median_mag_zp = np.median(mag_zps)

        scaling_factors = 10 ** (0.4 * (median_mag_zp - mag_zps))
        logger.debug(mag_zps)
        logger.debug(magerr_zps)
        logger.debug(median_mag_zp)
        logger.debug(scaling_factors)
        zpmask = np.abs(mag_zps - median_mag_zp) < 0.5

        ukirt_images = ukirt_images[zpmask]
        scaling_factors = scaling_factors[zpmask]
        for ind, image in enumerate(ukirt_images):
            image["FLXSCALE"] = scaling_factors[ind]

        ukirt_image_batch = ImageBatch(list(ukirt_images))
        resampler = self.swarp_resampler(
            center_ra=ra_cent,
            center_dec=dec_cent,
            include_scamp=False,
            combine=True,
            x_imgpixsize=2270,
            y_imgpixsize=2270,
        )
        resampler.set_night(night_sub_dir="ir_refbuild/mock")
        resampled_batch = resampler.apply(ukirt_image_batch)

        resampled_image = resampled_batch[0]
        # phot_calibrator = self.phot_calibrator(resampled_image)
        # phot_calibrator.set_night(night_sub_dir="ir_refbuild/mock")
        # photcaled_batch = phot_calibrator.apply(resampled_batch)
        #
        # photcaled_image = photcaled_batch[0]
        reference_hdu = fits.PrimaryHDU()
        reference_hdu.header = resampled_image.get_header()
        reference_hdu.data = resampled_image.get_data()
        return reference_hdu
        # #  then re-phot-caled)
        # if np.logical_and(self.stack_multiple_images, len(ukirt_images) > 1):
        #     # stack if multiple images are returned
        #     raise NotImplementedError
        #
        # if (len(ukirt_image_urls) == 1) | self.stack_multiple_images:
        #     pass
        #
        # return ukirt_image_hdulists[0]
