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
from astropy.io.fits import PrimaryHDU
from astropy.units import Quantity
from astropy.wcs import WCS
from astroquery.ukidss import UkidssClass
from astroquery.utils.commons import FileContainer
from astrosurveyutils import get_known_ukirt_surveys
from astrosurveyutils.surveys import MOCSurvey

from mirar.data import Image, ImageBatch
from mirar.errors import ProcessorError
from mirar.io import open_fits
from mirar.paths import (
    BASE_NAME_KEY,
    COADD_KEY,
    LATEST_SAVE_KEY,
    LATEST_WEIGHT_SAVE_KEY,
    core_fields,
    get_output_path,
)
from mirar.pipelines.winter.constants import winter_filters_map
from mirar.processors.astromatic.sextractor.sextractor import Sextractor
from mirar.processors.astromatic.swarp.swarp import Swarp
from mirar.processors.base_processor import ImageHandler
from mirar.processors.candidates.utils import (
    get_corners_ra_dec_from_header,
    get_image_center_wcs_coords,
)
from mirar.processors.photcal import PhotCalibrator
from mirar.processors.split import SUB_ID_KEY
from mirar.processors.sqldatabase.base_model import BaseDB
from mirar.processors.sqldatabase.database_exporter import DatabaseImageExporter
from mirar.references.base_reference_generator import BaseReferenceGenerator

logger = logging.getLogger(__name__)

ukirt_image_height = 90 * u.arcmin
ukirt_image_width = 90 * u.arcmin


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

    ra_list, dec_list = wcs.all_pix2world(xcrd_list, ycrd_list, 1)
    ra_list[ra_list < 0] = ra_list[ra_list < 0] + 360
    return ra_list, dec_list


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
    ) = get_ukirt_file_identifiers_from_url(url)
    basename = (
        f"{multiframeid}_{extension_id}_{frame_lx}_"
        f"{frame_hx}_{frame_ly}_{frame_hy}.fits"
    )

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
    image.header["MULTIFRAME_ID"] = multiframeid
    image.header["EXTENSION_ID"] = extension_id
    image.header["LX"] = frame_lx
    image.header["HX"] = frame_hx
    image.header["LY"] = frame_ly
    image.header["HY"] = frame_hy

    image.header["UKIRPATH"] = ukirt_filename

    return image


def get_ukirt_file_identifiers_from_url(url: str) -> list:
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


def check_query_exists_locally(
    query_ra, query_dec, query_filt, query_table, components_table
):
    """
    Function to check if component images exist locally based on the query_ra
    and query_dec
    Args:
        query_ra:
        query_dec:
        db_table:

    Returns:

    """
    results = query_table.sql_model().select_query(
        select_keys="compid",
        compare_values=[query_ra, query_dec, query_filt],
        compare_keys=["query_ra", "query_dec", "query_filt"],
        comparators=["__eq__", "__eq__", "__eq__"],
    )
    logger.debug(results)
    savepaths = []
    # TODO - run a joint query to get the savepaths
    if len(results) > 0:
        savepaths = []
        compids = [x[0] for x in results]
        for compid in compids:
            comp_results = components_table.sql_model().select_query(
                select_keys="savepath",
                compare_values=[compid],
                compare_keys=["compid"],
                comparators=["__eq__"],
            )
            if len(comp_results) == 0:
                raise ValueError(f"Component {compid} not found in database")
            savepaths.append(comp_results[0][0])
    return savepaths


def check_multiframe_exists_locally(
    db_table, multiframe_id, extension_id, frame_lx, frame_hx, frame_ly, frame_hy
):
    """
    Function to query database to check if a multiframe exists locally
    Args:
        db_table:
        multiframe_id:
        extension_id:
        frame_lx:
        frame_hx:
        frame_ly:
        frame_hy:

    Returns:

    """
    results = db_table.sql_model().select_query(
        select_keys=["savepath"],
        compare_values=[
            multiframe_id,
            extension_id,
            frame_lx,
            frame_hx,
            frame_ly,
            frame_hy,
        ],
        compare_keys=["multiframe_id", "extension_id", "lx", "hx", "ly", "hy"],
        comparators=["__eq__", "__eq__", "__eq__", "__eq__", "__eq__", "__eq__"],
    )
    logger.debug(results)
    if len(results) == 0:
        savepaths = []
    else:
        savepaths = [x[0] for x in results]
    return savepaths


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
        sextractor_generator: Callable[..., Sextractor] = None,
        phot_calibrator_generator: Callable[..., PhotCalibrator] = None,
        check_local_database: bool = True,
        write_to_db: bool = True,
        query_table: Type[BaseDB] = None,
        components_table: Type[BaseDB] = None,
        write_db_table: Type[BaseDB] = None,
        component_image_dir: str = None,
        night_sub_dir: str = None,
    ):
        super().__init__(
            filter_name=filter_name, write_to_db=write_to_db, db_table=write_db_table
        )
        self.stack_multiple_images = stack_multiple_images
        self.swarp_resampler = swarp_resampler
        self.num_query_points = num_query_points
        self.sextractor_generator = sextractor_generator
        self.phot_calibrator_generator = phot_calibrator_generator
        if np.logical_and(self.stack_multiple_images, self.swarp_resampler is None):
            err = (
                "Stacking mutiple reference images requested, but no swarp resampler "
                "suppied"
            )
            raise UKIRTRefError(err)

        self.check_local_database = check_local_database
        self.components_table = components_table
        self.query_table = query_table
        if np.logical_and(self.check_local_database, self.components_table is None):
            err = (
                "You have requested checking locally, but no database "
                "has been specified"
            )
            raise UKIRTRefError(err)

        self.component_image_dir = component_image_dir
        self.night_sub_dir = night_sub_dir

    def get_reference(self, image: Image) -> tuple[PrimaryHDU, PrimaryHDU]:
        header = image.get_header()

        query_ra_list, query_dec_list = get_query_coordinates_from_header(
            header, numpoints=self.num_query_points
        )

        query_crds = SkyCoord(ra=query_ra_list, dec=query_dec_list, unit=(u.deg, u.deg))
        logger.debug(f"Querying around {query_crds}")

        ukirt_query = UkidssClass()

        query_ra_cent, query_dec_cent = get_image_center_wcs_coords(image, origin=1)
        logger.debug(f"Center RA: {query_ra_cent} Dec: {query_dec_cent}")
        ukirt_surveys = find_ukirt_surveys(
            query_ra_cent, query_dec_cent, self.filter_name
        )
        if len(ukirt_surveys) == 0:
            err = "Coordinates not in UKIRT surveys"
            raise NotinUKIRTError(err)

        # Sort surveys in descending order of limiting mags
        lim_mags = [x.lim_mag for x in ukirt_surveys]
        ukirt_surveys = ukirt_surveys[np.argsort(lim_mags)[::-1]]
        logger.debug(f"Surveys are {[x.survey_name for x in ukirt_surveys]}")
        (
            ukirt_image_urls,
            ukirt_query_ras,
            ukirt_query_decs,
            ukirt_query_exists_locally_list,
        ) = ([], [], [], [])

        for survey in ukirt_surveys:
            ukirt_query.database = survey.wfau_dbname
            # url_list = ukirt_query.get_image_list(query_crds,
            #                                       image_width=image_width,
            #                                       image_height=image_height
            #                                       )
            # This seems to pull multiextension images, which we don't really want.
            # url_list = ukirt_query.get_image_list(center_crds, radius=search_radius)
            url_list, ukirt_qra_list, ukirt_qdec_list, query_exists_list = (
                [],
                [],
                [],
                [],
            )
            for crd in query_crds:
                # Need to add a cache and check there.
                url, query_exists = [], False
                if self.check_local_database:
                    url = check_query_exists_locally(
                        query_ra=crd.ra.deg,
                        query_dec=crd.dec.deg,
                        query_filt=self.filter_name,
                        query_table=self.query_table,
                        components_table=self.components_table,
                    )
                    logger.debug(f"Found {len(url)} images locally.")
                    query_exists = len(url) > 0

                if len(url) == 0:
                    url = ukirt_query.get_image_list(
                        crd,
                        image_width=ukirt_image_width,
                        image_height=ukirt_image_height,
                        waveband=self.filter_name,
                    )

                qexists_list = [query_exists] * len(url)
                ra_list, dec_list = [crd.ra.deg] * len(url), [crd.dec.deg] * len(url)
                url_list += url
                ukirt_qra_list += ra_list
                ukirt_qdec_list += dec_list
                query_exists_list += qexists_list

            ukirt_image_urls += url_list
            ukirt_query_ras += ukirt_qra_list
            ukirt_query_decs += ukirt_qdec_list
            ukirt_query_exists_locally_list += query_exists_list

        logger.debug(
            f"UKIRT image url length {len(ukirt_image_urls)}. "
            f"List {ukirt_image_urls}"
        )

        if len(ukirt_image_urls) == 0:
            err = "No image found at the given coordinates in the UKIRT database"
            raise UKIRTRefNotFoundError(err)

        if not self.stack_multiple_images:
            # Choose the image from the most sensitive survey,
            # if you do not wish to stack
            ukirt_image_urls = [ukirt_image_urls[0]]

        ukirt_images = []
        for url, ra, dec, qexists in zip(
            ukirt_image_urls,
            ukirt_query_ras,
            ukirt_query_decs,
            ukirt_query_exists_locally_list,
        ):
            if "http" in url:
                if self.check_local_database:
                    (
                        _,
                        multiframe_id,
                        extension_id,
                        frame_lx,
                        frame_hx,
                        frame_ly,
                        frame_hy,
                    ) = get_ukirt_file_identifiers_from_url(url)
                    imgpath = check_multiframe_exists_locally(
                        db_table=self.components_table,
                        multiframe_id=multiframe_id,
                        extension_id=extension_id,
                        frame_lx=frame_lx,
                        frame_hx=frame_hx,
                        frame_ly=frame_ly,
                        frame_hy=frame_hy,
                    )
                    if len(imgpath) > 0:
                        url = imgpath[0]
                        # TODO: Decide whether we should make an entry in the database
                        # Probably not, otherwise there will be multiple entries per
                        # image.

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
                    url=url,
                )
                savepath = get_output_path(
                    ukirt_image[BASE_NAME_KEY], dir_root=self.component_image_dir
                )
                ukirt_image["QUERY_RA"] = ra
                ukirt_image["QUERY_DEC"] = dec
                ukirt_image[LATEST_SAVE_KEY] = savepath.as_posix()
                if self.check_local_database:
                    dbexporter = DatabaseImageExporter(
                        db_table=self.components_table,
                        duplicate_protocol=self.duplicate_protocol,
                        q3c_bool=self.q3c_bool,
                    )
                    ukirt_db_batch = dbexporter.apply(ImageBatch([ukirt_image]))
                    ukirt_image = ukirt_db_batch[0]

                self.save_fits(ukirt_image, savepath)
            else:
                logger.debug(f"Reading image from {url}")
                with fits.open(url, ignore_missing_simple=True) as ukirt_hdulist:
                    ukirt_image = Image(
                        header=ukirt_hdulist[0].header,  # pylint: disable=no-member
                        data=ukirt_hdulist[0].data,  # pylint: disable=no-member
                    )

            if self.check_local_database & ~qexists:
                new = self.query_table(
                    query_ra=ra,
                    query_dec=dec,
                    query_filt=self.filter_name,
                    compid=ukirt_image["COMPID"],
                )
                new.insert_entry()

            ukirt_images.append(ukirt_image)

        # Only keep unique images at the end of image collection
        ukirt_images = np.array(ukirt_images)
        ukirt_image_names = [x[BASE_NAME_KEY] for x in ukirt_images]
        _, unique_image_indexes = np.unique(ukirt_image_names, return_index=True)
        ukirt_images = ukirt_images[unique_image_indexes]

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
        for ind, image_to_resamp in enumerate(ukirt_images):
            image_to_resamp["FLXSCALE"] = scaling_factors[ind]

        ukirt_image_batch = ImageBatch(list(ukirt_images))
        resampler = self.swarp_resampler(
            center_ra=query_ra_cent,
            center_dec=query_dec_cent,
            include_scamp=False,
            combine=True,
            calculate_dims_in_swarp=True,
            # x_imgpixsize=25 * 60 / 0.40,
            # y_imgpixsize=25 * 60 / 0.40,
        )
        resampler.set_night(night_sub_dir=self.night_sub_dir)
        resampled_batch = resampler.apply(ukirt_image_batch)

        resampled_image = resampled_batch[0]
        resamp_ra_cent, resamp_dec_cent = get_image_center_wcs_coords(
            image=resampled_image, origin=1
        )

        resampled_image["RA_CENT"] = resamp_ra_cent
        resampled_image["DEC_CENT"] = resamp_dec_cent

        compids = [str(x.header["compid"]) for x in ukirt_images]
        resampled_image["COMPIDS"] = ",".join(compids)

        ref_sextractor = self.sextractor_generator(resampled_image)
        ref_sextractor.set_night(night_sub_dir=self.night_sub_dir)
        resampled_batch = ref_sextractor.apply(resampled_batch)
        reference_weight_path = resampled_batch[0].header[LATEST_WEIGHT_SAVE_KEY]
        reference_weight_data, reference_weight_header = open_fits(
            reference_weight_path
        )

        phot_calibrator = self.phot_calibrator_generator(resampled_batch[0])
        phot_calibrator.set_night(night_sub_dir=self.night_sub_dir)
        phot_calibrator.set_preceding_steps([ref_sextractor])
        photcaled_batch = phot_calibrator.apply(resampled_batch)

        photcaled_image = photcaled_batch[0]
        reference_hdu = fits.PrimaryHDU()
        reference_hdu.header = photcaled_image.get_header()
        reference_hdu.data = photcaled_image.get_data()

        reference_weight_hdu = fits.PrimaryHDU()
        reference_weight_hdu.header = reference_weight_header
        reference_weight_hdu.data = reference_weight_data

        ra_cent, dec_cent = get_image_center_wcs_coords(image=photcaled_image, origin=1)

        (
            (ra0_0, dec0_0),
            (ra0_1, dec0_1),
            (ra1_0, dec1_0),
            (ra1_1, dec1_1),
        ) = get_corners_ra_dec_from_header(reference_hdu.header)

        if self.write_to_db:
            # for key in ["FIELDID", "SUBDETID"]:
            #     if key not in image.header:
            #         image.header[key] = 0
            stackid = (
                f"{str(image.header['FIELDID']).rjust(5, '0')}"
                f"{str(image.header[SUB_ID_KEY]).rjust(2, '0')}"
                f"{str(winter_filters_map[image.header['FILTER']])}"
            )
            reference_hdu.header["STACKID"] = int(stackid)
        reference_hdu.header["RA_CENT"] = ra_cent
        reference_hdu.header["DEC_CENT"] = dec_cent
        reference_hdu.header["RA0_0"] = ra0_0
        reference_hdu.header["RA0_1"] = ra0_1
        reference_hdu.header["RA1_0"] = ra1_0
        reference_hdu.header["RA1_1"] = ra1_1
        reference_hdu.header["DEC0_0"] = dec0_0
        reference_hdu.header["DEC0_1"] = dec0_1
        reference_hdu.header["DEC1_0"] = dec1_0
        reference_hdu.header["DEC1_1"] = dec1_1

        return reference_hdu, reference_weight_hdu
