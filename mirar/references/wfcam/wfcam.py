"""
Module for querying reference images from the UKIRT survey
"""
import logging
from collections.abc import Callable
from pathlib import Path
from typing import Type

import astropy.units as u
import numpy as np
from astropy.coordinates import SkyCoord
from astropy.io import fits
from astropy.units import Quantity
from astroquery.utils.commons import FileContainer
from astroquery.wfau import BaseWFAUClass
from astrosurveyutils.surveys import MOCSurvey
from tqdm import tqdm

from mirar.data import Image, ImageBatch
from mirar.errors import ProcessorError
from mirar.io import open_raw_image, save_fits
from mirar.paths import (
    BASE_NAME_KEY,
    LATEST_SAVE_KEY,
    ZP_KEY,
    ZP_STD_KEY,
    get_output_dir,
    get_output_path,
)
from mirar.processors.astromatic.sextractor.sextractor import Sextractor
from mirar.processors.astromatic.swarp import Swarp
from mirar.processors.base_processor import ImageHandler
from mirar.processors.photcal import PhotCalibrator
from mirar.processors.sqldatabase.base_model import BaseDB
from mirar.processors.sqldatabase.database_exporter import DatabaseImageExporter
from mirar.references.base_reference_generator import BaseStackReferenceGenerator
from mirar.references.wfcam.utils import (
    get_image_center_wcs_coords,
    get_query_coordinates_from_header,
    get_wfcam_file_identifiers_from_url,
    make_wfcam_image_from_hdulist,
)

logger = logging.getLogger(__name__)

wfau_image_height = 90 * u.arcmin
wfau_image_width = 90 * u.arcmin


class WFAURefError(ProcessorError):
    """
    Base UKIRTRef error
    """


class NotinWFAUError(ProcessorError):
    """
    Error when the coordinates are not in UKIRT footprint
    """


class WFAURefNotFoundError(ProcessorError):
    """
    Error when UKIRT ref is not found for some reason
    """


def check_query_exists_locally(
    query_ra, query_dec, query_filt, query_table, components_table
):
    """
    Function to check if component images exist locally based on the query_ra
    and query_dec
    Args:
        query_ra: ra that was queried
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


def get_locally_existing_overlap_images(
    query_ra, query_dec, query_filt, components_table
):
    """
    Function to get the locally existing images that overlap with the given coordinates
    Args:
        query_ra:
        query_dec:
        query_filt:
        components_table:

    Returns:

    """
    savepaths = []
    results = components_table.sql_model().select_query(
        select_keys=["savepath"],
        compare_values=[query_ra, query_ra, query_dec, query_dec, query_filt],
        compare_keys=["ra0_0", "ra1_1", "dec0_0", "dec1_1", "filter"],
        comparators=["__le__", "__ge__", "__ge__", "__le__", "__eq__"],
    )
    logger.debug(results)
    if len(results) > 0:
        savepaths = [x[0] for x in results]
    return savepaths


def check_multiframe_exists_locally(
    db_table, multiframe_id, extension_id, frame_lx, frame_hx, frame_ly, frame_hy
) -> list:
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


def default_filter_wfau_images(image_batch: ImageBatch) -> ImageBatch:
    """
    Function to filter WFAU images
    Args:
        image_batch:

    Returns:

    """
    image_array = np.array([x for x in image_batch])

    mag_zps = np.array(
        [
            x["MAGZPT"]
            + 2.5 * np.log10(x["EXPTIME"])
            - x["EXTINCT"] * ((x["AMSTART"] + x["AMEND"]) / 2)
            for x in image_batch
        ]
    )
    # magerr_zps = np.array([x["MAGZRR"] for x in ukirt_images])
    median_mag_zp = np.median(mag_zps)
    seeings = np.array([x["SEEING"] for x in image_batch])
    zpmask = np.abs(mag_zps - median_mag_zp) < 0.4
    seeingmask = (seeings < 2.5) & (seeings > 0)

    image_array = image_array[zpmask & seeingmask]

    return ImageBatch(image_array.tolist())


def download_wfcam_archive_images(
    crd: SkyCoord,
    wfau_query,
    survey_name: str,
    waveband: str,
    save_dir_path: Path,
    image_width: Quantity = wfau_image_width,
    image_height: Quantity = wfau_image_height,
    use_local_database: bool = False,
    components_table: Type[BaseDB] = None,
    duplicate_protocol: str = "ignore",
    q3c_bool: bool = False,
) -> list[Path]:
    """
    Download the image from UKIRT server. Optionally, check if the image exists locally
    and ingest it into a database.
    """
    # ukirt_query = UkidssClass()
    wfau_query.database = survey_name

    # First get a list with details of the images that overlap with the
    # coordinates.
    url_list = wfau_query.get_image_list(
        crd,
        image_width=image_width,
        image_height=image_height,
        waveband=waveband,
    )

    imagepaths = []
    for url in url_list:
        local_imagepaths = []
        if use_local_database:
            # Check if the image exists locally.
            (
                _,
                multiframe_id,
                extension_id,
                frame_lx,
                frame_hx,
                frame_ly,
                frame_hy,
            ) = get_wfcam_file_identifiers_from_url(url)

            local_imagepaths = check_multiframe_exists_locally(
                db_table=components_table,
                multiframe_id=multiframe_id,
                extension_id=extension_id,
                frame_lx=frame_lx,
                frame_hx=frame_hx,
                frame_ly=frame_ly,
                frame_hy=frame_hy,
            )

        image_exists_locally = len(local_imagepaths) > 0

        if image_exists_locally:
            imagepath = local_imagepaths[0]
        else:
            # Download the actual image. This is copied from what happens in
            # astroquery
            obj = FileContainer(
                url,
                encoding="binary",
                remote_timeout=wfau_query.TIMEOUT,
                show_progress=True,
            )
            wfcam_img_hdulist = obj.get_fits()

            # UKIRT ref images are stored as multiHDU files, need to combine the
            # hdus so no info from the headers is lost. This also adds in core_fields.
            wfcam_image = make_wfcam_image_from_hdulist(
                wfcam_img_hdulist,
                url=url,
            )
            imagepath = get_output_path(
                wfcam_image[BASE_NAME_KEY], dir_root=save_dir_path.as_posix()
            )
            wfcam_image["QUERY_RA"] = crd.ra.deg
            wfcam_image["QUERY_DEC"] = crd.dec.deg
            wfcam_image["QUERY_FILT"] = waveband
            wfcam_image[LATEST_SAVE_KEY] = imagepath.as_posix()

            if use_local_database:
                dbexporter = DatabaseImageExporter(
                    db_table=components_table,
                    duplicate_protocol=duplicate_protocol,
                    q3c_bool=q3c_bool,
                )
                wfcam_db_batch = dbexporter.apply(ImageBatch([wfcam_image]))
                wfcam_image = wfcam_db_batch[0]

            save_fits(wfcam_image, imagepath)
            logger.debug(f"Saved UKIRT image to {local_imagepaths}")

        imagepaths.append(imagepath)

    return np.unique(imagepaths).tolist()


class WFAUQuery:
    def __init__(
        self,
        num_query_points: int = 4,
        query_coords_function: Callable[
            [fits.Header, int], (list[float], list[float])
        ] = get_query_coordinates_from_header,
        use_db_for_component_queries: bool = False,
        components_db_table: Type[BaseDB] = None,
        query_db_table: Type[BaseDB] = None,
        skip_online_query: bool = False,
        filter_images: Callable[[ImageBatch], ImageBatch] = default_filter_wfau_images,
        num_images_per_query: int = None,
    ):
        self.num_query_points = num_query_points
        self.components_db_table = components_db_table
        self.query_db_table = query_db_table
        self.use_db_for_component_queries = use_db_for_component_queries
        self.query_coords_function = query_coords_function
        self.skip_online_query = skip_online_query
        self.filter_images = filter_images
        self.num_images_per_query = num_images_per_query

        if self.use_db_for_component_queries:
            if self.components_db_table is None:
                raise ValueError(
                    "components_table must be provided if "
                    "check_local_database is True"
                )
            if self.query_db_table is None:
                raise ValueError(
                    "query_table must be provided if check_local_database is True"
                )

    def get_query_crds(self, header, num_query_points: int):
        return self.query_coords_function(header, num_query_points)

    def run_wfau_query(
        self,
        wfau_query: BaseWFAUClass,
        query_crds: SkyCoord,
        wfau_survey_names: list[str],
        filter_name: str,
        savedir: Path,
    ) -> ImageBatch:
        (
            wfcam_image_paths,
            wfau_query_ras,
            wfau_query_decs,
            wfau_query_exists_locally_list,
        ) = ([], [], [], [])
        for survey in wfau_survey_names:
            wfau_query.database = survey
            paths_list, wfau_qra_list, wfau_qdec_list, query_exists_list = (
                [],
                [],
                [],
                [],
            )
            for ind in tqdm(range(len(query_crds))):
                crd = query_crds[ind]
                # Need to add a cache and check there.
                imagepaths, query_exists = [], False
                if self.use_db_for_component_queries:
                    # First, check if the exact coordinates have been queried to UKIRT
                    # server before.
                    imagepaths = check_query_exists_locally(
                        query_ra=crd.ra.deg,
                        query_dec=crd.dec.deg,
                        query_filt=filter_name,
                        query_table=self.query_db_table,
                        components_table=self.components_db_table,
                    )
                    logger.debug(f"Found {len(imagepaths)} images locally.")

                    # If no query found, check if the coordinates overlap with any of
                    # the component images present in the database. This is a hack to
                    # avoid failures when the server is down
                    if len(imagepaths) == 0:
                        imagepaths = get_locally_existing_overlap_images(
                            query_ra=crd.ra.deg,
                            query_dec=crd.dec.deg,
                            query_filt=filter_name,
                            components_table=self.components_db_table,
                        )
                        logger.debug(
                            f"Found {len(imagepaths)} component images containing "
                            f"the coordinates locally."
                        )

                    query_exists = len(imagepaths) > 0

                # If no query found locally, download from the UKIRT server.
                # This runs only is skip_online_query is False, again, as a safeguard
                # against cases where the server is out for long times.
                if (len(imagepaths) == 0) and (not self.skip_online_query):
                    imagepaths = download_wfcam_archive_images(
                        crd,
                        wfau_query=wfau_query,
                        survey_name=survey,
                        waveband=filter_name,
                        save_dir_path=savedir,
                        use_local_database=self.use_db_for_component_queries,
                        components_table=self.components_db_table,
                        duplicate_protocol="ignore",
                    )

                    # Make an entry in the queries table
                    if self.use_db_for_component_queries:
                        downloaded_images = [
                            open_raw_image(imagepath) for imagepath in imagepaths
                        ]
                        dbexporter = DatabaseImageExporter(db_table=self.query_db_table)
                        _ = dbexporter.apply(ImageBatch(downloaded_images))

                qexists_list = [query_exists] * len(imagepaths)
                ra_list, dec_list = [crd.ra.deg] * len(imagepaths), [crd.dec.deg] * len(
                    imagepaths
                )

                paths_list += imagepaths
                wfau_qra_list += ra_list
                wfau_qdec_list += dec_list
                query_exists_list += qexists_list

            wfcam_image_paths += paths_list
            wfau_query_ras += wfau_qra_list
            wfau_query_decs += wfau_qdec_list
            wfau_query_exists_locally_list += query_exists_list

        wfcam_image_paths = list(set(wfcam_image_paths))
        logger.debug(
            f"UKIRT image url length {len(wfcam_image_paths)}. "
            f"List {wfcam_image_paths}"
        )
        if len(wfcam_image_paths) == 0:
            err = "No image found at the given coordinates in the UKIRT database"
            raise WFAURefNotFoundError(err)

        wfcam_images = ImageBatch([open_raw_image(url) for url in wfcam_image_paths])

        return wfcam_images


class WFCAMStackedRefOnline(BaseStackReferenceGenerator, ImageHandler, WFAUQuery):
    """
    Class to query UKIRT images from the WFAU archive and stack them together
    """

    abbreviation = "wfcam_stack_online_ref"

    def __init__(
        self,
        filter_name: str,
        image_resampler_generator: Callable[..., Swarp],
        num_query_points: int = 4,
        write_stack_to_db: bool = False,
        stacks_db_table: Type[BaseDB] = None,
        use_db_for_component_queries: bool = False,
        components_db_table: Type[BaseDB] = None,
        query_db_table: Type[BaseDB] = None,
        component_image_sub_dir: str = None,
        references_base_subdir_name: str = "references",
        skip_online_query: bool = False,
        stack_image_annotator: Callable[[Image], Image] = None,
        photcal_stack: bool = False,
        sextractor_generator: Callable[..., Sextractor] = None,
        phot_calibrator_generator: Callable[..., PhotCalibrator] = None,
    ):
        BaseStackReferenceGenerator.__init__(
            self,
            filter_name=filter_name,
            write_stack_to_db=write_stack_to_db,
            stacks_db_table=stacks_db_table,
            references_base_subdir_name=references_base_subdir_name,
            image_resampler_generator=image_resampler_generator,
            stack_image_annotator=stack_image_annotator,
            photcal_stack=photcal_stack,
            sextractor_generator=sextractor_generator,
            phot_calibrator_generator=phot_calibrator_generator,
        )
        WFAUQuery.__init__(
            self,
            use_db_for_component_queries=use_db_for_component_queries,
            skip_online_query=skip_online_query,
            components_db_table=components_db_table,
            query_db_table=query_db_table,
            num_query_points=num_query_points,
        )

        if self.use_db_for_component_queries:
            if (self.components_db_table is None) | (self.query_db_table is None):
                err = (
                    "You have requested checking locally, but no database "
                    "has been specified"
                )
                raise WFAURefError(err)

        self.component_image_dir = component_image_sub_dir
        self.references_base_subdir_name = references_base_subdir_name

    def get_surveys(self, ra: float, dec: float) -> list[MOCSurvey]:
        raise NotImplementedError

    def get_query_class(self) -> BaseWFAUClass:
        raise NotImplementedError

    def get_component_images(self, image: Image) -> ImageBatch:
        header = image.get_header()

        query_ra_list, query_dec_list = self.get_query_crds(
            header, self.num_query_points
        )

        query_crds = SkyCoord(ra=query_ra_list, dec=query_dec_list, unit=(u.deg, u.deg))
        logger.debug(f"Querying around {query_crds}")

        query_ra_cent, query_dec_cent = get_image_center_wcs_coords(image, origin=1)
        logger.debug(f"Center RA: {query_ra_cent} Dec: {query_dec_cent}")

        # Get different surveys by the telescope
        surveys = self.get_surveys(query_ra_cent, query_dec_cent)
        if len(surveys) == 0:
            err = "Coordinates not in any survey"
            raise NotinWFAUError(err)
            # Sort surveys in descending order of limiting mags
        lim_mags = [x.lim_mag for x in surveys]
        surveys = surveys[np.argsort(lim_mags)[::-1]]
        logger.debug(f"Surveys are {[x.survey_name for x in surveys]}")

        # Get the query class
        wfau_query = self.get_query_class()
        component_dir_path = get_output_dir(
            dir_root=self.component_image_dir, sub_dir=self.references_base_subdir_name
        )
        wfau_survey_names = [x.wfau_dbname for x in surveys]
        wfau_images = self.run_wfau_query(
            query_crds=query_crds,
            wfau_survey_names=wfau_survey_names,
            wfau_query=wfau_query,
            filter_name=self.filter_name,
            savedir=component_dir_path,
        )

        wfau_images = self.filter_images(wfau_images)

        # change BASENAME to gel well with parallel processing
        for ind, ref_img in enumerate(wfau_images):
            new_basename = (
                f"{ref_img[BASE_NAME_KEY].strip('.fits')}" f"_{image[BASE_NAME_KEY]}"
            )
            ref_img[BASE_NAME_KEY] = new_basename

        # Get the scaling factors
        mag_zps = np.array(
            [
                x["MAGZPT"]
                + 2.5 * np.log10(x["EXPTIME"])
                - x["EXTINCT"] * ((x["AMSTART"] + x["AMEND"]) / 2)
                for x in wfau_images
            ]
        )

        median_mag_zp = np.median(mag_zps)
        scaling_factors = 10 ** (0.4 * (median_mag_zp - mag_zps))

        for ind, image_to_resamp in enumerate(wfau_images):
            image_to_resamp["FLXSCALE"] = scaling_factors[ind]
            image_to_resamp[ZP_KEY] = median_mag_zp
            image_to_resamp[ZP_STD_KEY] = np.std(mag_zps)

        wfau_image_batch = ImageBatch(list(wfau_images))

        return wfau_image_batch
