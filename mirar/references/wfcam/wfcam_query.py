"""
Module to query for WFCAM images.
You can either query the online WFAU archive, or query a local database to get
component images.
"""
import logging
from pathlib import Path
from typing import Callable, Type

import numpy as np
from astropy import units as u
from astropy.coordinates import SkyCoord
from astropy.io import fits
from astropy.units import Quantity
from astroquery.ukidss import UkidssClass
from astroquery.utils.commons import FileContainer
from astroquery.wfau import BaseWFAUClass
from astrosurveyutils.surveys import MOCSurvey
from tqdm import tqdm

from mirar.data import Image, ImageBatch
from mirar.data.utils import get_image_center_wcs_coords
from mirar.errors import ProcessorError
from mirar.io import open_raw_image, save_fits
from mirar.paths import BASE_NAME_KEY, LATEST_SAVE_KEY, get_output_dir, get_output_path
from mirar.processors.sqldatabase.base_model import BaseDB
from mirar.processors.sqldatabase.database_exporter import DatabaseImageExporter
from mirar.references.wfcam.utils import (
    find_ukirt_surveys,
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


class BaseWFCAMQuery:
    """
    Base class for querying WFCAM images
    """

    def __init__(
        self,
        filter_name: str,
        num_query_points: int = 4,
        query_coords_function: Callable[
            [fits.Header, int], tuple[list[float], list[float]]
        ] = get_query_coordinates_from_header,
        component_image_subdir: str | Path = None,
    ):
        """
        num_query_points: Number of points to use to define the query region. The
            image is divided into np.sqrt(num_query_points) x np.sqrt(num_query_points)
            regions.
        query_coords_function: Function to get the query coordinates from the header.
            The function should take a header and the number of query points as input,
            and return a list of tuples of coordinates.
        component_image_subdir: Subdirectory to save component images to.
        """
        self.num_query_points = num_query_points
        self.query_coords_function = query_coords_function
        self.component_image_subdir = component_image_subdir
        if self.component_image_subdir is not None:
            self.savedir = get_output_dir(
                dir_root=self.component_image_subdir,
            )
        else:
            self.savedir = None
        if isinstance(self.savedir, str):
            self.savedir = Path(self.savedir)
        self.filter_name = filter_name

    def run_query(self, image: Image):
        raise NotImplementedError


class WFAUQuery(BaseWFCAMQuery):
    """
    Class to handle queries to the online WFAU archive. To reduce the number of
    queries to the server, the user can optionally choose to set up databases. If
    this is chosen, this script currently assumes the following database structure ->
    1. Two tables : query_db_table and components_db_table
    a. query_db_table : This table stores the details of the queries. The following
    columns are required - query_ra, query_dec, query_filt, compid (primary key of the
    table entry of the image downloaded).
    b. components_db_table : This table stores the details of the individual downloaded
    images. The following columns are required - compid (primary_key),
    savepath (saved path) , multiframeid, extension_id, frame_lx,
    frame_ly, frame_hx, frame_hy (paramters used to uniquely identify a WFCAM image).
    It is recommended to use the database model files from
    mirar/pipelines/winter/models/ref_queries.py and
    mirar/pipelines/winter/models/ref_components.py to set up the tables in your
    database.

    1. The user specifies an image and the filter to query and optionally the
    database details.
    2. The image is broken down into user-specified number of coordinates to get
    overlapping images from the archive.
    3. If the user has specified database details, each coordinate is checked against
    the query database to see if it has been queried before. If so, the corresponding
    component images from the comoponent_db_table are used.
    4. If not, the query is made to the WFAU server to get the URLs of the images. The
    details of each image are extracted from the URL.
    5. If the user has specified database details, the image details are xmatched to the
    database to see if the image has been downloaded before. If so, the corresponding
    image is used.
    6. If not, the image is downloaded and saved to the user-specified path.
    7. If the user has specified database details, the details of the downloaded image
    and the query are ingested into the respective tables.
    """

    def __init__(
        self,
        filter_name: str,
        num_query_points: int = 4,
        query_coords_function: Callable[
            [fits.Header, int],
            tuple[list[float], list[float]],
        ] = get_query_coordinates_from_header,
        component_image_subdir: str = "wfau_components",
        use_db_for_component_queries: bool = False,
        components_db_table: Type[BaseDB] = None,
        query_db_table: Type[BaseDB] = None,
        skip_online_query: bool = False,
    ):
        """
        Parameters:
            query_coords_function: Function to use to get the query coordinates from the
            header.
            use_db_for_component_queries: Whether to use local databases to perform
            queries. This is useful if you want to reduce the number of queries to the
            online database. If set, the code assumes that you are storing the
            individual images in a `components_db_table` and also the details of every
            query to a `query_db_table`.
            components_db_table: Table with the details of the individual WFCAM
            single extension images. The following keys need to be present in the db, as
            they uniquely determine a WFCAM image:
            multiframe_id, extension_id, frame_lx, frame_hx, frame_ly, frame_hy, the
            primary key should be compid.
            query_db_table: Table with the details of the queries to the WFCAM database.
            The following keys need to be present in the db : query_ra, query_dec,
            query_filt and compid.
            skip_online_query: Whether to skip the online query and only use the local
            databases.
        """
        super().__init__(
            num_query_points=num_query_points,
            query_coords_function=query_coords_function,
            component_image_subdir=component_image_subdir,
            filter_name=filter_name,
        )
        self.components_db_table = components_db_table
        self.query_db_table = query_db_table
        self.use_db_for_component_queries = use_db_for_component_queries
        self.skip_online_query = skip_online_query

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

    def get_query_class(self) -> BaseWFAUClass:
        """
        Get the class that will be used to query the WFAU database, e.g. VSAClass or
        UKIDSSClass
        """
        raise NotImplementedError

    def get_surveys(self, ra: float, dec: float) -> list[MOCSurvey]:
        """
        Get the surveys that are available at the given coordinates
        """
        raise NotImplementedError

    def get_query_crds(
        self, header: fits.Header, num_query_points: int
    ) -> tuple[list[float], list[float]]:
        """
        Get the query coordinates from the header.
        header: Header of the image.
        num_query_points: Number of points to use to define the query region. The image
        is divided into np.sqrt(num_query_points) x np.sqrt(num_query_points) regions.
        """
        return self.query_coords_function(header, num_query_points)

    def run_wfau_query(
        self,
        image: Image,
    ) -> ImageBatch:
        """
        Run the query to the WFAU database.
        wfau_query: WFAU query object.
        query_crds: Coordinates to query.
        wfau_survey_names: WFAU names of the survey databases to query.
        filter_name: Name of the filter to query.
        savedir: Directory to save the images.
        """
        query_ra_list, query_dec_list = self.get_query_crds(
            image.header, self.num_query_points
        )

        query_crds = SkyCoord(ra=query_ra_list, dec=query_dec_list, unit=(u.deg, u.deg))
        logger.debug(f"Querying around {query_crds}")

        query_ra_cent, query_dec_cent = get_image_center_wcs_coords(image, origin=1)
        logger.debug(f"Center RA: {query_ra_cent} Dec: {query_dec_cent}")

        (
            wfcam_image_paths,
            wfau_query_ras,
            wfau_query_decs,
            wfau_query_exists_locally_list,
        ) = ([], [], [], [])

        # Get different surveys by the telescope
        query_ra_cent = np.median(query_crds.ra.deg)
        query_dec_cent = np.median(query_crds.dec.deg)
        surveys = self.get_surveys(query_ra_cent, query_dec_cent)
        if len(surveys) == 0:
            err = "Coordinates not in any survey"
            raise NotinWFAUError(err)
            # Sort surveys in descending order of limiting mags
        lim_mags = [x.lim_mag for x in surveys]
        surveys = surveys[np.argsort(lim_mags)[::-1]]
        logger.debug(f"Surveys are {[x.survey_name for x in surveys]}")
        wfau_survey_names = [x.wfau_dbname for x in surveys]
        # Get the query class
        wfau_query = self.get_query_class()

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
                        query_filt=self.filter_name,
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
                            query_filt=self.filter_name,
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
                        waveband=self.filter_name,
                        save_dir_path=self.savedir,
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

    def run_query(self, image: Image):
        return self.run_wfau_query(image=image)


def download_wfcam_archive_images(
    crd: SkyCoord,
    wfau_query: BaseWFAUClass,
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
            imagepath = Path(local_imagepaths[0])
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


def check_query_exists_locally(
    query_ra, query_dec, query_filt, query_table, components_table
):
    """
    Function to check if component images exist locally based on the query_ra
    and query_dec
    Args:
        query_ra: ra that was queried
        query_dec: dec that was queried
        query_filt: filter that was queried
        query_table: table with query details
        components_table: table with component image details
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
    query_ra: float, query_dec: float, query_filt: str, components_table: Type[BaseDB]
) -> list[str]:
    """
    Function to get the locally existing images that overlap with the given coordinates
    Args:
        query_ra: ra that was queried
        query_dec: dec that was queried
        query_filt: filter that was queried
        components_table: table with component image details

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
    db_table: Type[BaseDB],
    multiframe_id: int,
    extension_id: int,
    frame_lx: int,
    frame_hx: int,
    frame_ly: int,
    frame_hy: int,
) -> list[str]:
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


class UKIRTOnlineQuery(WFAUQuery):
    """
    Class to query the UKIRT online database at the WFAU.
    This is a subclass of the WFAUQuery.
    """

    def get_surveys(self, ra: float, dec: float) -> list[MOCSurvey]:
        return find_ukirt_surveys(ra=ra, dec=dec, band=self.filter_name)

    def get_query_class(self) -> BaseWFAUClass:
        return UkidssClass()
