"""
Module to query for WFCAM images.
You can either query the online WFAU archive, or query a local database to get
component images.
"""
import logging
import warnings
from pathlib import Path
from typing import Callable, Type

import numpy as np
import pandas as pd
from astropy import units as u
from astropy.coordinates import SkyCoord
from astropy.io import fits
from astropy.units import Quantity
from astropy.wcs import FITSFixedWarning
from astroquery.ukidss import UkidssClass
from astroquery.utils.commons import FileContainer
from astroquery.wfau import BaseWFAUClass
from astrosurveyutils.surveys import MOCSurvey

from mirar.data import Image, ImageBatch
from mirar.data.utils import check_coords_within_image, get_image_center_wcs_coords
from mirar.database.base_model import BaseDB
from mirar.database.constraints import DBQueryConstraints
from mirar.database.transactions import select_from_table
from mirar.errors import ProcessorError
from mirar.io import open_raw_image
from mirar.paths import BASE_NAME_KEY, LATEST_SAVE_KEY, get_output_dir, get_output_path
from mirar.processors.database import DatabaseImageInserter
from mirar.references.wfcam.files import wfcam_undeprecated_compid_file
from mirar.references.wfcam.utils import (
    COMPID_KEY,
    EXTENSION_ID_KEY,
    MULTIFRAME_ID_KEY,
    QUERY_DEC_KEY,
    QUERY_FILT_KEY,
    QUERY_RA_KEY,
    find_ukirt_surveys,
    get_query_coordinates_from_header,
    get_wfcam_file_identifiers_from_url,
    make_wfcam_image_from_hdulist,
    open_compressed_wfcam_fits,
    save_wfcam_as_compressed_fits,
)

logger = logging.getLogger(__name__)

wfau_image_height = 90 * u.arcmin
wfau_image_width = 90 * u.arcmin


class WFAURefError(ProcessorError):
    """
    Base UKIRTRef error
    """


class NotinWFCAMError(ProcessorError):
    """
    Error when the coordinates are not in WFAU footprint
    """


class WFAUQueryDBError(ProcessorError):
    """
    Error related to the databases while querying
    """


class WFCAMRefNotFoundError(ProcessorError):
    """
    Error when WFCAM ref is not found for some reason
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
        :param filter_name: Filter name to query for.
        :param num_query_points: Number of points to use to define the query region. The
            image is divided into np.sqrt(num_query_points) x np.sqrt(num_query_points)
            regions.
        :param query_coords_function: Function to get the query coordinates from the
        header.
            The function should take a header and the number of query points as input,
            and return a list of tuples of coordinates.
        :param component_image_subdir: Subdirectory to save component images to.
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

    def run_query(self, image: Image) -> ImageBatch:
        """
        Run the query for the given image
        :param image: Image to query for
        :return: ImageBatch containing the queried images
        """
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
    savepath (saved path) , mfid (multiframeid), xtnsnid (extension_id),
    (paramters used to uniquely identify a WFCAM image).
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
            :param query_coords_function: Function to use to get the query coordinates
            from the header.
            :param use_db_for_component_queries: Whether to use local databases to
            perform queries. This is useful if you want to reduce the number of queries
            to the
            online database. If set, the code assumes that you are storing the
            individual images in a `components_db_table` and also the details of every
            query to a `query_db_table`.
            :param components_db_table: Table with the details of the individual WFCAM
            single extension images. The following keys need to be present in the db, as
            they uniquely determine a WFCAM image:
            mfid, xtnsnid, compid the primary key should be compid.
            query_db_table: Table with the details of the queries to the WFCAM database.
            The following keys need to be present in the db : qry_ra, qry_dec,
            qry_filt and compid.
            :param skip_online_query: Whether to skip the online query and only use the
            local
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
        self.dbexporter = DatabaseImageInserter(db_table=self.query_db_table)

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

            required_components_db_keys = [
                MULTIFRAME_ID_KEY.lower(),
                EXTENSION_ID_KEY.lower(),
                COMPID_KEY.lower(),
            ]
            for key in required_components_db_keys:
                if key not in self.components_db_table.sql_model.__table__.columns:
                    raise ValueError(
                        f"{key} must be present in the components_db_table"
                    )
            required_query_db_keys = [
                QUERY_RA_KEY.lower(),
                QUERY_DEC_KEY.lower(),
                QUERY_FILT_KEY.lower(),
                COMPID_KEY.lower(),
            ]
            for key in required_query_db_keys:
                if key not in self.query_db_table.sql_model.__table__.columns:
                    raise ValueError(f"{key} must be present in the query_db_table")

    def get_query_class(self) -> BaseWFAUClass:
        """
        Get the class that will be used to query the WFAU database, e.g. VSAClass or
        UKIDSSClass
        :return: Class that will be used to query the WFAU database
        """
        raise NotImplementedError

    def get_surveys(self, ra: float, dec: float) -> list[MOCSurvey]:
        """
        Get the surveys that are available at the given coordinates
        :param ra: RA of the coordinates
        :param dec: Dec of the coordinates
        :return: List of surveys that are available at the given coordinates
        """
        raise NotImplementedError

    def get_query_crds(
        self, header: fits.Header, num_query_points: int
    ) -> tuple[list[float], list[float]]:
        """
        Get the query coordinates from the header.
        :param header: Header of the image.
        :param num_query_points: Number of points to use to define the query region.
        The image is divided into np.sqrt(num_query_points) x np.sqrt(num_query_points)
        regions.
        :return: Tuple of lists of RA and Dec coordinates.
        """
        return self.query_coords_function(header, num_query_points)

    def run_wfau_query(
        self,
        image: Image,
    ) -> ImageBatch:
        """
        Run the query to the WFAU database.
        :param image: Image to query around.
        :return: ImageBatch with the results of the query.
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
            raise NotinWFCAMError(err)
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
            for ind, crd in enumerate(query_crds):
                logger.debug(f"Running query {ind}/{len(query_crds)}")
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
                    query_exists = len(imagepaths) > 0

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
                if self.use_db_for_component_queries & (not query_exists):
                    queried_images = [
                        open_raw_image(imagepath, open_f=open_compressed_wfcam_fits)
                        for imagepath in imagepaths
                    ]
                    for img in queried_images:
                        img[QUERY_RA_KEY] = crd.ra.deg
                        img[QUERY_DEC_KEY] = crd.dec.deg
                        img[QUERY_FILT_KEY] = self.filter_name
                    self.dbexporter.apply(ImageBatch(queried_images))

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
            raise WFCAMRefNotFoundError(err)

        wfcam_images = ImageBatch(
            [
                open_raw_image(url, open_f=open_compressed_wfcam_fits)
                for url in wfcam_image_paths
            ]
        )

        return wfcam_images

    def run_query(self, image: Image) -> ImageBatch:
        """
        Run the query on the UKIRT server.
        :param image: Image object with the coordinates of the image.
        :return: ImageBatch object with the downloaded images.
        """
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
    undeprecated_compids_file: Path = wfcam_undeprecated_compid_file,
) -> list[Path]:
    """
    Download the image from UKIRT server. Optionally, check if the image exists locally
    and ingest it into a database.
    :param crd: SkyCoord object with the coordinates of the image.
    :param wfau_query: WFAU query object.
    :param survey_name: Name of the survey to query.
    :param waveband: Waveband of the image.
    :param save_dir_path: Path to the directory where the image will be saved.
    :param image_width: Width of the image to download.
    :param image_height: Height of the image to download.
    :param use_local_database: If True, check if the image exists locally and ingest it
    into a database.
    :param components_table: Table to use for the components database.
    :param duplicate_protocol: Protocol to follow if the image already exists locally.
    :param q3c_bool: Is q3c setup?
    :param undeprecated_compids_file: Path to the file with the list of undeprecated
    component ids.
    :return imagepaths: List of paths to the downloaded images.
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
    for url_ind, url in enumerate(url_list):
        logger.debug(f"Downloading {url_ind}/{len(url_list)}")
        local_imagepaths = []
        (
            ukirt_filename,
            multiframe_id,
            extension_id,
            _,
            _,
            _,
            _,
        ) = get_wfcam_file_identifiers_from_url(url)

        # Check if image is deprecated. If so, don't use it.
        if undeprecated_compids_file is not None:
            compid = int(f"{multiframe_id}{extension_id}")
            undeprecated_compids = pd.read_csv(undeprecated_compids_file)[
                "COMPID"
            ].values
            if compid not in undeprecated_compids:
                logger.debug(
                    f"File with multiframeid {multiframe_id} and "
                    f"extension {extension_id} is deprecated. Skipping."
                )
                continue

        if use_local_database:
            # Check if the image exists locally.
            local_imagepaths = check_multiframe_exists_locally(
                db_table=components_table,
                multiframe_id=multiframe_id,
                extension_id=extension_id,
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
                ukirt_hdulist=wfcam_img_hdulist,
                ukirt_filename=ukirt_filename,
                multiframeid=multiframe_id,
                extension_id=extension_id,
            )
            imagepath = get_output_path(
                wfcam_image[BASE_NAME_KEY], dir_root=save_dir_path.as_posix()
            )
            wfcam_image[QUERY_RA_KEY] = crd.ra.deg
            wfcam_image[QUERY_DEC_KEY] = crd.dec.deg
            wfcam_image[QUERY_FILT_KEY] = waveband
            wfcam_image[LATEST_SAVE_KEY] = imagepath.as_posix()

            if use_local_database:
                dbexporter = DatabaseImageInserter(
                    db_table=components_table,
                    duplicate_protocol=duplicate_protocol,
                )
                wfcam_db_batch = dbexporter.apply(ImageBatch([wfcam_image]))
                wfcam_image = wfcam_db_batch[0]

            save_wfcam_as_compressed_fits(wfcam_image, imagepath)
            logger.debug(f"Saved UKIRT image to {imagepath}")

        imagepaths.append(imagepath)

    return np.unique(imagepaths).tolist()


def check_query_exists_locally(
    query_ra: float,
    query_dec: float,
    query_filt: str,
    query_table: Type[BaseDB],
    components_table: Type[BaseDB],
) -> list[Path]:
    """
    Function to check if component images exist locally based on the query_ra
    and query_dec
    Args:
        :param query_ra: ra that was queried
        :param query_dec: dec that was queried
        :param query_filt: filter that was queried
        :param query_table: table with query details
        :param components_table: table with component image details
    Returns:
        :return: list of savepaths
    """

    constraints = DBQueryConstraints(
        columns=[QUERY_FILT_KEY],
        accepted_values=[query_filt],
    )
    constraints.add_q3c_constraint(
        ra=query_ra,
        dec=query_dec,
        crossmatch_radius_arcsec=10.0,
        ra_field_name=QUERY_RA_KEY,
        dec_field_name=QUERY_DEC_KEY,
    )
    results = select_from_table(
        db_constraints=constraints,
        sql_table=query_table.sql_model,
        output_columns=[COMPID_KEY.lower()],
    )

    logger.debug(results)
    savepaths = []
    if len(results) > 0:
        savepaths = []
        compids = results[COMPID_KEY.lower()].tolist()
        for compid in compids:
            constraints = DBQueryConstraints(
                columns=[COMPID_KEY],
                accepted_values=[compid],
            )
            comp_results = select_from_table(
                db_constraints=constraints,
                sql_table=components_table.sql_model,
                output_columns=["savepath"],
            )
            if len(comp_results) == 0:
                raise WFAUQueryDBError(
                    f"Component {compid} not found in database, but "
                    "a query corresponding to it exists. The query"
                    "table is likely out of sync with the component"
                )
            savepaths.append(Path(comp_results["savepath"].iloc[0]))
    return savepaths


def get_locally_existing_overlap_images(
    query_ra: float, query_dec: float, query_filt: str, components_table: Type[BaseDB]
) -> list[Path]:
    """
    Function to get the locally existing images that overlap with the given coordinates
    Args:
        :param query_ra: ra that was queried
        :param query_dec: dec that was queried
        :param query_filt: filter that was queried
        :param components_table: table with component image details

    Returns:
        :return: list of savepaths
    """

    # Get around RA=0/360 issue
    if (query_ra > 0.46) & (query_ra < 359.57):
        constraints = DBQueryConstraints(
            columns=["ramin", "ramax", "decmin", "decmax", "filter"],
            accepted_values=[
                (0.46, query_ra),
                (query_ra, 359.57),
                query_dec,
                query_dec,
                query_filt,
            ],
            comparison_types=["between", "between", "<=", ">=", "="],
        )
    else:
        constraints = DBQueryConstraints(
            columns=["ramin", "ramax", "decmin", "decmax", "filter"],
            accepted_values=[query_ra, query_ra, query_dec, query_dec, query_filt],
            comparison_types=[">=", "<=", ">=", "<=", "="],
        )
    logger.debug(f"Constraints: {constraints.parse_constraints()}")
    results = select_from_table(
        db_constraints=constraints,
        sql_table=components_table.sql_model,
        output_columns=["savepath"],
    )

    logger.debug(results)
    savepaths = []
    if len(results) > 0:
        savepaths = [Path(x) for x in results["savepath"].tolist()]
        # Confirm that the coordinates are in the image
        logger.debug(f"Checking WCS of {len(savepaths)} images for overlap")
        with warnings.catch_warnings():
            warnings.simplefilter("ignore", FITSFixedWarning)
            savepaths = [
                x
                for x in savepaths
                if check_coords_within_image(
                    header=fits.getheader(x, 1), ra=query_ra, dec=query_dec
                )
            ]
        logger.debug(f"{len(savepaths)} images confirmed to overlap")
    return savepaths


def check_multiframe_exists_locally(
    db_table: Type[BaseDB],
    multiframe_id: int,
    extension_id: int,
) -> list[Path]:
    """
    Function to query database to check if a multiframe exists locally
    Args:
        :param db_table: table with multiframe details
        :param multiframe_id: multiframe id
        :param extension_id: extension id

    Returns:
        :return: list of savepaths
    """

    db_constraints = DBQueryConstraints(
        columns=[MULTIFRAME_ID_KEY.lower(), EXTENSION_ID_KEY.lower()],
        accepted_values=[multiframe_id, extension_id],
    )

    results = select_from_table(
        db_constraints=db_constraints,
        sql_table=db_table.sql_model,
        output_columns=["savepath"],
    )

    logger.debug(results)
    if len(results) == 0:
        savepaths = []
    else:
        savepaths = [Path(x) for x in results["savepath"].tolist()]
    return savepaths


class UKIRTOnlineQuery(WFAUQuery):
    """
    Class to query the UKIRT online database at the WFAU.
    This is a subclass of the WFAUQuery.
    """

    def get_surveys(self, ra: float, dec: float) -> list[MOCSurvey]:
        """
        Function to get the surveys that overlap with the given coordinates
        Args:
            :param ra: ra that was queried
            :param dec: dec that was queried
        Returns:
            :return: list of surveys
        """
        return find_ukirt_surveys(ra=ra, dec=dec, band=self.filter_name)

    def get_query_class(self) -> BaseWFAUClass:
        """
        Function to get the query class
        Returns:
            :return: query class
        """
        return UkidssClass()
