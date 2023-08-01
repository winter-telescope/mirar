"""
Models for database and pydantic dataclass models
"""
import logging
from typing import Union

from psycopg import OperationalError
from sqlalchemy.orm import DeclarativeBase

from mirar.pipelines.winter.models._astrometry_stats import (
    AstrometryStat,
    AstrometryStatsTable,
)
from mirar.pipelines.winter.models._exposures import Exposure, ExposuresTable
from mirar.pipelines.winter.models._fields import (
    DEFAULT_FIELD,
    FieldEntry,
    FieldsTable,
    populate_fields,
)
from mirar.pipelines.winter.models._filters import (
    Filter,
    FiltersTable,
    populate_filters,
)
from mirar.pipelines.winter.models._img_type import (
    ALL_ITID,
    ImgType,
    ImgTypesTable,
    itid_dict,
    populate_itid,
)
from mirar.pipelines.winter.models._nights import Night, NightsTable
from mirar.pipelines.winter.models._programs import (
    DEFAULT_MAX_PRIORITY,
    LEN_PROG_KEY,
    Program,
    ProgramCredentials,
    ProgramsTable,
    default_program,
    populate_programs,
)
from mirar.pipelines.winter.models._raw import Raw, RawTable
from mirar.pipelines.winter.models._ref_components import (
    RefComponent,
    RefComponentsTable,
)
from mirar.pipelines.winter.models._ref_queries import RefQueriesTable, RefQuery
from mirar.pipelines.winter.models._ref_stacks import RefStack, RefStacksTable
from mirar.pipelines.winter.models._stack import Stack, StacksTable
from mirar.pipelines.winter.models._subdets import (
    Subdet,
    SubdetsTable,
    populate_subdets,
)
from mirar.pipelines.winter.models.base_model import WinterBase
from mirar.processors.sqldatabase.base_model import BaseTable
from mirar.processors.sqldatabase.postgres import (
    ADMIN_PASSWORD,
    ADMIN_USER,
    DB_PASSWORD,
    DB_USER,
    PostgresAdmin,
    PostgresUser,
)
from mirar.utils.sql import get_engine

logger = logging.getLogger(__name__)


def setup_database(db_base: Union[DeclarativeBase, BaseTable]):
    """
    Function to setup database

    :param db_base: BaseTable
    :return: None
    """
    db_name = db_base.db_name
    engine = get_engine(db_name=db_name)

    pg_user = PostgresUser()

    try:
        pg_user.validate_credentials()
    except OperationalError:
        logger.warning(
            f"Failed to credentials for user {DB_USER}. " f"Will try creating new user."
        )
        pg_admin = PostgresAdmin()
        pg_admin.validate_credentials()
        pg_admin.create_new_user(new_db_user=DB_USER, new_password=DB_PASSWORD)
        pg_user.validate_credentials()

    if not pg_user.check_if_db_exists(db_name=db_name):
        pg_user.create_db(db_name=db_name)

    db_base.metadata.create_all(engine)


if DB_USER is not None:
    setup_database(db_base=WinterBase)

    populate_fields()
    populate_itid()
    populate_filters()
    populate_programs()
    populate_subdets()
