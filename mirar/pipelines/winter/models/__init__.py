"""
Models for database and pydantic dataclass models
"""
import logging
from typing import Type, Union

from sqlalchemy.exc import OperationalError
from sqlalchemy.orm import DeclarativeBase

from mirar.database.base_model import BaseTable
from mirar.database.engine import get_engine
from mirar.database.q3c import create_q3c_extension
from mirar.database.credentials import DB_PASSWORD, DB_USER
from mirar.database.user import PostgresAdmin, PostgresUser
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

logger = logging.getLogger(__name__)


def set_up_q3c(db_name: str, db_table: BaseTable):
    """
    Function to setup q3c extension for a given table in db

    :param db_name: Name of database
    :param db_table: Table to setup q3c extension for
    :return:
    """
    create_q3c_extension(
        db_name=db_name,
        table_name=db_table.__tablename__,
        ra_column_name=db_table.ra_column_name,
        dec_column_name=db_table.dec_column_name,
    )


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
            f"Failed to credentials for user {DB_USER}. "
            f"Will try creating new user with this name using admin credentials."
        )
        pg_admin = PostgresAdmin()
        pg_admin.validate_credentials()
        pg_admin.create_new_user(new_db_user=DB_USER, new_password=DB_PASSWORD)
        pg_user.validate_credentials()

    pg_user.create_db(db_name=db_name)

    db_base.metadata.create_all(engine)

    if not pg_user.has_extension(extension_name="q3c", db_name=db_name):
        logger.info(
            f"No 'q3c' extension found. Creating it now for db "
            f"{db_name} with admin user."
        )

        pg_admin = PostgresAdmin()
        pg_admin.validate_credentials()
        pg_admin.create_extension(extension_name="q3c", db_name=db_name)


if DB_USER is not None:
    setup_database(db_base=WinterBase)

    for table in [ExposuresTable]:
        set_up_q3c(db_name=WinterBase.db_name, db_table=table)

    populate_fields()
    populate_itid()
    populate_filters()
    populate_programs()
    populate_subdets()
