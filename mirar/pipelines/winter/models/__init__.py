"""
Models for database and pydantic dataclass models
"""
from typing import Union

from sqlalchemy.orm import DeclarativeBase

from mirar.pipelines.winter.constants import NXSPLIT, NYSPLIT
from mirar.pipelines.winter.models._astrometry_stats import (
    AstrometryStats,
    AstrometryStatsTable,
)
from mirar.pipelines.winter.models._exposures import Exposures, ExposuresTable
from mirar.pipelines.winter.models._fields import (
    DEFAULT_FIELD,
    Fields,
    FieldsTable,
    populate_fields,
)
from mirar.pipelines.winter.models._filters import (
    Filters,
    FiltersTable,
    populate_filters,
)
from mirar.pipelines.winter.models._imgType import (
    ALL_ITID,
    ImgTypes,
    ImgTypesTable,
    itid_dict,
    populate_itid,
)
from mirar.pipelines.winter.models._nights import Nights, NightsTable
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
    RefComponents,
    RefComponentsTable,
)
from mirar.pipelines.winter.models._ref_queries import RefQueries, RefQueriesTable
from mirar.pipelines.winter.models._ref_stacks import RefStacks, RefStacksTable
from mirar.pipelines.winter.models._stack import Stacks, StacksTable
from mirar.pipelines.winter.models._subdets import (
    Subdets,
    SubdetsTable,
    populate_subdets,
)
from mirar.pipelines.winter.models.base_model import WinterBase
from mirar.processors.database.postgres import (
    ADMIN_PASSWORD,
    ADMIN_USER,
    DB_PASSWORD,
    DB_USER,
    PostgresAdmin,
)
from mirar.processors.sqldatabase.base_model import BaseTable
from mirar.utils.sql import get_engine


def setup_database(base: Union[DeclarativeBase, BaseTable]):
    """
    Function to setup database
    Args:
        base:
    Returns:
    """
    if DB_USER is not None:
        db_name = base.db_name
        admin_engine = get_engine(
            db_name=db_name, db_user=ADMIN_USER, db_password=ADMIN_PASSWORD
        )

        pg_admin = PostgresAdmin()

        if not pg_admin.check_if_db_exists(db_name=db_name):
            pg_admin.create_db(db_name=db_name)

        base.metadata.create_all(
            admin_engine
        )  # extensions need to be created as a superuser

        if not pg_admin.check_if_user_exists(user_name=DB_USER):
            pg_admin.create_new_user(new_db_user=DB_USER, new_password=DB_PASSWORD)

        pg_admin.grant_privileges(db_name=db_name, db_user=DB_USER)


if DB_USER is not None:
    setup_database(base=WinterBase)

    populate_fields()
    populate_itid()
    populate_filters()
    populate_programs()
    populate_subdets()
