"""
Models for database and pydantic dataclass models
"""
from typing import Union

from sqlalchemy.orm import DeclarativeBase

from winterdrp.pipelines.winter.models._exposures import Exposures, ExposuresTable
from winterdrp.pipelines.winter.models._fields import (
    Fields,
    FieldsTable,
    populate_fields,
)
from winterdrp.pipelines.winter.models._filters import (
    Filters,
    FiltersTable,
    populate_filters,
)
from winterdrp.pipelines.winter.models._imgType import (
    ALL_ITID,
    ImgTypes,
    ImgTypesTable,
    populate_itid,
)
from winterdrp.pipelines.winter.models._nights import Nights, NightsTable
from winterdrp.pipelines.winter.models._proc import Proc, ProcTable
from winterdrp.pipelines.winter.models._programs import (
    ProgramCredentials,
    Programs,
    ProgramsTable,
    populate_programs,
)
from winterdrp.pipelines.winter.models._raw import Raw, RawTable
from winterdrp.pipelines.winter.models._ref_components import (
    RefComponents,
    RefComponentsTable,
)
from winterdrp.pipelines.winter.models._ref_queries import RefQueries, RefQueriesTable
from winterdrp.pipelines.winter.models._ref_stacks import RefStacks, RefStacksTable
from winterdrp.pipelines.winter.models._subdets import Subdets, SubdetsTable
from winterdrp.pipelines.winter.models.base_model import WinterBase
from winterdrp.processors.database.postgres import (
    ADMIN_PASSWORD,
    ADMIN_USER,
    DB_PASSWORD,
    DB_USER,
    PostgresAdmin,
)
from winterdrp.processors.sqldatabase.base_model import BaseTable
from winterdrp.utils.sql import get_engine


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
