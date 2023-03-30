"""
Models for database and pydantic dataclass models
"""
from datetime import date, datetime

from sqlalchemy.orm import Session

from winterdrp.pipelines.summer.models._dithers import Dithers, DithersTable
from winterdrp.pipelines.summer.models._fields import (
    Fields,
    FieldsTable,
    populate_fields,
)
from winterdrp.pipelines.summer.models._filters import (
    Filters,
    FiltersTable,
    populate_filters,
)
from winterdrp.pipelines.summer.models._itid import (
    ALL_ITID,
    ITIDs,
    ITIDsTable,
    populate_itid,
)
from winterdrp.pipelines.summer.models._nights import Nights, NightsTable
from winterdrp.pipelines.summer.models._programs import (
    ProgramCredentials,
    Programs,
    ProgramsTable,
    default_program,
)
from winterdrp.pipelines.summer.models._raw import Raw, RawTable
from winterdrp.pipelines.summer.models.basemodel import SummerBase
from winterdrp.processors.database.postgres import (
    ADMIN_PASSWORD,
    ADMIN_USER,
    DB_PASSWORD,
    DB_USER,
    PostgresAdmin,
)
from winterdrp.utils.sql import get_engine

if DB_USER is not None:
    db_name = SummerBase.db_name
    admin_engine = get_engine(
        db_name=db_name, db_user=ADMIN_USER, db_password=ADMIN_PASSWORD
    )

    engine = get_engine(db_name=db_name)
    pg_admin = PostgresAdmin()

    if not pg_admin.check_if_db_exists(db_name=db_name):
        pg_admin.create_db(db_name=db_name)

    SummerBase.metadata.create_all(
        admin_engine
    )  # extensions need to be created as a superuser

    if not pg_admin.check_if_user_exists(user_name=DB_USER):
        pg_admin.create_new_user(new_db_user=DB_USER, new_password=DB_PASSWORD)

    pg_admin.grant_privileges(db_name=db_name, db_user=DB_USER)

    session = Session(bind=engine)
    if not default_program.exists():
        default_program.insert_entry()

    populate_fields()
    populate_itid()
    populate_filters()
