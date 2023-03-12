"""
Models for database and pydantic dataclass models
"""
from datetime import date, datetime

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
from winterdrp.pipelines.summer.models.basemodel import Base
from winterdrp.processors.database.postgres import DB_PASSWORD, DB_USER, PostgresAdmin
from winterdrp.utils.sql import get_engine

if DB_USER is not None:
    db_name = "summertest"
    engine = get_engine(db_name=db_name)

    pg_admin = PostgresAdmin()

    if not pg_admin.check_if_db_exists(db_name=db_name):
        pg_admin.create_db(db_name=db_name)

    if not pg_admin.check_if_user_exists(user_name=DB_USER):
        pg_admin.create_new_user(new_db_user=DB_USER, new_password=DB_PASSWORD)
        pg_admin.grant_privileges(db_name=db_name, db_user=DB_USER)

    Base.metadata.create_all(engine)

    if not default_program.exists():
        default_program.insert_entry()

    # populate_fields()
    populate_itid()
    populate_filters()

    night_id = date(2002, 1, 1)
    rawid = 1000001
    new = Raw(
        rawid=rawid,
        nightid=night_id,
        fid=1,
        expid=rawid,
        savepath=f"/Users/viraj/raw/{rawid}",
        obsdate=date(year=2020, month=10, day=2),
        timeutc=datetime(year=2020, month=10, day=2),
        obsID=1,
        itid=1,
        AExpTime=30,
        expMJD=59001.32,
        airmass=1.2,
        shutopen=datetime(year=2020, month=10, day=2),
        shutclsd=datetime(year=2020, month=10, day=2),
        altitude=30.0,
        azimuth=265.5,
        ra=284.66587,
        dec=45.5667,
    )
    #
    # new.insert_entry()
