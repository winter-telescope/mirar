"""
Models for database and pydantic dataclass models
"""
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
from winterdrp.processors.database.postgres import DB_USER
from winterdrp.utils.sql import get_engine
from winterdrp.processors.database.postgres import (
    PostgresAdmin,
)

if DB_USER is not None:
    db_name = "summer"
    engine = get_engine(db_name=db_name)

    pg_admin = PostgresAdmin()

    if not pg_admin.check_if_db_exists(db_name=db_name):
        pg_admin.create_db(db_name=db_name)

    Base.metadata.create_all(engine)

    if not default_program.exists():
        default_program.insert_entry()

    populate_fields()
    populate_itid()
    populate_filters()

    from datetime import date

    night_id = date(2002, 1, 1)

    new = Raw(
        nightid=night_id,
        # fid=1
    )
    new.insert_entry()
