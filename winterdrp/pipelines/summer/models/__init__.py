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
from winterdrp.pipelines.summer.models._programs import (
    ProgramCredentials,
    Programs,
    ProgramsTable,
    default_program,
)
from winterdrp.pipelines.summer.models.basemodel import Base
from winterdrp.processors.database.postgres import DB_USER
from winterdrp.utils.sql import get_engine

if DB_USER is not None:
    engine = get_engine()
    Base.metadata.create_all(engine)

    if not default_program.exists():
        default_program.insert_entry()

    populate_fields()
    populate_itid()
    populate_filters()
