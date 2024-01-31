"""
Models for database and pydantic dataclass models
"""

# pylint: disable=duplicate-code
from mirar.database.credentials import DB_USER
from mirar.database.setup import setup_database
from mirar.pipelines.summer.models._diff import Diff, DiffTable
from mirar.pipelines.summer.models._exposures import Exposure, ExposuresTable
from mirar.pipelines.summer.models._fields import (
    DEFAULT_FIELD,
    FieldEntry,
    FieldsTable,
    populate_fields,
)
from mirar.pipelines.summer.models._filters import (
    Filter,
    FiltersTable,
    populate_filters,
)
from mirar.pipelines.summer.models._img_type import (
    ALL_ITID,
    ImgType,
    ImgTypesTable,
    populate_itid,
)
from mirar.pipelines.summer.models._nights import (
    SUMMER_NIGHT_FORMAT,
    Night,
    NightsTable,
)
from mirar.pipelines.summer.models._proc import Proc, ProcTable
from mirar.pipelines.summer.models._programs import (
    DEFAULT_MAX_PRIORITY,
    LEN_PROG_KEY,
    Program,
    ProgramCredentials,
    ProgramsTable,
    default_program,
    populate_programs,
)
from mirar.pipelines.summer.models._raw import Raw, RawTable
from mirar.pipelines.summer.models._subdets import (
    SubDet,
    SubdetsTable,
    populate_subdets,
)
from mirar.pipelines.summer.models.base_model import SummerBase

if DB_USER is not None:
    setup_database(SummerBase)
    populate_fields()
    populate_itid()
    populate_filters()
    populate_programs()
    populate_subdets()
