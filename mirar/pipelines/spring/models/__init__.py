"""
Models for database and pydantic dataclass models
"""

import logging

from mirar.database.base_table import BaseTable
from mirar.database.credentials import DB_USER
from mirar.database.q3c import create_q3c_extension
from mirar.database.setup import setup_database
from mirar.pipelines.spring.models._diff import Diff, DiffsTable
from mirar.pipelines.spring.models._filters import (
    Filter,
    FiltersTable,
    populate_filters,
)
from mirar.pipelines.spring.models._img_type import (
    ALL_ITID,
    ImgType,
    ImgTypesTable,
    itid_dict,
    populate_itid,
)
from mirar.pipelines.spring.models._programs import (
    DEFAULT_MAX_PRIORITY,
    LEN_PROG_KEY,
    Program,
    ProgramCredentials,
    ProgramsTable,
    default_program,
    populate_programs,
)
from mirar.pipelines.spring.models._raw import Raw, RawsTable
from mirar.pipelines.spring.models._ref_components import (
    RefComponent,
    RefComponentsTable,
)
from mirar.pipelines.spring.models._ref_queries import RefQueriesTable, RefQuery
from mirar.pipelines.spring.models._ref_stacks import RefStack, RefStacksTable
from mirar.pipelines.spring.models._sources import Source, SourcesTable
from mirar.pipelines.spring.models._stack import Stack, StacksTable
from mirar.pipelines.spring.models.base_model import SPRINGBase
from mirar.pipelines.winter.models import set_up_q3c

logger = logging.getLogger(__name__)


def set_up_spring_databases():
    """
    Setup the spring databases

    :return: None
    """

    if DB_USER is not None:
        setup_database(db_base=SPRINGBase)

        for table in [
            RawsTable,
            StacksTable,
        ]:
            set_up_q3c(db_name=SPRINGBase.db_name, db_table=table)

        populate_itid()
        populate_filters()
        populate_programs()

    else:
        logger.warning("No database user provided. Skipping SPRING database setup.")
