"""
Models for database and pydantic dataclass models
"""
from winterdrp.pipelines.summer.models._program import (
    Program,
    ProgramCredentials,
    ProgramTable,
)
from winterdrp.pipelines.summer.models.basemodel import Base
from winterdrp.processors.database.postgres import DB_USER
from winterdrp.utils.sql import get_engine

if DB_USER is not None:
    engine = get_engine()
    Base.metadata.create_all(engine)


if __name__ == "__main__":
    p = Program(
        progname="2019A000",
        prog_key="???",
        piname="drwho",
        progid=1,
        progtitle="???",
        startdate="2019-01-01",
        enddate="2020-01-01",
        hours_allocated=100,
        hours_remaining=100,
        basepriority=100,
    )
    p.insert_entry()
