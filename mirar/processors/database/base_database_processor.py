"""
Module containing base database processor class
"""

import logging
from abc import ABC
from typing import Optional, Type

from mirar.database.base_model import BaseDB
from mirar.database.user import PostgresAdmin, PostgresUser
from mirar.processors.base_processor import BaseProcessor

logger = logging.getLogger(__name__)


class BaseDatabaseProcessor(BaseProcessor, ABC):
    """Base class for processors which interact with a postgres database"""

    max_n_cpu = 1

    def __init__(
        self,
        db_table: Type[BaseDB],
        pg_user: Optional[PostgresUser] = None,
        pg_admin: Optional[PostgresAdmin] = None,
    ):
        super().__init__()
        self.db_table = db_table
        self.db_name = self.db_table.sql_model.db_name

        # Instantiated here rather than as a default argument, so that
        # credentials are read from the environment at run time.
        self.pg_user = pg_user if pg_user is not None else PostgresUser()
        self._pg_admin = pg_admin if pg_admin is not None else PostgresAdmin()
