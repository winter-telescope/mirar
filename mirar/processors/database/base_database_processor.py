"""
Module containing base database processor class
"""
import logging
from abc import ABC

from mirar.database.constants import POSTGRES_DUPLICATE_PROTOCOLS
from mirar.database.user import PostgresUser
from mirar.paths import max_n_cpu
from mirar.processors.base_processor import BaseProcessor

logger = logging.getLogger(__name__)


class BaseDatabaseProcessor(BaseProcessor, ABC):
    """Base class for processors which interact with a postgres database"""

    max_n_cpu = min(max_n_cpu, 5)

    def __init__(
        self,
        db_name: str,
        db_table: str,
        pg_user: PostgresUser = PostgresUser(),
        duplicate_protocol: str = "fail",
        q3c_bool: bool = False,
    ):
        super().__init__()
        self.db_name = db_name
        self.db_table = db_table

        self.pg_user = pg_user

        self.duplicate_protocol = duplicate_protocol

        assert self.duplicate_protocol in POSTGRES_DUPLICATE_PROTOCOLS

        self.q3c = q3c_bool
