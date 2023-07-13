"""
Module containing base database processor class
"""
import logging
from abc import ABC
from pathlib import Path
from typing import Optional

from mirar.paths import max_n_cpu
from mirar.processors.base_processor import BaseProcessor
from mirar.processors.database.postgres import (
    POSTGRES_DUPLICATE_PROTOCOLS,
    DataBaseError,
    PostgresAdmin,
    PostgresUser,
)

logger = logging.getLogger(__name__)


class BaseDatabaseProcessor(BaseProcessor, ABC):
    """Base class for processors which interact with a postgres database"""

    max_n_cpu = min(max_n_cpu, 5)

    def __init__(
        self,
        db_name: str,
        db_table: str,
        schema_path: str | Path,
        pg_user: PostgresUser = PostgresUser(),
        pg_admin: PostgresAdmin = PostgresAdmin(),
        has_foreign_keys: bool = False,
        schema_dir: Optional[str | Path] = None,
        duplicate_protocol: str = "fail",
        q3c_bool: bool = False,
    ):
        super().__init__()
        self.db_name = db_name
        self.db_table = db_table

        self.schema_path = schema_path
        self.schema_dir = Path(schema_dir) if schema_dir is not None else None

        self.pg_user = pg_user
        self._pg_admin = pg_admin

        self.has_foreign_keys = has_foreign_keys
        self.db_check_bool = False
        self.duplicate_protocol = duplicate_protocol
        self.q3c = q3c_bool

        assert self.duplicate_protocol in POSTGRES_DUPLICATE_PROTOCOLS

    def db_exists(self) -> bool:
        """
        Checks if a database exists

        :return: boolean
        """
        return self._pg_admin.check_if_db_exists(db_name=self.db_name)

    def _make_db(self):
        """
        Creates a database

        :return:
        """
        self._pg_admin.create_db(db_name=self.db_name)

    def pg_user_exists(self):
        """
        Check if the specified db user exists

        :return: boolean
        """
        return self._pg_admin.check_if_user_exists(self.pg_user.db_user)

    def _make_pg_user(self):
        """
        Creates a new database user

        :return: None
        """
        self._pg_admin.create_new_user(
            new_db_user=self.pg_user.db_user, new_password=self.pg_user.db_password
        )

    def _grant_privileges_to_pg_user(self):
        """
        Grants db privilege to user

        :return: None
        """
        self._pg_admin.grant_privileges(self.db_name, db_user=self.pg_user.db_user)

    def table_exists(self):
        """
        Check if the database table exists

        :return: boolean
        """
        return self.pg_user.check_if_table_exists(
            db_name=self.db_name,
            db_table=self.db_table,
        )

    def _create_table(self, schema_path: str | Path):
        """
        Makes a database table using schema

        :param schema_path: Path of schema file
        :return: None
        """
        self.pg_user.create_table(
            schema_path,
            db_name=self.db_name,
        )

    def check_prerequisites(
        self,
    ):
        if not self.db_check_bool:
            self.check_database_setup()
            self.db_check_bool = True

    def check_database_setup(self):
        """
        Checks the database is set up, and if not, creates the database and tables

        :return: None
        """
        # If not specified, admin credentials equal user credentials
        # Check that the credentials work at all

        self._pg_admin.validate_credentials()

        # Create pg user for the first time

        if not self.pg_user_exists():
            self._make_pg_user()
        self.pg_user.validate_credentials()

        # Check there is a db

        if not self.db_exists():
            self._make_db()
            self._grant_privileges_to_pg_user()

            if self.has_foreign_keys:
                if self.schema_dir is None:
                    self.schema_dir = self.schema_path.parent
                    logger.warning(
                        f"Warning, integrated db setup with foreign keys requested, "
                        f"but no schema directory specified. "
                        f"Will search for schema files in {self.schema_dir} directory."
                    )

                logger.debug(
                    f"Looking in {self.schema_dir} directory to search for schema files"
                )

                self.pg_user.create_tables_from_schema(self.schema_dir, self.db_name)

                if self.q3c:
                    q3c_dir = self.schema_dir.joinpath("q3c")
                    q3c_indexes_file = q3c_dir.joinpath("q3c_indexes.sql")
                    self.pg_user.run_sql_command_from_file(
                        file_path=q3c_indexes_file, db_name=self.db_name
                    )
                    logger.debug("Created q3c_bool indexes")

        if not self.table_exists():
            self._create_table(self.schema_path)

            if self.q3c:
                q3c_dir = self.schema_dir.joinpath("q3c")
                table_q3c_path = q3c_dir.joinpath(f"q3c_{self.db_table}.sql")

                if not table_q3c_path.exists():
                    err = (
                        f"q3c_bool extension requested but no "
                        f"{table_q3c_path} file found. Please add it in."
                    )
                    logger.error(err)
                    raise DataBaseError(err)

                self.pg_user.run_sql_command_from_file(
                    file_path=table_q3c_path, db_name=self.db_name
                )
