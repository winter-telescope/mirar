"""
Module containing base database processor class
"""
import logging
import os
from abc import ABC
from pathlib import Path
from typing import Optional

import numpy as np

from winterdrp.processors.base_processor import BaseProcessor
from winterdrp.processors.database.postgres import (
    PG_ADMIN_PWD_KEY,
    PG_ADMIN_USER_KEY,
    POSTGRES_DUPLICATE_PROTOCOLS,
    DataBaseError,
    check_if_db_exists,
    check_if_table_exists,
    check_if_user_exists,
    create_db,
    create_new_user,
    create_table,
    create_tables_from_schema,
    grant_privileges,
    run_sql_command_from_file,
)

logger = logging.getLogger(__name__)


class BaseDatabaseProcessor(BaseProcessor, ABC):
    """Base class for processors which interact with a postgres database"""

    def __init__(
        self,
        db_name: str,
        db_table: str,
        schema_path: str | Path,
        db_user: str = os.environ.get("DB_USER"),
        db_password: str = os.environ.get("DB_PWD"),
        full_setup: bool = False,
        schema_dir: Optional[str | Path] = None,
        duplicate_protocol: str = "fail",
        q3c_bool: bool = False,
    ):
        super().__init__()
        self.db_name = db_name
        self.db_table = db_table
        self.db_user = db_user
        self.db_password = db_password
        self.full_setup = full_setup
        self.schema_path = schema_path
        self.schema_dir = Path(schema_dir) if schema_dir is not None else None
        self.db_check_bool = False
        self.duplicate_protocol = duplicate_protocol
        self.q3c = q3c_bool
        assert self.duplicate_protocol in POSTGRES_DUPLICATE_PROTOCOLS

    def db_exists(self) -> bool:
        """
        Checks if a database exists

        :return: boolean
        """
        return check_if_db_exists(db_name=self.db_name)

    def make_db(self):
        """
        Creates a database

        :return:
        """
        create_db(db_name=self.db_name)

    def user_exists(self):
        """
        Check if the specified db user exists

        :return: boolean
        """
        return check_if_user_exists(self.db_user)

    def make_user(self):
        """
        Creates a new database user

        :return: None
        """
        create_new_user(new_db_user=self.db_user, new_password=self.db_password)

    def grant_privileges(self):
        """
        Grants db privilege to user

        :return: None
        """
        grant_privileges(self.db_name, self.db_user)

    def table_exists(self):
        """
        Check if the database table exists

        :return: boolean
        """
        return check_if_table_exists(
            db_name=self.db_name,
            db_table=self.db_table,
            db_user=self.db_user,
            password=self.db_password,
        )

    def make_table(self, schema_path: str | Path):
        """
        Makes a database table using schema

        :param schema_path: Path of schema file
        :return: None
        """
        create_table(
            schema_path,
            db_name=self.db_name,
            db_user=self.db_user,
            password=self.db_password,
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

        admin_user = os.environ.get(PG_ADMIN_USER_KEY)
        admin_password = os.environ.get(PG_ADMIN_PWD_KEY)

        if np.logical_and(self.db_exists(), np.invert(self.user_exists())):
            err = "Database exists but user does not exist"
            logger.error(err)
            raise DataBaseError(err)

        if not self.db_exists():
            self.make_db()

            if not self.user_exists():
                self.make_user()

            self.grant_privileges()

            if self.full_setup:

                if self.schema_dir is None:
                    self.schema_dir = self.schema_path.parent
                    logger.warning(
                        f"Warning, full db setup requested, "
                        f"but no schema directory specified. "
                        f"Will search for schema files in {self.schema_dir} directory."
                    )

                logger.info(
                    f"Looking in {self.schema_dir} directory to search for schema files"
                )

                create_tables_from_schema(
                    self.schema_dir, self.db_name, self.db_user, self.db_password
                )

                if self.q3c:
                    q3c_dir = self.schema_dir.joinpath("q3c")
                    q3c_indexes_file = q3c_dir.joinpath("q3c_indexes.sql")
                    run_sql_command_from_file(
                        file_path=q3c_indexes_file,
                        db_name=self.db_name,
                        db_user=admin_user,
                        password=admin_password,
                    )
                    logger.info("Created q3c_bool indexes")

        if not self.table_exists():

            self.make_table(self.schema_path)

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

                run_sql_command_from_file(
                    file_path=table_q3c_path,
                    db_name=self.db_name,
                    db_user=admin_user,
                    password=admin_password,
                )
