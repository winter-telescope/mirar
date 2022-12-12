import logging
import os
from abc import ABC

import numpy as np

from winterdrp.processors.base_processor import BaseProcessor
from winterdrp.processors.database.postgres import (
    PG_ADMIN_PWD_KEY,
    PG_ADMIN_USER_KEY,
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
    def __init__(
        self,
        db_name: str,
        db_table: str,
        schema_path: str,
        db_user: str = os.environ.get("DB_USER"),
        db_password: str = os.environ.get("DB_PWD"),
        full_setup: bool = False,
        schema_dir: str = None,
        duplicate_protocol: str = "fail",
        q3c: bool = False,
    ):
        super().__init__()
        self.db_name = db_name
        self.db_table = db_table
        self.db_user = db_user
        self.db_password = db_password
        self.full_setup = full_setup
        self.schema_path = schema_path
        self.schema_dir = schema_dir
        self.db_check = False
        self.duplicate_protocol = duplicate_protocol
        self.q3c = q3c
        assert self.duplicate_protocol in ["fail", "ignore", "replace"]

    def db_exists(self):
        return check_if_db_exists(db_name=self.db_name)

    def make_db(self):
        create_db(db_name=self.db_name)

    def user_exists(self):
        return check_if_user_exists(self.db_user)

    def make_user(self):
        return create_new_user(new_db_user=self.db_user, new_password=self.db_password)

    def grant_privileges(self):
        return grant_privileges(self.db_name, self.db_user)

    def table_exists(self):
        return check_if_table_exists(
            db_name=self.db_name,
            db_table=self.db_table,
            db_user=self.db_user,
            password=self.db_password,
        )

    def make_table(self, schema_path: str):
        create_table(
            schema_path,
            db_name=self.db_name,
            db_user=self.db_user,
            password=self.db_password,
        )

    def check_prerequisites(
        self,
    ):
        if not self.db_check:
            self.check_database_setup()
            self.db_check = True

    def check_database_setup(self):

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
                    self.schema_dir = os.path.dirname(self.schema_path)
                    logger.warning(
                        f"Warning, full db setup requested, but no schema directory specified. \
                    Will search for schema files in {self.schema_dir} directory."
                    )
                logger.info(
                    f"Looking in {self.schema_dir} directory to search for schema files"
                )

                create_tables_from_schema(
                    self.schema_dir, self.db_name, self.db_user, self.db_password
                )

                if self.q3c:
                    admin_user = os.environ.get(PG_ADMIN_USER_KEY)
                    admin_password = os.environ.get(PG_ADMIN_PWD_KEY)
                    q3c_dir = os.path.join(self.schema_dir, "q3c")
                    q3c_indexes_file = os.path.join(q3c_dir, "q3c_indexes.sql")
                    run_sql_command_from_file(
                        file_path=q3c_indexes_file,
                        db_name=self.db_name,
                        db_user=admin_user,
                        password=admin_password,
                    )
                    logger.info(f"Created q3c indexes")

        if not self.table_exists():
            self.make_table(self.schema_path)
            if self.q3c:
                q3c_dir = os.path.join(self.schema_dir, "q3c")
                table_q3c_path = os.path.join(q3c_dir, f"q3c_{self.db_table}.sql")
                admin_user = os.environ.get(PG_ADMIN_USER_KEY)
                admin_password = os.environ.get(PG_ADMIN_PWD_KEY)
                if not os.path.exists(table_q3c_path):
                    err = f"q3c extension requested but no {table_q3c_path} file found. Please add it in."
                    logger.error(err)
                    raise DataBaseError(err)
                else:
                    run_sql_command_from_file(
                        file_path=table_q3c_path,
                        db_name=self.db_name,
                        db_user=admin_user,
                        password=admin_password,
                    )
