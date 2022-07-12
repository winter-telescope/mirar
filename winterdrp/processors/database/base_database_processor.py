import numpy as np
import os
from winterdrp.processors.base_processor import BaseProcessor
import logging
from abc import ABC
from winterdrp.processors.database.postgres import check_if_db_exists, check_if_user_exists, check_if_table_exists,\
    create_db, create_table, create_new_user, grant_privileges, create_tables_from_schema, DataBaseError, \
    default_db_user


logger = logging.getLogger(__name__)


class BaseDatabaseProcessor(BaseProcessor, ABC):

    def __init__(
            self,
            db_name: str,
            db_table: str,
            schema_path: str,
            db_user: str = os.environ.get('PG_DEFAULT_USER'),
            db_password: str = os.environ.get('PG_DEFAULT_PWD'),
            full_setup: bool = False,
            schema_dir: str = None,
            *args,
            **kwargs
    ):
        super().__init__(*args, **kwargs)
        self.db_name = db_name
        self.db_table = db_table
        self.db_user = db_user
        self.db_password = db_password
        self.full_setup = full_setup
        self.schema_path = schema_path
        self.schema_dir = schema_dir

        # logger.error("Here...")

        if np.logical_and(self.db_exists(), np.invert(self.user_exists())):
            err = 'Database exists but user does not exist'
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
                    logger.warning(f'Warning, full db setup requested, but no schema directory specified. \
                    Will search for schema files in {self.schema_dir} directory.')
                logger.info(f'Looking in {self.schema_dir} directory to search for schema files')

                create_tables_from_schema(self.schema_dir, self.db_name, self.db_user, self.db_password)

        if not self.table_exists():
            self.make_table(schema_path)

    def db_exists(self):
        return check_if_db_exists(
            db_name=self.db_name
        )

    def make_db(self):
        create_db(
            db_name=self.db_name,
            db_user=self.db_user,
            password=self.db_password
        )

    def user_exists(self):
        return check_if_user_exists(self.db_user)

    @staticmethod
    def make_user():
        return create_new_user

    def grant_privileges(self):
        return grant_privileges(self.db_name, self.db_user)

    def table_exists(self):
        return check_if_table_exists(db_name=self.db_name,
                                     db_table=self.db_table,
                                     db_user=self.db_user,
                                     password=self.db_password)

    def make_table(self, schema_path: str):
        create_table(
            schema_path,
            db_name=self.db_name,
            db_user=self.db_user,
            password=self.db_password
        )


