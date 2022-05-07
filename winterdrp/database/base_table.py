import astropy.table
import os
import logging
import astropy.io.fits
import psycopg2

logger = logging.getLogger(__name__)


class BaseTable:

    @property
    def abbreviation(self):
        raise NotImplementedError()

    def __init__(
            self,
            schema: str,
            user: str,
            password: str,
            db_name: str
    ):
        self.schema = schema
        self.user = user
        self.password = password
        self.db_name = db_name

        self.conn = psycopg2.connect(database=self.db_name, user=self.db_user, password=self.password)
        self.conn.autocommit = True
        self.cursor = self.conn.cursor()

    @staticmethod
    def create_table(self):
        self.cursor.execute(open(self.schema, "r").read())

    def ingest_database(
            self,
            header: astropy.io.fits.Header
    ):
        raise NotImplementedError()
