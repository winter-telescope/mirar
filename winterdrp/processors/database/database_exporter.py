import astropy.io.fits
import numpy as np
import os
from winterdrp.processors.base_processor import BaseProcessor
import logging
from winterdrp.processors.database.setup import check_if_db_exists, create_db, export_to_db, create_table

logger = logging.getLogger(__name__)


class DatabaseExporter(BaseProcessor):

    base_key = "db"

    def __init__(
            self,
            db_name: str,
            db_table: str,
            schema_path: str,
            db_user: str = os.path.basename(os.environ["HOME"]),
            db_password: str = "FIXME",
            *args,
            **kwargs
    ):
        super().__init__(*args, **kwargs)
        self.db_name = db_name
        self.db_table = db_table
        self.db_user = db_user
        self.db_password = db_password

        logger.error("Here...")

        if not self.db_exists():
            self.make_db()

        self.make_table(schema_path)

    def db_exists(self):
        return check_if_db_exists(
            db_name=self.db_name,
            db_user=self.db_user,
            password=self.db_password
        )

    def make_db(self):
        create_db(
            db_name=self.db_name,
            db_user=self.db_user,
            password=self.db_password
        )

    def make_table(self, schema_path: str):
        create_table(
            schema_path,
            db_name=self.db_name,
            db_user=self.db_user,
            password=self.db_password
        )

    def _apply_to_images(
            self,
            images: list[np.ndarray],
            headers: list[astropy.io.fits.Header],
    ) -> tuple[list[np.ndarray], list[astropy.io.fits.Header]]:

        for header in headers:
            export_to_db(
                header,
                db_name=self.db_name,
                db_table=self.db_table,
                db_user=self.db_user,
                password=self.db_password
            )
        return images, headers
