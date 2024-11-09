"""
Models for the 'fields' table
"""

import time
from typing import ClassVar

from pydantic import Field
from sqlalchemy import REAL, Column, Insert, Integer
from sqlalchemy.orm import Mapped, relationship
from tqdm import tqdm
from wintertoo.data import all_winter_fields

from mirar.database.base_model import BaseDB, dec_field, ra_field
from mirar.database.engine import get_engine
from mirar.database.transactions.select import is_populated
from mirar.pipelines.winter.models.base_model import WinterBase

DEFAULT_FIELD = 999999999


class FieldsTable(WinterBase):  # pylint: disable=too-few-public-methods
    """
    Field table in database
    """

    __tablename__ = "fields"

    fieldid = Column(Integer, primary_key=True)
    ra = Column(REAL, nullable=True)
    dec = Column(REAL, nullable=True)
    ebv = Column(REAL, nullable=True)
    gall = Column(REAL, nullable=True)
    galb = Column(REAL, nullable=True)
    exposures: Mapped["ExposuresTable"] = relationship(back_populates="field")


fieldid_field: int = Field(ge=0, default=DEFAULT_FIELD)


class FieldEntry(BaseDB):
    """
    A pydantic model for a fields database entry
    """

    sql_model: ClassVar = FieldsTable
    fieldid: int = fieldid_field
    ra: float = ra_field
    dec: float = dec_field
    ebv: float = Field(ge=1)
    gall: float = ra_field
    galb: float = dec_field
    # Need to fix formatting of fields file in wintertoo before including these.


def populate_fields():
    """
    Downloads a field grid (text file) and imports it in chunks into the database

    :return: None
    """

    engine = get_engine(db_name=FieldsTable.db_name)

    if not is_populated(FieldsTable):
        chunk = 10000

        all_winter_fields["fieldid"] = all_winter_fields["ID"]
        all_winter_fields["ra"] = all_winter_fields["RA"]
        all_winter_fields["dec"] = all_winter_fields["Dec"]
        all_winter_fields["ebv"] = all_winter_fields["EBV"]
        all_winter_fields["gall"] = all_winter_fields["Gal_Long"]
        all_winter_fields["galb"] = all_winter_fields["Gal_Lat"]

        keys = list(FieldEntry.__fields__)

        idx = list(range(0, len(all_winter_fields), chunk)) + [len(all_winter_fields)]

        for k, i in tqdm(enumerate(idx[:-1]), total=len(idx) - 1):
            j = idx[k + 1]

            res = all_winter_fields[i:j]

            stmt = Insert(FieldsTable).values(
                res[keys].to_dict(orient="records"),
            )

            with engine.connect() as conn:
                with conn.begin():
                    conn.execute(stmt)

            time.sleep(1)

        stmt = Insert(FieldsTable).values(fieldid=DEFAULT_FIELD, ra=None, dec=None)

        with engine.connect() as conn:
            with conn.begin():
                conn.execute(stmt)
