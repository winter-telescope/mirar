"""
Models for the 'field' table
"""
import time
import urllib.request
from typing import ClassVar

import pandas as pd
from pydantic import Field
from sqlalchemy import REAL, Column, Insert, Integer, Select
from tqdm import tqdm

from winterdrp.pipelines.summer.models.basemodel import Base, BaseDB, _exists, dec, ra
from winterdrp.utils.sql import get_engine

DEFAULT_FIELD = 999999999


class FieldsTable(Base):  # pylint: disable=too-few-public-methods
    """
    Field table in database
    """

    __tablename__ = "fields"

    fieldid = Column(Integer, primary_key=True)
    ra = Column(REAL, nullable=True)
    dec = Column(REAL, nullable=True)


fieldid_field: int = Field(ge=0, default=DEFAULT_FIELD)


class Fields(BaseDB):
    """
    A pydantic model for a fields database entry
    """

    sql_model: ClassVar = FieldsTable
    fieldid: int = fieldid_field
    ra: float = ra
    dec: float = dec


_SUMMER_FIELDS_URL = (
    "https://github.com/winter-telescope/wintertoo/raw/"
    "main/wintertoo/data/summer_fields.txt"
)


def populate_fields(url=_SUMMER_FIELDS_URL):
    """
    Downloads a field grid (text file) and imports it in chunks into the database

    :param url: url of grid
    :return: None
    """

    if not _exists(Select(FieldsTable)):
        engine = get_engine()

        with urllib.request.urlopen(url) as url_s:
            full_res = pd.read_csv(url_s, sep=r"\s+")

        chunk = 10000

        full_res["fieldid"] = full_res["#ID"]
        full_res["ra"] = full_res["RA"]
        full_res["dec"] = full_res["Dec"]

        keys = list(Fields.__fields__)

        idx = list(range(0, len(full_res), chunk)) + [len(full_res)]

        for k, i in tqdm(enumerate(idx[:-1]), total=len(idx) - 1):
            j = idx[k + 1]

            res = full_res[i:j]

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
