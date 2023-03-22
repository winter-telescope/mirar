"""
Base class for models
"""
from datetime import date
from typing import ClassVar

from pydantic import BaseModel, Extra, Field, root_validator, validator
from sqlalchemy import Insert, Select, Table, inspect, select
from sqlalchemy.orm import DeclarativeBase

from winterdrp.utils.sql import get_engine

db_name = "summertest"


class PydanticBase(BaseModel):
    """
    Base code pydantic model (no extra colunns!)
    """

    class Config:
        """
        Config
        """

        extra = Extra.forbid


class BaseDB(PydanticBase):
    """
    Base Database Table model, requiring an associated SQLalchemy table
    """

    sql_model: ClassVar[Table]

    @root_validator(pre=True)
    def sql_model_exists_check(cls, values):
        """
        Validator to ensure an sql model has been specified in the child class

        :param values: values
        :return: values
        """

        if "sql_model" not in cls.__annotations__:
            raise ValueError("sql model must be set in class")
        return values

    @validator("*")
    @classmethod
    def validate_sql(cls, value, field):
        """
        Validator to ensure that the field names of a pydantic model
        match the database table

        :param value: value
        :param field: field
        :return: value
        """
        if field.name not in cls.sql_model.__table__.columns.keys():
            err = f"Field '{field.name}' not duplicated in {cls.sql_model}"
            raise ValueError(err)
        return value

    def _insert_entry(self):
        """
        Insert the pydantic-ified data into the corresponding sql database

        :return: None
        """
        stmt = Insert(self.sql_model).values(**self.dict())

        engine = get_engine(db_name=self.sql_model.db_name)
        with engine.connect() as conn:
            res = conn.execute(stmt)
            conn.commit()

        return res

    def insert_entry(self):
        """
        Insert the pydantic-ified data into the corresponding sql database

        :return: None
        """
        self._insert_entry()

    # def exists(self, key=None) -> bool:
    #     """
    #     Checks if the pydantic-ified data exists the corresponding sql database
    #
    #     :return: bool
    #     """
    #     if key is None:
    #         key = self.sql_model().get_primary_key()
    #
    #     if key not in self.__dict__.keys():
    #         err = f"Key {key} not present in the {self.__class__} fields." \
    #               f"Either provide a different key, or consider using the exists " \
    #               f"function in {self.sql_model} instead"
    #         raise ValueError(err)
    #     print(f"Pydantic existential check")
    #     return self.sql_model().exists(value=self.__dict__[key], key=key)


class Base(DeclarativeBase):
    """
    Parent class for databases
    """

    db_name = db_name

    def get_primary_key(self) -> str:
        return inspect(self.__class__).primary_key[0].name

    def exists(self, value, key: str = None) -> bool:
        """
        Function to query a table and see whether an entry with key==value exists.
        If key is None, key will default to the table's primary key.
        Args:
            value: = Value to match
            key: str = column name

        Returns:
            boolean
        """
        engine = get_engine(db_name=self.db_name)
        if key is None:
            key = self.get_primary_key()
        else:
            columns = inspect(self.__class__).mapper.column_attrs
            column_names = [c.key for c in columns]
            if key not in column_names:
                err = f"{self.__tablename__} has no column {key}"
                raise ValueError(err)

        print(f"Checking if {key}={value} exists in {self.__tablename__}")
        stmt = select(self.__class__).where(self.__class__.__dict__[key] == value)

        return _exists(stmt, engine=engine)


def _execute(stmt, engine):
    """
    Checks if the pydantic-ified data exists the corresponding sql database

    :return: bool
    """
    with engine.connect() as conn:
        with conn.begin():
            res = conn.execute(stmt).fetchall()
    return res


def _exists(stmt: Select, engine) -> bool:
    """
    Checks if the pydantic-ified data exists the corresponding sql database

    :return: bool
    """
    return len(_execute(stmt, engine)) > 0


ra_field: float = Field(title="RA (degrees)", ge=0.0, le=360.0)
dec_field: float = Field(title="Dec (degrees)", ge=-90.0, le=90.0)
alt_field: float = Field(title="Alt (degrees)", ge=0.0, le=90.0)
az_field: float = Field(title="Az (degrees)", ge=0.0, le=360.0)
date_field: date = Field()
