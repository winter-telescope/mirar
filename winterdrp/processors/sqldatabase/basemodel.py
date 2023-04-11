"""
Module to define PSQL database tables using sqlalchemy
"""
import logging
from typing import ClassVar

import numpy as np
from pydantic import BaseModel, Extra, root_validator, validator
from sqlalchemy import Insert, Select, Table, Update, inspect, select

from winterdrp.processors.sqldatabase.postgres import get_sequence_key_names_from_table
from winterdrp.utils.sql import get_engine

logger = logging.getLogger(__name__)


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

    def _insert_entry(self, returning_keys: list = None) -> tuple:
        """
        Insert the pydantic-ified data into the corresponding sql database

        :return: sequence_key_names, sequence_key_values of the updated entry
        """
        if returning_keys is None:
            returning_keys = self.get_sequence_keys()
        logger.debug(f"Returning keys {returning_keys}")
        stmt = Insert(self.sql_model).values(**self.dict()).returning(*returning_keys)
        engine = get_engine(db_name=self.sql_model.db_name)
        with engine.connect() as conn:
            res = conn.execute(stmt)
            conn.commit()

        if len(returning_keys) == 0:
            returning_values = []
        else:
            returning_values = res.fetchall()[0]
        logger.debug(returning_values)

        returning_key_names = [x.key for x in returning_keys]
        assert len(returning_values) == len(returning_key_names)
        return returning_key_names, returning_values

    def insert_entry(self, returning_keys=None):
        """
        Insert the pydantic-ified data into the corresponding sql database

        :return: sequence_key_names, sequence_key_values of the updated entry
        """
        result = self._insert_entry(returning_keys=returning_keys)
        logger.debug(f"Return result {result}")
        return result

    def get_sequence_keys(self) -> list:
        """
        Get sequence column instances of the sql_model table
        Returns:
        list of sequence keys
        """
        sequence_key_names = get_sequence_key_names_from_table(
            db_name=self.sql_model.db_name, db_table=self.sql_model.__tablename__
        )
        sequence_keys = self.get_table_keys_from_names(key_names=sequence_key_names)
        return sequence_keys

    def get_table_keys_from_names(self, key_names: list | np.ndarray) -> list:
        """
        Function to get column instances from names of sql table columns
        Args:
            key_names:

        Returns:
            list
        """
        return [self.sql_model.__dict__[x] for x in key_names]

    def _update_entry(self, primary_key_val, update_key_names=None) -> tuple:
        """
        Update database entry
        Args:
            primary_key_val: value of primary key to find the db entry
            update_key_names: names of keys to be updates, if None, will update all keys

        Returns:
            sequence_key_names, sequence_key_values of the updated entry
        """
        primary_key = inspect(self.sql_model).primary_key[0]

        if update_key_names is None:
            update_key_names = [x for x in self.dict() if x != primary_key.name]

        update_keys = self.get_table_keys_from_names(update_key_names)
        update_vals = [self.dict()[x] for x in update_key_names]
        logger.debug(f"{update_keys}, {update_vals}")

        returning_keys = self.get_sequence_keys()

        update_dict = {}
        for x in range(len(update_keys)):
            update_dict[update_keys[x].key] = update_vals[x]

        logger.debug(update_dict)
        stmt = (
            Update(self.sql_model)
            .values(**update_dict)
            .where(primary_key == primary_key_val)
            .returning(*returning_keys)
        )

        engine = get_engine(db_name=self.sql_model.db_name)
        with engine.connect() as conn:
            res = conn.execute(stmt)
            conn.commit()

        if len(returning_keys) == 0:
            returning_values = []
        else:
            returning_values = res.fetchall()[0]
        returning_key_names = [x.key for x in returning_keys]
        logger.debug(f"{returning_key_names}, {returning_values}")
        assert len(returning_values) == len(returning_key_names)
        return returning_key_names, returning_values

    def update_entry(self, primary_key_val, update_keys=None) -> tuple:
        """
        Wrapper to update database entry. Users should override this function.
        Args:
            primary_key_val:
            update_keys:

        Returns:
            sequence_key_names, sequence_key_values of the updated entry
        """
        return self._update_entry(primary_key_val, update_keys)


class BaseTable:
    """
    Parent class for database tables. Tables should inherit from this
    and DeclarativeBase.
    """

    @property
    def __tablename__(self):
        raise NotImplementedError

    @property
    def db_name(self):
        """
        Name of the database.

        :return: None
        """
        raise NotImplementedError

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

        logger.debug(f"Checking if {key}={value} exists in {self.__tablename__}")
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
