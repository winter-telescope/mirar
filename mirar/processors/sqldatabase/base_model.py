"""
Module to define PSQL database tables using sqlalchemy
"""
import logging
from datetime import date
from typing import ClassVar

import numpy as np
from pydantic import BaseModel, Extra, Field, root_validator, validator
from sqlalchemy import Insert, Select, Table, Update, and_, inspect, select

from mirar.processors.sqldatabase.postgres_utils import (
    get_sequence_key_names_from_table,
)
from mirar.utils.sql import get_engine

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

    def _insert_entry(self, returning_key_names: list = None) -> tuple:
        """
        Insert the pydantic-ified data into the corresponding sql database
        :return: sequence_key_names, sequence_key_values of the updated entry
        """
        if returning_key_names is None:
            returning_keys = self.get_sequence_keys()
        else:
            if not isinstance(returning_key_names, list):
                returning_key_names = [returning_key_names]

            returning_keys = self.get_table_keys_from_names(returning_key_names)

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

    def insert_entry(self, returning_key_names: str | list = None):
        """
        Insert the pydantic-ified data into the corresponding sql database
        :return: sequence_key_names, sequence_key_values of the updated entry
        """
        result = self._insert_entry(returning_key_names=returning_key_names)
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

        returning_keys = self.get_sequence_keys()

        update_dict = {}
        for x, key in enumerate(update_keys):
            update_dict[key.key] = update_vals[x]

        stmt = (
            Update(self.sql_model)
            .values(**update_dict)
            .where(primary_key == primary_key_val)
            .returning(*returning_keys)
        )

        logger.debug(stmt)
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
        """
        Function to get primary key of table
        Returns:
        primary key
        """
        return inspect(self.__class__).primary_key[0].name

    def construct_select_statement(
        self,
        compare_values: list,
        select_keys: str | list[str] = None,
        compare_keys: str | list[str] = None,
        comparators: str | list = None,
    ):
        """
        Function to construct a select statement
        Args:
            compare_values: Values to compare
            select_keys: Keys to select. If None, defaults to primary key
            compare_keys: Keys to compare. If None, defaults to primary key
            comparators: Comparators : '__eq__', '__gt__', '__lt__', etc.
            If not specified, defaults to '__eq__
        Returns:
        """

        if compare_keys is None:
            compare_keys = [self.get_primary_key()]
        else:
            if not isinstance(compare_keys, list):
                compare_keys = [compare_keys]
            if not isinstance(compare_values, list):
                compare_values = [compare_values]
            assert len(compare_keys) == len(compare_values)
            columns = inspect(self.__class__).mapper.column_attrs
            column_names = [c.key for c in columns]
            for key in compare_keys:
                if key not in column_names:
                    err = f"{self.__tablename__} has no column {key}"
                    raise ValueError(err)

        if comparators is None:
            comparators = ["__eq__"] * len(compare_keys)

        if not isinstance(comparators, list):
            comparators = [comparators]

        assert len(compare_keys) == len(comparators)

        logger.debug(
            f"Checking if {compare_keys} {comparators} {compare_values} "
            f"exists in {self.__tablename__}"
        )

        if select_keys is None:
            select_items = [self.__class__]
        else:
            if not isinstance(select_keys, list):
                select_keys = [select_keys]
            select_items = [getattr(self.__class__, k) for k in select_keys]
        stmt = select(*select_items).where(
            and_(
                *[
                    getattr(self.__class__, key).__getattribute__(comparator)(value)
                    for key, comparator, value in zip(
                        compare_keys, comparators, compare_values
                    )
                ]
            )
        )

        return stmt

    def select_query(
        self,
        compare_values: list,
        select_keys: str | list[str] = None,
        compare_keys: str | list[str] = None,
        comparators: str | list[str] = None,
    ):
        """
        Run a select query
        Args:
            compare_values: Values to compare
            select_keys: Keys to select. If None, defaults to primary key
            compare_keys: Keys to compare. If None, defaults to primary key
            comparators: Comparators : '__eq__', '__gt__', '__lt__', etc.
            If None, defaults to '__eq__

        Returns:
            result of select query
        """
        engine = get_engine(db_name=self.db_name)
        stmt = self.construct_select_statement(
            compare_values=compare_values,
            select_keys=select_keys,
            compare_keys=compare_keys,
            comparators=comparators,
        )
        with engine.connect() as conn:
            res = conn.execute(stmt)
            conn.commit()
        return res.fetchall()

    def exists(self, values, keys: str | list = None) -> bool:
        """
        Function to query a table and see whether an entry with key==value exists.
        If key is None, key will default to the table's primary key.
        Args:
            values: = Value to match
            keys: str = column name

        Returns:
            boolean
        """
        if not isinstance(keys, list):
            keys = [keys]

        if not isinstance(values, list):
            values = [values]
        assert len(values) == len(keys)
        engine = get_engine(db_name=self.db_name)

        logger.debug(f"Checking if {keys}=={values} exists in {self.__tablename__}")
        stmt = self.construct_select_statement(
            compare_values=values, compare_keys=keys, comparators=["__eq__"] * len(keys)
        )
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
