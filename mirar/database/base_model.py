"""
Module to define PSQL database tables using sqlalchemy
"""

import logging
from datetime import date
from typing import Any, ClassVar, Type

import pandas as pd
from psycopg import errors
from pydantic import (
    BaseModel,
    ConfigDict,
    Field,
    FieldValidationInfo,
    field_validator,
    model_validator,
)
from sqlalchemy import Column, Table, inspect
from sqlalchemy.exc import IntegrityError

from mirar.database.constants import POSTGRES_DUPLICATE_PROTOCOLS
from mirar.database.constraints import DBQueryConstraints
from mirar.database.transactions import select_from_table
from mirar.database.transactions.insert import _insert_in_table
from mirar.database.transactions.update import _update_database_entry

logger = logging.getLogger(__name__)


class PydanticBase(BaseModel):
    """
    Base code pydantic model (no extra colunns!)
    """

    model_config = ConfigDict(extra="ignore")


class BaseDB(PydanticBase):
    """
    Base Database Table model, requiring an associated SQLalchemy table
    """

    sql_model: ClassVar[Type[Table]]

    @model_validator(mode="before")
    @classmethod
    def sql_model_exists_check(cls, data: Any) -> Any:
        """
        Validator to ensure an sql model has been specified in the child class
        :param data: data to validate
        :return: data
        """

        if "sql_model" not in cls.__annotations__:
            raise ValueError("sql model must be set in class")
        return data

    @field_validator("*")
    @classmethod
    def validate_sql(cls, value: Any, info: FieldValidationInfo) -> Any:
        """
        Validator to ensure that the field names of a pydantic model
        match the database table

        :param value: value
        :param info: field info
        :return: value
        """
        if info.field_name not in cls.sql_model.__table__.columns.keys():
            err = f"Field '{info.field_name}' not duplicated in {cls.sql_model}"
            raise ValueError(err)
        return value

    def _insert_entry(
        self,
        returning_key_names: str | list[str] = None,
        duplicate_protocol: str = "ignore",
    ) -> pd.DataFrame:
        """
        Insert the pydantic-ified data into the corresponding sql database

        :param returning_key_names: names of keys to return
        :param duplicate_protocol: protocol to follow if duplicate entry is found
        :return: sequence_key dataframe
        """

        assert duplicate_protocol in POSTGRES_DUPLICATE_PROTOCOLS

        if returning_key_names is None:
            returning_key_names = self.get_primary_key()

        if not isinstance(returning_key_names, list):
            returning_key_names = [returning_key_names]

        try:
            res = _insert_in_table(
                new_entry=self.model_dump(),
                sql_table=self.sql_model,
                returning_keys=returning_key_names,
            )
        except IntegrityError as exc:
            if not isinstance(exc.orig, errors.UniqueViolation):
                raise exc

            db_name = self.sql_model.db_name

            if duplicate_protocol == "fail":
                err = (
                    f"Duplicate error, entry with {self.model_dump()} "
                    f"already exists in {self.sql_model.name}."
                )
                logger.error(err)
                raise errors.UniqueViolation from exc

            if duplicate_protocol == "replace":
                logger.debug(f"Conflict at {exc.orig.diag.constraint_name}")
                logger.debug(
                    f"Found duplicate entry in {db_name} - "
                    f"{str(exc)}."
                    f"Replacing with a new entry."
                )
                self.update_entry()

            else:
                logger.debug(
                    f"Found duplicate entry in - "
                    f"{str(exc)}."
                    f"Ignoring, no new entry made."
                )

            present_unique_keys = self.get_available_unique_keys()

            assert len(present_unique_keys) > 0

            res = None

            for key in present_unique_keys:
                constraint = DBQueryConstraints(
                    columns=[key.name],
                    accepted_values=[self.model_dump()[key.name]],
                )
                new_res = select_from_table(
                    sql_table=self.sql_model,
                    db_constraints=constraint,
                    output_columns=returning_key_names,
                )

                if len(new_res) > 0:
                    if res is None:
                        res = new_res
                    elif not new_res.equals(res):
                        raise ValueError(
                            f"Multiple matches found: {new_res} and {res}"
                        ) from exc

            if res is None:
                raise ValueError("No results found") from exc

        return res

    def get_primary_key(self) -> str:
        """
        Get the primary key of the table

        :return: primary key
        """
        primary_key = inspect(self.sql_model).primary_key[0]
        return primary_key.name

    def get_unique_keys(self) -> list[Column]:
        """
        Get the unique key of the table

        :return: unique keys
        """
        cols = [x for x in self.sql_model.__table__.columns if x.unique]
        return cols

    def get_available_unique_keys(self) -> list[Column]:
        """
        Get the unique keys of the table which are present in the data

        :return: unique keys
        """
        return [x for x in self.get_unique_keys() if x.name in self.model_fields]

    def insert_entry(
        self, returning_key_names: str | list[str] | None = None
    ) -> pd.DataFrame:
        """
        Insert the pydantic-ified data into the corresponding sql database

        :param returning_key_names: names of the keys to return
        :return: dataframe of the sequence keys
        """
        result = self._insert_entry(returning_key_names=returning_key_names)
        logger.debug(f"Return result {result}")
        return result

    def _update_entry(self, update_key_names: list[str] | str | None = None):
        """
        Update database entry

        :param update_key_names: names of keys to be updates,
            if None, will update all keys
        :return: None
        """

        available_unique_keys = self.get_available_unique_keys()

        assert len(available_unique_keys) > 0

        constraints = DBQueryConstraints(
            columns=[x.name for x in available_unique_keys],
            accepted_values=[self.model_dump()[x.name] for x in available_unique_keys],
        )

        full_dict = self.model_dump()

        update_dict = {key: full_dict[key] for key in update_key_names}

        _update_database_entry(
            update_dict=update_dict,
            sql_table=self.sql_model,
            db_constraints=constraints,
        )

    def update_entry(self, update_keys=None):
        """
        Wrapper to update database entry. Users should override this function.

        :param update_keys: keys to update
        :return: None
        """
        self._update_entry(update_keys)

    @classmethod
    def _exists(cls, values, keys: str | list = None) -> bool:
        """
        Function to query a table and see whether an entry with key==value exists.
        If key is None, key will default to the table's primary key.

        :param values: values to query
        :param keys: keys to query
        :return: True if entry exists, False otherwise
        """
        db_constraints = DBQueryConstraints(
            columns=keys, accepted_values=values, comparison_types="="
        )

        match = select_from_table(
            db_constraints=db_constraints,
            sql_table=cls.sql_model,
        )
        return len(match) > 0


ra_field: float = Field(title="RA (degrees)", ge=0.0, le=360.0)
dec_field: float = Field(title="Dec (degrees)", ge=-90.0, le=90.0)
alt_field: float = Field(title="Alt (degrees)", ge=0.0, le=90.0)
az_field: float = Field(title="Az (degrees)", ge=0.0, le=360.0)
date_field: date = Field()
