"""
Base class for models
"""
from typing import ClassVar

from pydantic import BaseModel, Extra, root_validator, validator
from sqlalchemy import Insert, Table
from sqlalchemy.orm import declarative_base

from winterdrp.utils.sql import get_engine


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

    def insert_entry(self):
        """
        Insert the pydantic-ified data into the corresponding sql database

        :return: None
        """
        engine = get_engine()
        with engine.connect() as conn:
            with conn.begin():
                stmt = Insert(self.sql_model).values(**self.dict())
                conn.execute(stmt)


Base = declarative_base()
