"""
Module to define PSQL database tables using sqlalchemy
"""

import logging

from sqlalchemy import inspect

logger = logging.getLogger(__name__)


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
