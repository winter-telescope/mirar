"""
Module to export data to a database table
"""

import logging
from typing import Type

from pydantic import ValidationError
from sqlalchemy import inspect

from mirar.data import DataBlock
from mirar.database.base_model import BaseDB
from mirar.database.constants import POSTGRES_DUPLICATE_PROTOCOLS
from mirar.database.errors import DataBaseError

logger = logging.getLogger(__name__)


def export_to_db(
    value_dict: dict | DataBlock,
    db_table: Type[BaseDB],
    duplicate_protocol: str = "fail",
) -> tuple[list, list]:
    """
    Export a list of fields in value dict to a batabase table

    :param value_dict: dictionary/DataBlock/other dictonary-like object to export
    :param db_table: table of DB to export to
    :param duplicate_protocol: protocol for handling duplicates,
        in "fail"/"ignore"/"replace"
    :return:
    """

    assert duplicate_protocol in POSTGRES_DUPLICATE_PROTOCOLS

    column_names = [x for x in db_table.__dict__["__annotations__"] if x != "sql_model"]

    column_dict = {}
    for column in column_names:
        column_dict[column] = value_dict[column]

    try:
        new = db_table(**column_dict)
    except ValidationError as err:
        logger.error(err)
        raise DataBaseError from err

    db_name = new.sql_model.db_name
    primary_key = inspect(db_table.sql_model).primary_key[0]

    sequence_key_names, sequence_values = [], []
    try:
        sequence_key_names, sequence_values = new.insert_entry()

    except IntegrityError as exc:
        if not isinstance(exc.orig, errors.UniqueViolation):
            raise exc

        if duplicate_protocol == "fail":
            err = (
                f"Duplicate error, entry with {column_dict} "
                f"already exists in {db_name}."
            )
            logger.error(err)
            raise errors.UniqueViolation from exc

        if duplicate_protocol == "ignore":
            logger.debug(
                f"Found duplicate entry in {db_name} - "
                f"{str(exc)}."
                f"Ignoring, no new entry made."
            )
            primary_key_val = value_dict[primary_key.name]
            sequence_keys = new.get_sequence_keys()
            sequence_key_names = [k.name for k in sequence_keys]
            sequence_values = []
            if len(sequence_keys) > 0:
                ret = new.sql_model().select_query(
                    compare_values=[primary_key_val],
                    compare_keys=[primary_key.name],
                    select_keys=sequence_key_names,
                )
                sequence_values = [x[0] for x in ret]

        if duplicate_protocol == "replace":
            logger.debug(f"Conflict at {exc.orig.diag.constraint_name}")
            logger.debug(
                f"Found duplicate entry in {db_name} - "
                f"{str(exc)}."
                f"Replacing with a new entry."
            )
            primary_key_val = value_dict[primary_key.name]
            sequence_key_names, sequence_values = new.update_entry(
                primary_key_val=primary_key_val
            )

    return sequence_key_names, sequence_values
