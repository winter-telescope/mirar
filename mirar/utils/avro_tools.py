"""
Module with tools for reading avro schemas
"""

import logging
from pathlib import Path

from fastavro.schema import load_schema

logger = logging.getLogger(__name__)


def load_avro_schema(schema_path: Path) -> dict:
    """
    Load an avro record schema from a file

    fastavro takes the path as a str, and types the result loosely as
    str | list | dict, so this narrows it to the record dict callers expect.

    :param schema_path: Path of the avro schema file
    :return: Schema, as a dictionary
    """
    schema = load_schema(str(schema_path))

    if not isinstance(schema, dict):
        err = (
            f"Expected a record schema in {schema_path}, "
            f"but found a {type(schema).__name__}"
        )
        logger.error(err)
        raise ValueError(err)

    return schema
