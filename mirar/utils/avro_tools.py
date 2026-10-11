"""
Module with tools for reading avro schemas
"""

import logging
from pathlib import Path

from fastavro.schema import load_schema

logger = logging.getLogger(__name__)


def load_avro_schema(schema_path: Path) -> dict:
    """Load an avro record schema, which fastavro types loosely and takes as a str"""
    schema = load_schema(str(schema_path))

    if not isinstance(schema, dict):
        err = (
            f"Expected a record schema in {schema_path}, "
            f"but found a {type(schema).__name__}"
        )
        logger.error(err)
        raise ValueError(err)

    return schema
