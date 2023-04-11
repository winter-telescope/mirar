"""
Module for summer-specific db creation with psycopg
"""
import os

summer_schema_dir = os.path.dirname(__file__)


def get_summer_schema_path(db_name: str) -> str:
    return os.path.join(summer_schema_dir, f"{db_name}.sql")
