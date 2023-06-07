"""
Central module for interactions with postgres databases
"""
from mirar.processors.database.base_database_processor import BaseDatabaseProcessor
from mirar.processors.database.postgres import (
    POSTGRES_DUPLICATE_PROTOCOLS,
    PostgresAdmin,
    PostgresUser,
)
