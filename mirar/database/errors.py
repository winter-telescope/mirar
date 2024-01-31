"""
Database errors
"""

from mirar.errors.exceptions import ProcessorError


class DataBaseError(ProcessorError):
    """Error relating to postgres interactions"""
