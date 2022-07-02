from winterdrp.errors import ProcessorError
import os
import logging

logger = logging.getLogger(__name__)
pg_default_user = os.getenv('PG_DEFAULT_USER')
if pg_default_user is None:
    err = "No default psql user for database (postgres) specified. Run 'export PG_DEFAULT_USER=user (usually home or " \
          "pgadmin)' "
    logger.error(err)
    raise ValueError(err)

pg_default_pwd = os.getenv('PG_DEFAULT_PWD')
if pg_default_user is None:
    err = "No password specified for the default psql user. Run 'export PG_DEFAULT_PWD=password'"
    logger.error(err)
    raise ValueError(err)


class DataBaseError(ProcessorError):
    pass
