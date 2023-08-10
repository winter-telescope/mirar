"""
Module to select database entries
"""
import pandas as pd
from sqlalchemy import Select, text

from mirar.database.base_table import BaseTable
from mirar.database.constraints import DBQueryConstraints
from mirar.database.engine import get_engine


def run_select(
    query: Select,
    sql_table: BaseTable,
    columns: list[str] | None = None,
) -> pd.DataFrame:
    """
    Run a select query

    :param query: select query to run
    :param sql_table: table to run query on
    :param columns: columns to output (default: all)
    :return: results
    """

    engine = get_engine(db_name=sql_table.db_name)

    with engine.connect() as conn:
        res = pd.read_sql(query, conn, columns=columns)

    if (columns is not None) & len(res) > 0:
        res = res[columns]

    return res


def select_from_table(
    db_constraints: DBQueryConstraints,
    sql_table: BaseTable,
    output_columns: list[str] | None = None,
    max_num_results: int | None = None,
) -> pd.DataFrame:
    """
    Select database entries

    :param db_constraints: database query constraints
    :param sql_table: database SQL table
    :param output_columns: columns to output (default: all)
    :param max_num_results: maximum number of results to return (default: all)
    :return: results
    """
    query = Select(sql_table).where(text(db_constraints.parse_constraints()))
    if max_num_results is not None:
        query = query.limit(max_num_results)
    res = run_select(
        query=query,
        sql_table=sql_table,
        columns=output_columns,
    )

    return res


def check_table_exists(
    sql_table: BaseTable,
) -> bool:
    """
    Check if a table exists

    :param sql_table: database SQL table
    :return: True if table exists, False otherwise
    """

    engine = get_engine(db_name=sql_table.db_name)

    with engine.connect() as conn:
        res = conn.execute(
            text(
                f"""
                SELECT EXISTS (
                    SELECT FROM information_schema.tables
                    WHERE  table_schema = 'public'
                    AND    table_name   = '{sql_table.__tablename__}'
                );
                """
            )
        )

    return res.fetchone()[0]


def is_populated(
    sql_table: BaseTable,
) -> bool:
    """
    Function to check if a table is populated

    :param sql_table: database SQL table
    :return: boolean
    """

    engine = get_engine(db_name=sql_table.db_name)

    with engine.connect() as conn:
        res = conn.execute(
            text(
                f"""
                SELECT EXISTS (
                    SELECT FROM {sql_table.__tablename__}
                );
                """
            )
        )

    return res.fetchone()[0]
