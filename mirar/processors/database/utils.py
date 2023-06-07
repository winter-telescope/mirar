"""
Module containing postgres util functions
"""
from pathlib import Path

import numpy as np


def get_column_names_from_schema(schema_file_path: str | Path) -> list[str]:
    """
    Get column names from a schema file

    :param schema_file_path: file to read
    :return: list of columns
    """
    with open(schema_file_path, "r", encoding="utf8") as schema_file:
        dat = schema_file.read()
    dat = dat.split(");")[0]
    dat = dat.split("\n")[1:-1]
    pkstrip = [x.strip(",").split("PRIMARY KEY")[0].strip() for x in dat]
    fkstrip = [x.strip(",").split("FOREIGN KEY")[0].strip() for x in pkstrip]
    colnames = [x.split(" ")[0].strip('"') for x in fkstrip]
    return colnames


def get_foreign_tables_list(schema_files: list[str]) -> np.ndarray:
    """
    Returns a list of foreign tables

    :param schema_files: List of schema files to read
    :return: Returns list of foreign keys in schema
    """
    foreign_tables_list = []
    for schema_file_path in schema_files:
        table_names = []
        with open(schema_file_path, "r", encoding="utf8") as schema_file:
            schema = schema_file.read()
        if "FOREIGN KEY" in schema:
            schema = schema.replace("\n", "")
            schema = schema.replace("\t", "")
            schema_split = np.array(schema.split(","))
            fk_rows = np.array(["FOREIGN KEY" in x for x in schema_split])
            for row in schema_split[fk_rows]:
                words = np.array(row.split(" "))
                refmask = np.array(["REFERENCES" in x for x in words])
                idx = np.where(refmask)[0][0] + 1
                tablename = words[idx].split("(")[0]
                table_names.append(tablename)
        foreign_tables_list.append(np.array(table_names))
    return np.array(foreign_tables_list)


def get_ordered_schema_list(schema_files: list[str]) -> list[str]:
    """
    Returns an ordered list of schema, ensuring connected database tables
    are created in the right order

    :param schema_files: List of schema files to read
    :return: Returns list of foreign keys in schema
    """
    foreign_tables_list = get_foreign_tables_list(schema_files)
    ordered_schema_list = []
    tables_created = []
    schema_table_names = [x.split("/")[-1].split(".sql")[0] for x in schema_files]
    while len(tables_created) < len(schema_files):
        for ind, schema_file in enumerate(schema_files):
            table_name = schema_table_names[ind]
            if table_name in tables_created:
                pass
            else:
                foreign_tables = foreign_tables_list[ind]
                if len(foreign_tables) == 0:
                    ordered_schema_list.append(schema_file)
                    tables_created.append(table_name)
                else:
                    if np.all(np.isin(foreign_tables, tables_created)):
                        ordered_schema_list.append(schema_file)
                        tables_created.append(table_name)

    return ordered_schema_list
