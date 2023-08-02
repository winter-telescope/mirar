"""
Central module for all DB transaction types.
"""
from mirar.database.transactions.insert import insert_in_table
from mirar.database.transactions.select import select_from_table
from mirar.database.transactions.update import update_database_entry
