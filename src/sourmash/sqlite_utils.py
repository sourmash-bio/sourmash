"""
Common utility functions for handling sqlite3 databases.
"""
import os
import sqlite3
from .logging import debug_literal


def open_sqlite_db(filename):
    """
    Is this a pre-existing sqlite3 database? Return connection object if so.
    """
    # 
    if not os.path.exists(filename) or os.path.getsize(filename) == 0:
        return None

    # can we connect to it?
    try:
        conn = sqlite3.connect(filename)
    except (sqlite3.OperationalError, sqlite3.DatabaseError):
        debug_literal("open_sqlite_db: cannot connect.")
        return None

    # grab a schema dump.
    cursor = conn.cursor()
    try:
        cursor.execute('SELECT DISTINCT key, value FROM sourmash_internal')
    except (sqlite3.OperationalError, sqlite3.DatabaseError):
        debug_literal("open_sqlite_db: cannot read sourmash_internal.")

        # is this a taxonomy DB?
        try:
            cursor.execute('SELECT * FROM taxonomy LIMIT 1')
        except (sqlite3.OperationalError, sqlite3.DatabaseError):
            debug_literal("open_sqlite_db: cannot read 'taxonomy', either.")
            return None

    return conn
