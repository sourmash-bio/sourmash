"""
Common utility functions for handling sqlite3 databases.
"""
import os
import sqlite3
from .logging import debug_literal


def open_sqlite_db(filename):
    """
    Is this a pre-existing sqlite3 database? Return connection object if so.

    Otherwise, return None.
    """
    # does it already exist/is it non-zero size?

    # note: sqlite3.connect creates the file if it doesn't exist, which
    # we don't want in this function.
    if not os.path.exists(filename) or os.path.getsize(filename) == 0:
        return None

    # can we connect to it?
    try:
        conn = sqlite3.connect(filename)
    except (sqlite3.OperationalError, sqlite3.DatabaseError):
        debug_literal("open_sqlite_db: cannot connect.")
        return None

    # check for the 'sourmash_internal' table.
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


def add_sourmash_internal(cursor, use_type, version):
    cursor.execute("""
    CREATE TABLE IF NOT EXISTS sourmash_internal (
       key TEXT,
       value TEXT
    )
    """)

    cursor.execute('SELECT DISTINCT key, value FROM sourmash_internal')
    d = dict(cursor)

    val = d.get(use_type)
    if val is not None:
        # do version compatibility foo here?
        if version != val:
            raise Exception(f"sqlite problem: for {use_type}, want version {version}, got version {val}")
    else:
        # @CTB supply unique constraints?
        cursor.execute("""
        INSERT INTO sourmash_internal (key, value) VALUES (?, ?)
        """, (use_type, version))
