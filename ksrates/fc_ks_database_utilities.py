"""
Plumbing shared by every table of the Ks database, which is served remotely by a long-running sqld
server (self-hosted libSQL server): address file parsing, connection, table creation/migration,
existence checks, atomic partial upserts, and the KsDatabaseError raised on any database failure.

The data-type-specific code lives in fc_ks_database_paralog.py (table "paralog_ks", one row per
species) and fc_ks_database_ortholog.py (table "ortholog_ks", one row per species pair).

Every function raises KsDatabaseError if the database cannot be reached or queried: there is no
silent fallback to the legacy TSV files once the database has been requested in the expert
configuration file (use_ks_database_server = yes).
"""
import sys
import pickle
import zlib
import logging
from contextlib import contextmanager
from dataclasses import dataclass
import libsql_client


# ======================================================================================
# (A) GENERIC PLUMBING
# ======================================================================================

class KsDatabaseError(Exception):
	"""Raised when the Ks database server can't be reached, its address file is invalid, or a query fails."""


@dataclass(frozen=True)
class TableSpec:
	"""
	Description of one database table: its name, its primary key column, and all the other columns
	with their SQLite type. Used both to create a brand new table and to add any missing columns to
	a table created by an older version of the schema (plain CREATE TABLE IF NOT EXISTS would not
	add columns to an existing table).
	"""
	name: str
	key_column: str
	columns: dict


def read_address_file(db_path):
	"""
	Parse the sqld server address file: "host:port" on the first line and, if the server was
	started with JWT auth enabled, the shared auth token on the second line.

	:param db_path: path to the address file
	:return: (address, token) tuple; token is None if the file has no second (non-blank) line
	:raises ValueError: if the file is empty or its first line isn't a "host:port" address
	"""
	with open(db_path, "r") as f:
		lines = f.read().splitlines()
	address = lines[0].strip() if lines else ""
	if not address or ":" not in address:
		raise ValueError(f"Malformed sqld server address file [{db_path}]: {address!r}")
	token = lines[1].strip() if len(lines) > 1 and lines[1].strip() else None
	return address, token


@contextmanager
def connection(db_path):
	"""
	Open a connection to the Ks database, which is served remotely by a long-running sqld server
	rather than being a local file: this lets many independent cluster jobs, possibly on different
	compute nodes, safely read/write a shared database without relying on SQLite's WAL mode over a
	network filesystem (which SQLite's own docs say doesn't work reliably across hosts).

	Used as "with connection(db_path) as client: client.execute(...)". The connection itself is opened
	lazily by libsql_client, so an unreachable server only shows up on the first query.

	:param db_path: path to a small text file, written once at startup by the sqld server job,
	                containing the server's own address as "host:port" on the first line and,
	                if the server was started with JWT auth enabled, the shared auth token on
	                the second line
	:raises KsDatabaseError: on any error reading/parsing the address file, reaching the server or
	                         running a query inside the with-block
	"""
	try:
		address, token = read_address_file(db_path)
		# ws:// (not http://) since http:// explicitly can't support transactions; note that
		# libsql-client is an archived (frozen, unmaintained) package as of writing, but its
		# execute()/close() surface is small, stable, and verified to work against a real sqld
		# instance for our usage pattern.
		client = libsql_client.create_client_sync(f"ws://{address}", auth_token=token)
	except Exception as e:
		raise KsDatabaseError(f"Could not connect to the Ks database via address file [{db_path}]: {e}") from e
	try:
		with client:
			yield client
	except KsDatabaseError:
		raise
	except Exception as e:
		raise KsDatabaseError(f"Ks database [{db_path}] could not be queried: {e}") from e


def _check_columns(spec, columns):
	"""Guards the SQL built below against typos: only columns declared in the TableSpec are allowed."""
	for column in columns:
		if column != spec.key_column and column not in spec.columns:
			raise ValueError(f"Unknown column [{column}] for table [{spec.name}]")


def _initialize_table(client, spec):
	client.execute(f"CREATE TABLE IF NOT EXISTS {spec.name} ({spec.key_column} TEXT PRIMARY KEY)")
	# Ordinary SELECT against the pragma_table_info table-valued function, rather than a
	# bare "PRAGMA table_info(...)" statement: some PRAGMA statement forms are known to be
	# rejected over sqld's remote/Hrana protocol, while this SELECT form is not.
	result = client.execute(f"SELECT name FROM pragma_table_info('{spec.name}')")
	existing_columns = {row[0] for row in result.rows}
	for column, col_type in spec.columns.items():
		if column not in existing_columns:
			try:
				client.execute(f"ALTER TABLE {spec.name} ADD COLUMN {column} {col_type}")
			except Exception as e:
				# Ignore "duplicate column" errors from concurrent ALTER TABLE calls:
				# another process may have added this column between our pragma_table_info
				# query and this ALTER statement. Any other error is real and should propagate.
				if "duplicate column" not in str(e).lower():
					raise
	# No commit needed/possible here: outside of an explicit transaction, the sqld client
	# has no .commit() method, since each execute() above is already committed by the
	# server as soon as it returns.


def initialize_database(db_path):
	"""
	Create every table of the Ks database (paralog_ks and ortholog_ks) if it doesn't exist yet, and
	add any columns missing from a database created by an older version of the schema. Idempotent.
	Called once per command so that a server is always ready for either data type.

	:param db_path: path to the sqld server address file (see connection)
	:raises KsDatabaseError: if the server is unreachable or the address file is invalid
	"""
	# Imported here because those modules themselves import this one
	from ksrates.fc_ks_database_paralog import PARALOG_TABLE
	from ksrates.fc_ks_database_ortholog import ORTHOLOG_TABLE
	with connection(db_path) as client:
		for spec in (PARALOG_TABLE, ORTHOLOG_TABLE):
			_initialize_table(client, spec)


def initialize_or_exit(db_path):
	"""
	Initialize the Ks database for a command that needs it (use_ks_database_server = yes), stopping the
	command with a clear error message if the server can't be used: once the database has been
	requested, the analysis must not silently continue on the TSV files.

	:param db_path: path to the sqld server address file (see connection)
	"""
	try:
		initialize_database(db_path)
	except KsDatabaseError as e:
		logging.error(f"{e}")
		logging.error(f'The Ks database was requested in the expert configuration file (use_ks_database_server = yes) but it cannot be used.')
		logging.error(f'Please check that the sqld server is running and that "ks_database_server_address_path" points to its address file [{db_path}], or set "use_ks_database_server = no". Exiting.')
		sys.exit(1)


def row_exists(db_path, spec, key):
	"""
	Check whether a key already has a row in a table, regardless of which columns hold data.

	:param db_path: path to the sqld server address file (see connection)
	:param spec: TableSpec of the table to look into
	:param key: primary key value
	:return: True if the row exists, False otherwise
	:raises KsDatabaseError: if the database can't be read
	"""
	with connection(db_path) as client:
		result = client.execute(f"SELECT 1 FROM {spec.name} WHERE {spec.key_column} = ? LIMIT 1", (key,))
	return len(result.rows) > 0


def existing_columns(db_path, spec, key, columns):
	"""
	Check which columns of a row hold non-null data, without transferring the data itself (columns
	can hold large blobs).

	:param db_path: path to the sqld server address file (see connection)
	:param spec: TableSpec of the table to look into
	:param key: primary key value
	:param columns: list of column names to check
	:return: dict column -> True if populated; all False if the row doesn't exist yet
	:raises KsDatabaseError: if the database can't be read
	"""
	_check_columns(spec, columns)
	select = ", ".join(f"{column} IS NOT NULL" for column in columns)
	with connection(db_path) as client:
		result = client.execute(f"SELECT {select} FROM {spec.name} WHERE {spec.key_column} = ?", (key,))
	if not result.rows:
		return {column: False for column in columns}
	return {column: bool(value) for column, value in zip(columns, result.rows[0])}


def read_columns(db_path, spec, key, columns):
	"""
	Read some columns of one row, without loading the rest of the table.

	:param db_path: path to the sqld server address file (see connection)
	:param spec: TableSpec of the table to look into
	:param key: primary key value
	:param columns: list of column names to read
	:return: tuple of values in the same order as columns, or None if the row doesn't exist
	:raises KsDatabaseError: if the database can't be read
	"""
	_check_columns(spec, columns)
	with connection(db_path) as client:
		result = client.execute(f"SELECT {', '.join(columns)} FROM {spec.name} WHERE {spec.key_column} = ?", (key,))
	return tuple(result.rows[0]) if result.rows else None


def upsert(db_path, spec, key, values, clear=()):
	"""
	Write some columns of one row in a single atomic statement (INSERT ... ON CONFLICT DO UPDATE).
	Only the columns with a non-None value are written; the others are left untouched, so a partial
	write never wipes out columns that another (possibly concurrent) run already wrote.

	:param db_path: path to the sqld server address file (see connection)
	:param spec: TableSpec of the table to write to
	:param key: primary key value
	:param values: dict column -> value (None values are ignored)
	:param clear: columns to explicitly set to NULL in the same statement (the only way to erase a value)
	:return: True if something was written, False if there was nothing to write
	:raises KsDatabaseError: if the write fails
	"""
	values = {column: value for column, value in values.items() if value is not None}
	values.update({column: None for column in clear})
	if not values:
		return False
	_check_columns(spec, values.keys())
	columns = list(values.keys())
	placeholders = ", ".join(["?"] * (len(columns) + 1))
	set_clause = ", ".join(f"{column}=excluded.{column}" for column in columns)
	with connection(db_path) as client:
		client.execute(f"""
			INSERT INTO {spec.name} ({spec.key_column}, {', '.join(columns)})
			VALUES ({placeholders})
			ON CONFLICT({spec.key_column}) DO UPDATE SET {set_clause}
		""", (key, *[values[column] for column in columns]))
	return True


def list_keys(db_path, spec, not_null_column=None):
	"""
	List the primary keys of a table, in alphabetical order.

	:param db_path: path to the sqld server address file (see connection)
	:param spec: TableSpec of the table to list
	:param not_null_column: if given, only the keys whose row has non-null data in this column
	:return: list of keys
	:raises KsDatabaseError: if the database can't be read
	"""
	where = ""
	if not_null_column is not None:
		_check_columns(spec, [not_null_column])
		where = f"WHERE {not_null_column} IS NOT NULL"
	with connection(db_path) as client:
		result = client.execute(f"SELECT {spec.key_column} FROM {spec.name} {where} ORDER BY {spec.key_column}")
	return [row[0] for row in result.rows]


def delete_keys(db_path, spec, keys):
	"""
	Delete rows from a table.

	:param db_path: path to the sqld server address file (see connection)
	:param spec: TableSpec of the table to delete from
	:param keys: list of primary keys of the rows to delete
	:raises KsDatabaseError: if the deletion fails
	"""
	if not keys:
		return
	placeholders = ", ".join(["?"] * len(keys))
	with connection(db_path) as client:
		client.execute(f"DELETE FROM {spec.name} WHERE {spec.key_column} IN ({placeholders})", tuple(keys))


def compress(obj):
	return zlib.compress(pickle.dumps(obj))


def decompress(blob):
	return pickle.loads(zlib.decompress(blob))
