"""
Shared Ks database, served remotely by a long-running sqld server (self-hosted libSQL server),
holding both paralog Ks data (table "paralog_ks", one row per species) and ortholog Ks data (table
"ortholog_ks", one row per species pair).

The module has three sections:
  (A) generic plumbing, shared by every table: connection, table creation/migration, existence
      checks, atomic partial upserts, the KsDatabaseError raised on any database failure;
  (B) paralog-specific serialization/extraction (wgd TSV files, i-ADHoRe output files);
  (C) ortholog-specific serialization/extraction (peak + Ks list per species pair).

Every function raises KsDatabaseError if the database cannot be reached or queried: there is no
silent fallback to the legacy TSV files once the database has been requested in the expert
configuration file (use_ks_database = yes).
"""
import os
import sys
import pickle
import zlib
import logging
from contextlib import contextmanager
from dataclasses import dataclass
import libsql_client
from pandas import read_csv, DataFrame
import ksrates.fc_check_input as fcCheck


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
	with connection(db_path) as client:
		for spec in (PARALOG_TABLE, ORTHOLOG_TABLE):
			_initialize_table(client, spec)


def initialize_or_exit(db_path):
	"""
	Initialize the Ks database for a command that needs it (use_ks_database = yes), stopping the
	command with a clear error message if the server can't be used: once the database has been
	requested, the analysis must not silently continue on the TSV files.

	:param db_path: path to the sqld server address file (see connection)
	"""
	try:
		initialize_database(db_path)
	except KsDatabaseError as e:
		logging.error(f"{e}")
		logging.error(f'The Ks database was requested in the expert configuration file (use_ks_database = yes) but it cannot be used.')
		logging.error(f'Please check that the sqld server is running and that "ks_database_server_address_path" points to its address file [{db_path}], or set "use_ks_database = no". Exiting.')
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


def upsert(db_path, spec, key, values):
	"""
	Write some columns of one row in a single atomic statement (INSERT ... ON CONFLICT DO UPDATE).
	Only the columns with a non-None value are written; the others are left untouched, so a partial
	write never wipes out columns that another (possibly concurrent) run already wrote.

	:param db_path: path to the sqld server address file (see connection)
	:param spec: TableSpec of the table to write to
	:param key: primary key value
	:param values: dict column -> value (None values are ignored)
	:return: True if something was written, False if there was nothing to write
	:raises KsDatabaseError: if the write fails
	"""
	values = {column: value for column, value in values.items() if value is not None}
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


def _compress(obj):
	return zlib.compress(pickle.dumps(obj))


def _decompress(blob):
	return pickle.loads(zlib.decompress(blob))


# ======================================================================================
# (B) PARALOG Ks DATA (table paralog_ks, one row per species, key = latin name)
# ======================================================================================

# Columns needed to (1) recompute per-family/per-node weights later (via filter_compute_weights)
# for an arbitrary Ks/alignment-quality cutoff, and (2) identify which gene pair each Ks value
# belongs to (Paralog1, Paralog2), without having to re-read the original TSV file.
_REQUIRED_COLUMNS = ['Family', 'Node', 'Ks', 'AlignmentCoverage', 'AlignmentIdentity', 'AlignmentLength', 'Paralog1', 'Paralog2']

# Maps the extract_paralog_ks_from_tsv()/ks_data_dict keys to the SQLite column names
_ANALYSIS_TYPE_TO_COLUMN = {"paranome": "paranome", "anchors": "anchors", "recret": "reciprocally_retained"}

# i-ADHoRe output files needed by anchor clustering (cluster_anchor_ks.py), stored as raw text
# rather than parsed: the existing line-based parsers in fc_cluster_anchors.py can read from an
# io.StringIO wrapping this text exactly as they read from a real file, so no parser rewrite is
# needed to support the database as an alternative source.
IADHORE_FILES = ['anchorpoints', 'multiplicons', 'segments', 'list_elements', 'multiplicon_pairs']
IADHORE_COLUMN = {key: f"{key}_txt" for key in IADHORE_FILES}

PARALOG_TABLE = TableSpec(
	name="paralog_ks",
	key_column="latin_name",
	columns={
		"paranome": "BLOB",
		"anchors": "BLOB",
		"reciprocally_retained": "BLOB",
		# Declared TEXT since the first version of the schema (values are zlib-compressed bytes;
		# SQLite's flexible typing stores them fine): kept as is to match already-existing databases.
		**{col: "TEXT" for col in IADHORE_COLUMN.values()},
	},
)


def species_exists(db_path, latin_name):
	"""
	Check whether a species already has a row in the paralog Ks table, regardless of which
	analysis types it holds data for.

	:param db_path: path to the sqld server address file (see connection)
	:param latin_name: latin name of species of interest
	:return: True if the species has a row, False otherwise
	:raises KsDatabaseError: if the database can't be read
	"""
	return row_exists(db_path, PARALOG_TABLE, latin_name)


def existing_data_types(db_path, latin_name):
	"""
	Check which analysis types (and the i-ADHoRe file bundle) a species already has non-null data
	for, without loading the actual blobs/text. Used by the populate command to fill in only what's
	missing for a species already in the database (e.g. reciprocal retention finishes days after
	paranome/anchors were first stored), without touching data that's already there, and by the
	skip-recomputation checks in main.nf.

	:param db_path: path to the sqld server address file (see connection)
	:param latin_name: latin name of species of interest
	:return: dict with keys 'paranome', 'anchors', 'recret', 'iadhore', each True if already
	         populated; 'iadhore' is True only if all 5 files are present, matching how they're
	         always read/written as one bundle elsewhere. All False if the species has no row yet.
	:raises KsDatabaseError: if the database can't be read
	"""
	columns = list(_ANALYSIS_TYPE_TO_COLUMN.values()) + [IADHORE_COLUMN[key] for key in IADHORE_FILES]
	present = existing_columns(db_path, PARALOG_TABLE, latin_name, columns)
	return {
		"paranome": present["paranome"],
		"anchors": present["anchors"],
		"recret": present["reciprocally_retained"],
		"iadhore": all(present[IADHORE_COLUMN[key]] for key in IADHORE_FILES),
	}


def read_analysis_data(db_path, latin_name, analysis_type):
	"""
	Look up one species' stored Ks data for one analysis type, without loading the rest of the
	database. This is the single shared entry point all consumers should use to read from the
	paralog Ks table.

	:param db_path: path to the sqld server address file (see connection)
	:param latin_name: latin name of species of interest
	:param analysis_type: one of 'paranome', 'anchors', 'recret'
	:return: dict of lists (Family, Node, Ks, AlignmentCoverage, AlignmentIdentity, AlignmentLength,
	         Paralog1, Paralog2), or None if the species isn't in the database or has no data for
	         this analysis type
	:raises KsDatabaseError: if the database can't be read
	"""
	row = read_columns(db_path, PARALOG_TABLE, latin_name, [_ANALYSIS_TYPE_TO_COLUMN[analysis_type]])
	if row is None or row[0] is None:
		return None
	return _decompress(row[0])


def read_anchor_iadhore_files(db_path, latin_name):
	"""
	Look up one species' stored i-ADHoRe output file contents (anchorpoints.txt, multiplicons.txt,
	segments.txt, list_elements.txt, multiplicon_pairs.txt), without loading the rest of the database.
	These are stored as raw text rather than parsed, so they can be fed into the existing line-based
	parsers in fc_cluster_anchors.py via io.StringIO exactly as a real file would be.

	:param db_path: path to the sqld server address file (see connection)
	:param latin_name: latin name of species of interest
	:return: dict mapping each of 'anchorpoints', 'multiplicons', 'segments', 'list_elements',
	         'multiplicon_pairs' to its raw file text, or None (per key) if not present
	:raises KsDatabaseError: if the database can't be read
	"""
	row = read_columns(db_path, PARALOG_TABLE, latin_name, [IADHORE_COLUMN[key] for key in IADHORE_FILES])
	if row is None:
		return {key: None for key in IADHORE_FILES}
	return {key: (zlib.decompress(blob).decode('utf-8') if blob is not None else None)
	        for key, blob in zip(IADHORE_FILES, row)}


def write_anchor_iadhore_files(latin_name, iadhore_dict, db_path):
	"""
	Write the i-ADHoRe output file contents (raw text) for one species to the database. Like
	write_to_paralog_db, only the keys present (non-None) in iadhore_dict are overwritten; keys
	not present are left untouched (so a partial dict never wipes out previously stored files).
	Text is compressed with zlib before storage to avoid exceeding parameter size limits.

	:param latin_name: latin name of species of interest
	:param iadhore_dict: dict with any of the keys 'anchorpoints', 'multiplicons', 'segments',
	                     'list_elements', 'multiplicon_pairs', each mapping to raw file text or None
	:param db_path: path to the sqld server address file (see connection)
	:return: True if write succeeded, False if there was nothing to write
	:raises KsDatabaseError: if the write fails
	"""
	values = {IADHORE_COLUMN[key]: (zlib.compress(iadhore_dict[key].encode('utf-8')) if iadhore_dict.get(key) is not None else None)
	          for key in IADHORE_FILES}
	if all(v is None for v in values.values()):
		logging.warning("  No i-ADHoRe output files could be read.")
		return False
	return upsert(db_path, PARALOG_TABLE, latin_name, values)


def extract_anchor_iadhore_files(species_name):
	"""
	Read the raw text of the i-ADHoRe output files (anchorpoints.txt, multiplicons.txt, segments.txt,
	list_elements.txt, multiplicon_pairs.txt) generated by the wgd colinearity pipeline for one
	species, from their standard output location. Used to store them in the database right after
	they are generated (wgd_paralogs.py), so cluster_anchor_ks.py doesn't need to re-read them from
	disk later.

	:param species_name: species of interest (informal name)
	:return: dict mapping each of 'anchorpoints', 'multiplicons', 'segments', 'list_elements',
	         'multiplicon_pairs' to its raw file text, or None (per key) if the file is missing
	"""
	wgd_i_adhore_dir = os.path.join("paralog_distributions", f"wgd_{species_name}", f"{species_name}_i-adhore")
	iadhore_texts = {}
	for key in IADHORE_FILES:
		filename = f"{key}.txt"
		path = os.path.join(wgd_i_adhore_dir, filename)
		if os.path.isfile(path):
			with open(path, "r") as f:
				iadhore_texts[key] = f.read()
		else:
			logging.warning(f"  - {filename} not found [{path}].")
			iadhore_texts[key] = None
	return iadhore_texts


def _extract_columns_from_tsv(tsv_path):
	"""
	Read a wgd Ks TSV file and extract the columns needed to recompute weights later
	(Family, Node, Ks, AlignmentCoverage, AlignmentIdentity, AlignmentLength, Paralog1, Paralog2), as a dict of lists.
	All rows are kept unfiltered: filtering by Ks range or alignment quality happens at analysis time.

	:param tsv_path: path to a wgd output Ks TSV file
	:return: dict mapping each required column name to its list of values
	"""
	with open(tsv_path, "r") as f:
		tsv = read_csv(f, sep="\t")
	return {col: tsv[col].to_list() for col in _REQUIRED_COLUMNS}


def extract_paralog_ks_from_tsv(species_name, paranome_enabled=False, anchors_enabled=False,
                                 reciprocal_retention_enabled=False, num_gfs=None, rank_type="lambda", bottom=False):
	"""
	Extract paralog Ks data from TSV files (paranome, anchors, reciprocally_retained), keeping the
	columns needed to recompute weights later for an arbitrary Ks/alignment cutoff. All rows are
	extracted without filtering (filtering is applied at analysis/plotting time).

	:param species_name: species of interest (informal name)
	:param paranome_enabled: whether to extract paranome Ks values
	:param anchors_enabled: whether to extract anchor pair Ks values
	:param reciprocal_retention_enabled: whether to extract reciprocally retained Ks values
	:param num_gfs: number of top/bottom gene families (used only if reciprocal_retention_enabled is True)
	:param rank_type: ranking type for reciprocally retained gene families (default: "lambda")
	:param bottom: whether to use bottom gene families instead of top (default: False)
	:return: dictionary with keys 'paranome', 'anchors', 'recret', each either None (not extracted or
	         file not found) or a dict of lists with keys Family, Node, Ks, AlignmentCoverage,
	         AlignmentIdentity, AlignmentLength, Paralog1, Paralog2
	"""
	# Imported here (not at module top) so that the ortholog commands, which also import this
	# module, don't pull in the whole wgd machinery just for three filename patterns.
	from ksrates.fc_wgd import _OUTPUT_KS_FILE_PATTERN_PARA, _OUTPUT_KS_FILE_PATTERN_ANCHORS, _OUTPUT_KS_FILE_PATTERN_RR_OMCL

	ks_data = {'paranome': None, 'anchors': None, 'recret': None}
	paralog_dists_dir = os.path.join("paralog_distributions", "")

	# Extract paranome Ks data
	if paranome_enabled:
		default_path_paranome = os.path.join(paralog_dists_dir, f"wgd_{species_name}", _OUTPUT_KS_FILE_PATTERN_PARA.format(species_name))
		paranome_path = fcCheck.check_file_existence_and_content_in_default_paths(default_path_paranome, "Paralog (whole-paranome) Ks TSV file")
		if paranome_path != "":
			logging.info("  - Extracting paranome Ks data from TSV file")
			ks_data['paranome'] = _extract_columns_from_tsv(paranome_path)
		else:
			logging.warning(f"  - Paranome Ks TSV file not found [{_OUTPUT_KS_FILE_PATTERN_PARA.format(species_name)}].")

	# Extract anchor pairs Ks data
	if anchors_enabled:
		default_path_anchors = os.path.join(paralog_dists_dir, f"wgd_{species_name}", _OUTPUT_KS_FILE_PATTERN_ANCHORS.format(species_name))
		anchors_path = fcCheck.check_file_existence_and_content_in_default_paths(default_path_anchors, "Anchor pair Ks TSV file")
		if anchors_path != "":
			logging.info("  - Extracting anchor Ks data from TSV file")
			ks_data['anchors'] = _extract_columns_from_tsv(anchors_path)
		else:
			logging.warning(f"  - Anchor Ks TSV file not found [{_OUTPUT_KS_FILE_PATTERN_ANCHORS.format(species_name)}].")

	# Extract reciprocally retained Ks data
	if reciprocal_retention_enabled:
		wgd_species_dir = os.path.join(paralog_dists_dir, f"wgd_{species_name}")
		top_or_bottom = "bottom" if bottom else "top"
		default_path_recret = os.path.join(wgd_species_dir, _OUTPUT_KS_FILE_PATTERN_RR_OMCL.format(species_name, top_or_bottom, num_gfs))
		recret_path = fcCheck.check_file_existence_and_content_in_default_paths(default_path_recret, f"Reciprocally retained paralog Ks TSV file (top {num_gfs})")

		if recret_path != "":
			logging.info("  - Extracting reciprocally retained Ks data from TSV file")
			ks_data['recret'] = _extract_columns_from_tsv(recret_path)
		else:
			logging.warning(f"  - Reciprocally retained Ks TSV file not found [{_OUTPUT_KS_FILE_PATTERN_RR_OMCL.format(species_name, top_or_bottom, num_gfs)}].")

	return ks_data


def write_to_paralog_db(latin_name, ks_data_dict, db_path):
	"""
	Write consolidated paralog Ks data to the database. Weights are not stored: they are recomputed
	at analysis time from the stored Family/Node/alignment columns for whichever Ks/alignment cutoff
	is requested.

	If the species already has a row, only the analysis types present in ks_data_dict (non-None) are
	overwritten; analysis types not present here are left untouched. This matters because this function
	may be called multiple times for the same species with different analyses enabled (e.g. a first run
	with only paranome, a later run adding colinearity) - without this, a later call would otherwise wipe
	out data for analysis types it wasn't asked to extract.

	:param latin_name: latin name of species of interest
	:param ks_data_dict: dictionary from extract_paralog_ks_from_tsv with keys 'paranome', 'anchors', 'recret'
	                     each value is either None or a dict of lists (Family, Node, Ks, AlignmentCoverage,
	                     AlignmentIdentity, AlignmentLength, Paralog1, Paralog2)
	:param db_path: path to the sqld server address file (see connection)
	:return: True if write succeeded, False if there was nothing to write
	:raises KsDatabaseError: if the write fails
	"""
	values = {_ANALYSIS_TYPE_TO_COLUMN[analysis_type]: (_compress(ks_data_dict[analysis_type]) if ks_data_dict[analysis_type] is not None else None)
	          for analysis_type in _ANALYSIS_TYPE_TO_COLUMN}
	if all(v is None for v in values.values()):
		logging.warning(f"  No paralog Ks data could be extracted.")
		return False
	logging.info("  - Writing consolidated Ks lists to database")
	return upsert(db_path, PARALOG_TABLE, latin_name, values)


def consolidate_paralog_ks_lists(species_name, latin_name, db_path, paranome_enabled=False,
                                  anchors_enabled=False, reciprocal_retention_enabled=False,
                                  num_gfs=None, rank_type="lambda", bottom=False):
	"""
	Consolidates paralog Ks lists from three separate TSV files (paranome, anchors, reciprocally_retained)
	into a single database row. All Ks values are stored without filtering to allow flexible filtering
	at the analysis stage.

	:param species_name: species of interest (informal name)
	:param latin_name: latin name of species of interest
	:param db_path: path to the sqld server address file (see connection)
	:param paranome_enabled: whether paranome Ks values are present
	:param anchors_enabled: whether anchor pair Ks values are present
	:param reciprocal_retention_enabled: whether reciprocally retained Ks values are present
	:param num_gfs: number of top/bottom gene families (used only if reciprocal_retention_enabled is True)
	:param rank_type: ranking type for reciprocally retained gene families (default: "lambda")
	:param bottom: whether to use bottom gene families instead of top (default: False)
	:return: True if consolidation failed (nothing could be written), False otherwise
	:raises KsDatabaseError: if the write fails
	"""
	logging.info(f"{species_name} [{latin_name}]:")
	logging.info("- Consolidating paralog Ks lists")

	ks_data = extract_paralog_ks_from_tsv(species_name, paranome_enabled, anchors_enabled,
	                                      reciprocal_retention_enabled, num_gfs, rank_type, bottom)

	write_succeeded = write_to_paralog_db(latin_name, ks_data, db_path)
	logging.info("Done")
	return not write_succeeded


# ======================================================================================
# (C) ORTHOLOG Ks DATA (table ortholog_ks, one row per species pair, key = pair_key)
# ======================================================================================

ORTHOLOG_TABLE = TableSpec(
	name="ortholog_ks",
	key_column="pair_key",
	columns={
		"latin_sp1": "TEXT",   # alphabetically-first (case-insensitive) latin name of the pair
		"latin_sp2": "TEXT",   # the other one
		"mode": "REAL",        # bootstrap-estimated Ks distribution mode
		"mode_sd": "REAL",     # standard deviation of the mode across bootstrap replicates
		"ks_values": "BLOB",   # zlib-compressed pickled list of floats
	},
)

# Column names of the legacy ortholog TSV databases (peak_database_path / ks_list_database_path)
PEAK_TSV_COLUMNS = ['Species1', 'Species2', 'Mode', 'Mode_SD']
KS_LIST_TSV_COLUMNS = ['Species1', 'Species2', 'Ks_Values']


def pair_key(latin_sp1, latin_sp2):
	"""
	Primary key of a species pair: the two latin names sorted alphabetically (case-insensitive) and
	joined by an underscore. Identical to the index of the legacy ortholog TSV databases.

	:return: pair key string
	"""
	sorted_names = sorted([latin_sp1, latin_sp2], key=str.casefold)
	return f"{sorted_names[0]}_{sorted_names[1]}"


def _sorted_latin_names(latin_sp1, latin_sp2):
	return sorted([latin_sp1, latin_sp2], key=str.casefold)


def existing_ortholog_data(db_path, latin_sp1, latin_sp2):
	"""
	Check whether a species pair already has a peak and/or a Ks list, without loading the Ks list.

	:param db_path: path to the sqld server address file (see connection)
	:return: dict with keys 'peak' and 'ks_list', each True if already populated (both False if the
	         pair has no row yet)
	:raises KsDatabaseError: if the database can't be read
	"""
	present = existing_columns(db_path, ORTHOLOG_TABLE, pair_key(latin_sp1, latin_sp2), ["mode", "mode_sd", "ks_values"])
	return {"peak": present["mode"] and present["mode_sd"], "ks_list": present["ks_values"]}


def read_peak(db_path, latin_sp1, latin_sp2):
	"""
	:return: (mode, mode_sd) of the pair, or None if the pair has no stored peak
	:raises KsDatabaseError: if the database can't be read
	"""
	row = read_columns(db_path, ORTHOLOG_TABLE, pair_key(latin_sp1, latin_sp2), ["mode", "mode_sd"])
	if row is None or row[0] is None or row[1] is None:
		return None
	return row[0], row[1]


def read_ks_list(db_path, latin_sp1, latin_sp2):
	"""
	:return: list of floats (the pair's Ks values), or None if the pair has no stored Ks list
	:raises KsDatabaseError: if the database can't be read
	"""
	row = read_columns(db_path, ORTHOLOG_TABLE, pair_key(latin_sp1, latin_sp2), ["ks_values"])
	if row is None or row[0] is None:
		return None
	return _decompress(row[0])


def pair_keys_in_database(db_path, what):
	"""
	List the keys of the species pairs having a given kind of data, without loading the data
	themselves. This is all that is needed by the consumers that only ask "is this pair known yet?".

	:param db_path: path to the sqld server address file (see connection)
	:param what: 'peak' or 'ks_list'
	:return: set of pair keys
	:raises KsDatabaseError: if the database can't be read
	"""
	column = {"peak": "mode", "ks_list": "ks_values"}[what]
	return set(list_keys(db_path, ORTHOLOG_TABLE, not_null_column=column))


def read_all_peaks(db_path):
	"""
	Read all the stored ortholog peaks, as a DataFrame with the same shape as the legacy peak TSV
	database read via pandas.read_csv(path, sep="\\t", index_col=0): indexed by pair key, with columns
	Species1, Species2, Mode, Mode_SD. An empty DataFrame (with those columns) means "no peaks stored yet".

	:param db_path: path to the sqld server address file (see connection)
	:return: DataFrame
	:raises KsDatabaseError: if the database can't be read
	"""
	with connection(db_path) as client:
		result = client.execute("SELECT pair_key, latin_sp1, latin_sp2, mode, mode_sd FROM ortholog_ks "
		                        "WHERE mode IS NOT NULL AND mode_sd IS NOT NULL ORDER BY pair_key")
	df = DataFrame([list(row) for row in result.rows], columns=["pair_key"] + PEAK_TSV_COLUMNS).set_index("pair_key")
	df.index.name = None # same as the unnamed index column of the legacy TSV
	return df


def read_ks_lists(db_path, keys):
	"""
	Read the stored Ks lists of the given species pairs, as a DataFrame with the shape of the legacy
	Ks list TSV database (indexed by pair key, columns Species1, Species2, Ks_Values) except that the
	Ks_Values cells are already real Python lists, NOT strings to be passed to ast.literal_eval.
	Pairs without a stored Ks list are simply absent from the result. Rows are fetched one at a time
	because Ks lists can be large (several can exceed the websocket message size limit together).

	:param db_path: path to the sqld server address file (see connection)
	:param keys: iterable of pair keys to read
	:return: DataFrame
	:raises KsDatabaseError: if the database can't be read
	"""
	rows = []
	with connection(db_path) as client:
		for key in keys:
			result = client.execute("SELECT pair_key, latin_sp1, latin_sp2, ks_values FROM ortholog_ks "
			                        "WHERE pair_key = ? AND ks_values IS NOT NULL", (key,))
			for row in result.rows:
				rows.append([row[0], row[1], row[2], _decompress(row[3])])
	df = DataFrame(rows, columns=["pair_key"] + KS_LIST_TSV_COLUMNS).set_index("pair_key")
	df.index.name = None # same as the unnamed index column of the legacy TSV
	return df


def write_to_ortholog_db(db_path, latin_sp1, latin_sp2, peak=None, ks_list=None):
	"""
	Write the peak and/or the Ks list of a species pair in one atomic statement. Only what is given
	is written: whatever the row already holds (e.g. written by another concurrent run) is left untouched.

	:param db_path: path to the sqld server address file (see connection)
	:param latin_sp1, latin_sp2: latin names of the two species (in any order)
	:param peak: (mode, mode_sd) tuple, or None
	:param ks_list: list of Ks values, or None
	:return: True if something was written, False if there was nothing to write
	:raises KsDatabaseError: if the write fails
	"""
	if peak is None and ks_list is None:
		return False
	name1, name2 = _sorted_latin_names(latin_sp1, latin_sp2)
	values = {"latin_sp1": name1, "latin_sp2": name2}
	if peak is not None:
		values["mode"], values["mode_sd"] = float(peak[0]), float(peak[1])
	if ks_list is not None:
		values["ks_values"] = _compress([float(ks) for ks in ks_list])
	return upsert(db_path, ORTHOLOG_TABLE, pair_key(name1, name2), values)


def extract_from_legacy_tsvs(peak_tsv_path, ks_list_tsv_path):
	"""
	Read the legacy ortholog TSV databases, for the populate command only. Already-computed values
	are copied forward as they are: they are never recomputed, since the bootstrap peak estimate is
	random and recomputing it would silently change the stored values every time. Pairs present in
	only one of the two files are returned with the other kind of data set to None.

	:param peak_tsv_path: path to the legacy peak TSV (or None / missing file: no peaks)
	:param ks_list_tsv_path: path to the legacy Ks list TSV (or None / missing file: no Ks lists)
	:return: dict pair_key -> {'latin_sp1', 'latin_sp2', 'peak': (mode, mode_sd) or None, 'ks_list': list or None}
	"""
	from ast import literal_eval
	pairs = {}

	def _entry(key, sp1, sp2):
		return pairs.setdefault(key, {"latin_sp1": sp1, "latin_sp2": sp2, "peak": None, "ks_list": None})

	if peak_tsv_path and os.path.isfile(peak_tsv_path):
		peaks = read_csv(peak_tsv_path, sep="\t", index_col=0)
		for key, row in peaks.iterrows():
			_entry(key, row["Species1"], row["Species2"])["peak"] = (float(row["Mode"]), float(row["Mode_SD"]))
	if ks_list_tsv_path and os.path.isfile(ks_list_tsv_path):
		ks_lists = read_csv(ks_list_tsv_path, sep="\t", index_col=0)
		for key, row in ks_lists.iterrows():
			try:
				ks_list = literal_eval(row["Ks_Values"])
			except Exception:
				logging.warning(f"  Skipping Ks list of [{key}]: unreadable (truncated?) list in [{ks_list_tsv_path}]")
				continue
			_entry(key, row["Species1"], row["Species2"])["ks_list"] = ks_list
	return pairs
