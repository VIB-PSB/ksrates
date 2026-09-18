import os
import pickle
import zlib
import logging
import libsql_client
from pandas import read_csv
import ksrates.fc_check_input as fcCheck
from ksrates.fc_wgd import _OUTPUT_KS_FILE_PATTERN_PARA, _OUTPUT_KS_FILE_PATTERN_ANCHORS, _OUTPUT_KS_FILE_PATTERN_RR_OMCL

# Columns needed to (1) recompute per-family/per-node weights later (via filter_compute_weights)
# for an arbitrary Ks/alignment-quality cutoff, and (2) identify which gene pair each Ks value
# belongs to (Paralog1, Paralog2), without having to re-read the original TSV file.
_REQUIRED_COLUMNS = ['Family', 'Node', 'Ks', 'AlignmentCoverage', 'AlignmentIdentity', 'AlignmentLength', 'Paralog1', 'Paralog2']

_TABLE = "paralog_ks"
# Maps the extract_paralog_ks_from_tsv()/ks_data_dict keys to the SQLite column names
_ANALYSIS_TYPE_TO_COLUMN = {"paranome": "paranome", "anchors": "anchors", "recret": "reciprocally_retained"}

# i-ADHoRe output files needed by anchor clustering (cluster_anchor_ks.py), stored as raw text
# rather than parsed: the existing line-based parsers in fc_cluster_anchors.py can read from an
# io.StringIO wrapping this text exactly as they read from a real file, so no parser rewrite is
# needed to support the database as an alternative source.
_IADHORE_FILES = ['anchorpoints', 'multiplicons', 'segments', 'list_elements', 'multiplicon_pairs']
_IADHORE_COLUMN = {key: f"{key}_txt" for key in _IADHORE_FILES}

# All columns beyond the latin_name primary key, with their SQLite type: used both to create a
# brand new database and to add any missing columns to a database created by an older version of
# this schema (plain CREATE TABLE IF NOT EXISTS would not add columns to an existing table).
_ALL_COLUMNS = {
	"paranome": "BLOB",
	"anchors": "BLOB",
	"reciprocally_retained": "BLOB",
	**{col: "TEXT" for col in _IADHORE_COLUMN.values()},
}


def _read_address_file(db_path):
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


def _connect(db_path):
	"""
	Open a connection to the paralog Ks database, which is served remotely by a long-running
	sqld server (self-hosted libSQL server) rather than being a local file: this lets many
	independent cluster jobs, possibly on different compute nodes, safely read/write a shared
	database without relying on SQLite's WAL mode over a network filesystem (which SQLite's own
	docs say doesn't work reliably across hosts).

	:param db_path: path to a small text file, written once at startup by the sqld server job,
	                containing the server's own address as "host:port" on the first line and,
	                if the server was started with JWT auth enabled, the shared auth token on
	                the second line
	:return: a libsql_client sync client connected to the sqld server
	:raises: any error reading/parsing the address file, or (on the first query issued against
	         the returned client, since the connection itself is opened lazily) any error
	         reaching the server, exactly as sqlite3.connect() used to raise on a bad local path.
	         Every caller already wraps its database calls in a broad except Exception, so this
	         requires no changes on the caller side.
	"""
	address, token = _read_address_file(db_path)
	# ws:// (not http://) since http:// explicitly can't support transactions; note that
	# libsql-client is an archived (frozen, unmaintained) package as of writing, but its
	# execute()/close() surface is small, stable, and verified to work against a real sqld
	# instance for our usage pattern.
	return libsql_client.create_client_sync(f"ws://{address}", auth_token=token)


def initialize_paralog_db(db_path):
	"""
	Initialize the paralog Ks database if it doesn't exist yet (creates the table), and add any
	columns missing from a database created by an older version of this schema.

	:param db_path: path to the sqld server address file (see _connect)
	:return: True if the database is ready to use, False if it could not be created/opened
	         (e.g. server unreachable, malformed address file). Callers are responsible for
	         acting on failure.
	"""
	try:
		client = _connect(db_path)
		client.execute(f"CREATE TABLE IF NOT EXISTS {_TABLE} (latin_name TEXT PRIMARY KEY)")
		# Ordinary SELECT against the pragma_table_info table-valued function, rather than a
		# bare "PRAGMA table_info(...)" statement: some PRAGMA statement forms are known to be
		# rejected over sqld's remote/Hrana protocol, while this SELECT form is not.
		result = client.execute(f"SELECT name FROM pragma_table_info('{_TABLE}')")
		existing_columns = {row[0] for row in result.rows}
		for column, col_type in _ALL_COLUMNS.items():
			if column not in existing_columns:
				try:
					client.execute(f"ALTER TABLE {_TABLE} ADD COLUMN {column} {col_type}")
				except Exception as e:
					# Ignore "duplicate column" errors from concurrent ALTER TABLE calls:
					# another process may have added this column between our pragma_table_info
					# query and this ALTER statement. Any other error is real and should propagate.
					if "duplicate column" not in str(e).lower():
						raise
		# No commit needed/possible here: outside of an explicit transaction, the sqld client
		# has no .commit() method, since each execute() above is already committed by the
		# server as soon as it returns.
		client.close()
		return True
	except Exception as e:
		logging.error(f"Could not create/open paralog Ks database at [{db_path}]: {str(e)}")
		return False


def species_exists(db_path, latin_name):
	"""
	Check whether a species already has a row in the paralog Ks database, regardless of which
	analysis types it holds data for.

	:param db_path: path to the sqld server address file (see _connect)
	:param latin_name: latin name of species of interest
	:return: True if the species has a row, False otherwise (including on any read error)
	"""
	try:
		client = _connect(db_path)
		result = client.execute(f"SELECT 1 FROM {_TABLE} WHERE latin_name = ? LIMIT 1", (latin_name,))
		client.close()
		return len(result.rows) > 0
	except Exception:
		return False


def existing_data_types(db_path, latin_name):
	"""
	Check which analysis types (and the i-ADHoRe file bundle) a species already has non-null data
	for, without loading the actual blobs/text. Used by populate_paralog_ks_db_batch to fill in
	only what's missing for a species already in the database (e.g. reciprocal retention finishes
	days after paranome/anchors were first stored), without touching data that's already there.

	:param db_path: path to the sqld server address file (see _connect)
	:param latin_name: latin name of species of interest
	:return: dict with keys 'paranome', 'anchors', 'recret', 'iadhore', each True if already
	         populated; 'iadhore' is True only if all 5 files are present, matching how they're
	         always read/written as one bundle elsewhere. All False if the species has no row yet
	         or the database couldn't be read.
	"""
	empty = {"paranome": False, "anchors": False, "recret": False, "iadhore": False}
	columns = ["paranome", "anchors", "reciprocally_retained"] + [_IADHORE_COLUMN[key] for key in _IADHORE_FILES]
	try:
		client = _connect(db_path)
		result = client.execute(f"SELECT {', '.join(columns)} FROM {_TABLE} WHERE latin_name = ?", (latin_name,))
		client.close()
	except Exception as e:
		logging.warning(f"Could not read from paralog Ks database [{db_path}]: {str(e)}")
		return empty
	if not result.rows:
		return empty
	paranome, anchors, recret, *iadhore_values = result.rows[0]
	return {
		"paranome": paranome is not None,
		"anchors": anchors is not None,
		"recret": recret is not None,
		"iadhore": all(v is not None for v in iadhore_values),
	}


def read_analysis_data(db_path, latin_name, analysis_type):
	"""
	Look up one species' stored Ks data for one analysis type, without loading the rest of the
	database. This is the single shared entry point all consumers should use to read from the
	paralog Ks database.

	:param db_path: path to the sqld server address file (see _connect)
	:param latin_name: latin name of species of interest
	:param analysis_type: one of 'paranome', 'anchors', 'recret'
	:return: dict of lists (Family, Node, Ks, AlignmentCoverage, AlignmentIdentity, AlignmentLength),
	         or None if the species isn't in the database, has no data for this analysis type, or the
	         database couldn't be read
	"""
	column = _ANALYSIS_TYPE_TO_COLUMN[analysis_type]
	try:
		client = _connect(db_path)
		result = client.execute(f"SELECT {column} FROM {_TABLE} WHERE latin_name = ?", (latin_name,))
		client.close()
	except Exception as e:
		logging.warning(f"Could not read from paralog Ks database [{db_path}]: {str(e)}")
		return None
	row = result.rows[0] if result.rows else None
	if row is None or row[0] is None:
		return None
	return pickle.loads(zlib.decompress(row[0]))


def read_anchor_iadhore_files(db_path, latin_name):
	"""
	Look up one species' stored i-ADHoRe output file contents (anchorpoints.txt, multiplicons.txt,
	segments.txt, list_elements.txt, multiplicon_pairs.txt), without loading the rest of the database.
	These are stored as raw text rather than parsed, so they can be fed into the existing line-based
	parsers in fc_cluster_anchors.py via io.StringIO exactly as a real file would be.

	:param db_path: path to the sqld server address file (see _connect)
	:param latin_name: latin name of species of interest
	:return: dict mapping each of 'anchorpoints', 'multiplicons', 'segments', 'list_elements',
	         'multiplicon_pairs' to its raw file text, or None (per key) if not present or on read error
	"""
	columns = [_IADHORE_COLUMN[key] for key in _IADHORE_FILES]
	try:
		client = _connect(db_path)
		result = client.execute(f"SELECT {', '.join(columns)} FROM {_TABLE} WHERE latin_name = ?", (latin_name,))
		client.close()
	except Exception as e:
		logging.warning(f"Could not read from paralog Ks database [{db_path}]: {str(e)}")
		return {key: None for key in _IADHORE_FILES}
	row = result.rows[0] if result.rows else None
	if row is None:
		return {key: None for key in _IADHORE_FILES}
	# Decompress each i-ADHoRe file's text
	decompressed = {}
	for key, blob in zip(_IADHORE_FILES, row):
		if blob is not None:
			decompressed[key] = zlib.decompress(blob).decode('utf-8')
		else:
			decompressed[key] = None
	return decompressed


def write_anchor_iadhore_files(latin_name, iadhore_dict, db_path):
	"""
	Write the i-ADHoRe output file contents (raw text) for one species to the database. Like
	write_to_paralog_db, only the keys present (non-None) in iadhore_dict are overwritten; keys
	not present are left untouched (so a partial dict never wipes out previously stored files).
	Text is compressed with zlib before storage to avoid exceeding parameter size limits.

	:param latin_name: latin name of species of interest
	:param iadhore_dict: dict with any of the keys 'anchorpoints', 'multiplicons', 'segments',
	                       'list_elements', 'multiplicon_pairs', each mapping to raw file text or None
	:param db_path: path to the sqld server address file (see _connect)
	:return: True if write succeeded, False if there was nothing to write
	"""
	values = {key: iadhore_dict.get(key) for key in _IADHORE_FILES}
	if all(v is None for v in values.values()):
		logging.warning("  No i-ADHoRe output files could be read.")
		return False

	columns = [_IADHORE_COLUMN[key] for key in _IADHORE_FILES]
	placeholders = ', '.join(['?'] * len(columns))
	set_clause = ", ".join(f"{col}=COALESCE(excluded.{col}, {_TABLE}.{col})" for col in columns)

	# Compress each i-ADHoRe file's text before storage
	compressed_values = []
	for key in _IADHORE_FILES:
		val = values[key]
		if val is not None:
			compressed_values.append(zlib.compress(val.encode('utf-8')))
		else:
			compressed_values.append(None)

	client = _connect(db_path)
	client.execute(f"""
		INSERT INTO {_TABLE} (latin_name, {', '.join(columns)})
		VALUES (?, {placeholders})
		ON CONFLICT(latin_name) DO UPDATE SET {set_clause}
	""", (latin_name, *compressed_values))
	client.close()
	return True


def extract_anchor_iadhore_files(species_name):
	"""
	Read the raw text of the i-ADHoRe output files (anchorpoints.txt, multiplicons.txt, segments.txt,
	list_elements.txt, multiplicon_pairs.txt) generated by the wgd colinearity pipeline for one
	species, from their standard output location. Used to opportunistically store them in the
	database right after they are generated (wgd_paralogs.py), so cluster_anchor_ks.py doesn't need
	to re-read them from disk later.

	:param species_name: species of interest (informal name)
	:return: dict mapping each of 'anchorpoints', 'multiplicons', 'segments', 'list_elements',
	         'multiplicon_pairs' to its raw file text, or None (per key) if the file is missing
	"""
	wgd_i_adhore_dir = os.path.join("paralog_distributions", f"wgd_{species_name}", f"{species_name}_i-adhore")
	filenames = {
		'anchorpoints': 'anchorpoints.txt',
		'multiplicons': 'multiplicons.txt',
		'segments': 'segments.txt',
		'list_elements': 'list_elements.txt',
		'multiplicon_pairs': 'multiplicon_pairs.txt',
	}
	iadhore_texts = {}
	for key in _IADHORE_FILES:
		path = os.path.join(wgd_i_adhore_dir, filenames[key])
		if os.path.isfile(path):
			with open(path, "r") as f:
				iadhore_texts[key] = f.read()
		else:
			logging.warning(f"  - {filenames[key]} not found [{path}].")
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
	         AlignmentIdentity, AlignmentLength
	"""
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
	with only paranome, a later run adding colinearity) — without this, a later call would otherwise wipe
	out data for analysis types it wasn't asked to extract.

	:param latin_name: latin name of species of interest
	:param ks_data_dict: dictionary from extract_paralog_ks_from_tsv with keys 'paranome', 'anchors', 'recret'
	                     each value is either None or a dict of lists (Family, Node, Ks, AlignmentCoverage,
	                     AlignmentIdentity, AlignmentLength)
	:param db_path: path to the sqld server address file (see _connect)
	:return: True if write succeeded, False otherwise
	"""
	ks_paranome = ks_data_dict['paranome']
	ks_anchors = ks_data_dict['anchors']
	ks_recret = ks_data_dict['recret']

	if ks_paranome is None and ks_anchors is None and ks_recret is None:
		logging.warning(f"  No paralog Ks data could be extracted.")
		return False

	logging.info("  - Writing consolidated Ks lists to database")
	blob_paranome = zlib.compress(pickle.dumps(ks_paranome)) if ks_paranome is not None else None
	blob_anchors = zlib.compress(pickle.dumps(ks_anchors)) if ks_anchors is not None else None
	blob_recret = zlib.compress(pickle.dumps(ks_recret)) if ks_recret is not None else None

	client = _connect(db_path)
	client.execute(f"""
		INSERT INTO {_TABLE} (latin_name, paranome, anchors, reciprocally_retained)
		VALUES (?, ?, ?, ?)
		ON CONFLICT(latin_name) DO UPDATE SET
			paranome=COALESCE(excluded.paranome, {_TABLE}.paranome),
			anchors=COALESCE(excluded.anchors, {_TABLE}.anchors),
			reciprocally_retained=COALESCE(excluded.reciprocally_retained, {_TABLE}.reciprocally_retained)
	""", (latin_name, blob_paranome, blob_anchors, blob_recret))
	client.close()
	return True


def consolidate_paralog_ks_lists(species_name, latin_name, ks_list_paralog_db_path, paranome_enabled=False,
                                  anchors_enabled=False, reciprocal_retention_enabled=False,
                                  num_gfs=None, rank_type="lambda", bottom=False):
	"""
	Consolidates paralog Ks lists from three separate TSV files (paranome, anchors, reciprocally_retained)
	into a single database row. All Ks values are stored without filtering to allow flexible filtering
	at the analysis stage.

	:param species_name: species of interest (informal name)
	:param ks_list_paralog_db_path: path to the sqld server address file (see _connect)
	:param paranome_enabled: whether paranome Ks values are present
	:param anchors_enabled: whether anchor pair Ks values are present
	:param reciprocal_retention_enabled: whether reciprocally retained Ks values are present
	:param num_gfs: number of top/bottom gene families (used only if reciprocal_retention_enabled is True)
	:param rank_type: ranking type for reciprocally retained gene families (default: "lambda")
	:param bottom: whether to use bottom gene families instead of top (default: False)
	:return: True if consolidation failed, False otherwise
	"""
	logging.info(f"{species_name} [{latin_name}]:")
	logging.info("- Consolidating paralog Ks lists")

	ks_data = extract_paralog_ks_from_tsv(species_name, paranome_enabled, anchors_enabled,
										reciprocal_retention_enabled, num_gfs, rank_type, bottom)

	write_to_paralog_db(latin_name, ks_data, ks_list_paralog_db_path)
	return
