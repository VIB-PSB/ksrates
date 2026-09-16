import os
import glob
import pickle
import sqlite3
import logging
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


def _connect(db_path):
	"""
	Open a connection to the paralog Ks SQLite database. Uses WAL mode and a generous busy timeout
	so that multiple processes (e.g. one per species in a batch/cluster run) writing to the same
	shared database file don't fail outright on lock contention.

	:param db_path: path to the paralog Ks database file
	:return: sqlite3 connection
	"""
	conn = sqlite3.connect(db_path, timeout=30)
	conn.execute("PRAGMA journal_mode=WAL")
	return conn


def initialize_paralog_db(db_path):
	"""
	Initialize the paralog Ks SQLite database if it doesn't exist yet (creates the file and table).

	:param db_path: path to the paralog Ks database file
	:return: True if the database is ready to use, False if it could not be created/opened
	         (e.g. parent directory missing, permission denied). Callers are responsible for acting on failure.
	"""
	try:
		conn = _connect(db_path)
		conn.execute(f"""
			CREATE TABLE IF NOT EXISTS {_TABLE} (
				latin_name TEXT PRIMARY KEY,
				paranome BLOB,
				anchors BLOB,
				reciprocally_retained BLOB
			)
		""")
		conn.commit()
		conn.close()
		return True
	except Exception as e:
		logging.error(f"Could not create/open paralog Ks database at [{db_path}]: {str(e)}")
		return False


def species_exists(db_path, latin_name):
	"""
	Check whether a species already has a row in the paralog Ks database, regardless of which
	analysis types it holds data for.

	:param db_path: path to the paralog Ks database file
	:param latin_name: latin name of species of interest
	:return: True if the species has a row, False otherwise (including on any read error)
	"""
	try:
		conn = _connect(db_path)
		row = conn.execute(f"SELECT 1 FROM {_TABLE} WHERE latin_name = ? LIMIT 1", (latin_name,)).fetchone()
		conn.close()
		return row is not None
	except Exception:
		return False


def read_analysis_data(db_path, latin_name, analysis_type):
	"""
	Look up one species' stored Ks data for one analysis type, without loading the rest of the
	database. This is the single shared entry point all consumers should use to read from the
	paralog Ks database.

	:param db_path: path to the paralog Ks database file
	:param latin_name: latin name of species of interest
	:param analysis_type: one of 'paranome', 'anchors', 'recret'
	:return: dict of lists (Family, Node, Ks, AlignmentCoverage, AlignmentIdentity, AlignmentLength),
	         or None if the species isn't in the database, has no data for this analysis type, or the
	         database couldn't be read
	"""
	column = _ANALYSIS_TYPE_TO_COLUMN[analysis_type]
	try:
		conn = _connect(db_path)
		row = conn.execute(f"SELECT {column} FROM {_TABLE} WHERE latin_name = ?", (latin_name,)).fetchone()
		conn.close()
	except Exception as e:
		logging.warning(f"Could not read from paralog Ks database [{db_path}]: {str(e)}")
		return None
	if row is None or row[0] is None:
		return None
	return pickle.loads(row[0])


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
	:param db_path: path to the paralog Ks database file
	:return: True if write succeeded, False otherwise
	"""
	ks_paranome = ks_data_dict['paranome']
	ks_anchors = ks_data_dict['anchors']
	ks_recret = ks_data_dict['recret']

	if ks_paranome is None and ks_anchors is None and ks_recret is None:
		logging.warning(f"  No paralog Ks data could be extracted.")
		return False

	logging.info("  - Writing consolidated Ks lists to database")
	blob_paranome = pickle.dumps(ks_paranome) if ks_paranome is not None else None
	blob_anchors = pickle.dumps(ks_anchors) if ks_anchors is not None else None
	blob_recret = pickle.dumps(ks_recret) if ks_recret is not None else None

	conn = _connect(db_path)
	conn.execute(f"""
		INSERT INTO {_TABLE} (latin_name, paranome, anchors, reciprocally_retained)
		VALUES (?, ?, ?, ?)
		ON CONFLICT(latin_name) DO UPDATE SET
			paranome=COALESCE(excluded.paranome, {_TABLE}.paranome),
			anchors=COALESCE(excluded.anchors, {_TABLE}.anchors),
			reciprocally_retained=COALESCE(excluded.reciprocally_retained, {_TABLE}.reciprocally_retained)
	""", (latin_name, blob_paranome, blob_anchors, blob_recret))
	conn.commit()
	conn.close()
	return True


def consolidate_paralog_ks_lists(species_name, latin_name, ks_list_paralog_db_path, paranome_enabled=False,
                                  anchors_enabled=False, reciprocal_retention_enabled=False,
                                  num_gfs=None, rank_type="lambda", bottom=False):
	"""
	Consolidates paralog Ks lists from three separate TSV files (paranome, anchors, reciprocally_retained)
	into a single database row. All Ks values are stored without filtering to allow flexible filtering
	at the analysis stage.

	:param species_name: species of interest (informal name)
	:param ks_list_paralog_db_path: filename/path to the consolidated paralog Ks list database
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