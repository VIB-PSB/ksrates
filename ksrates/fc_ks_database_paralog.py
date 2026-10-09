"""
Paralog part of the Ks database: table "paralog_ks", one row per species (key = latin name),
holding the paranome, anchor-pair and reciprocally retained Ks data, and the i-ADHoRe output files.
"""
import os
import zlib
import logging
from pandas import read_csv
import ksrates.fc_check_input as fcCheck
from ksrates.fc_ks_database_utilities import (TableSpec, row_exists, existing_columns, read_columns, upsert,
                                              compress, decompress)


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
	return decompress(row[0])


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
	values = {_ANALYSIS_TYPE_TO_COLUMN[analysis_type]: (compress(ks_data_dict[analysis_type]) if ks_data_dict[analysis_type] is not None else None)
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

