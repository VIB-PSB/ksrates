"""
Ortholog part of the Ks database: table "ortholog_ks", one row per species pair (key = pair_key),
holding the bootstrap-estimated Ks distribution peak and the Ks list it was estimated from.
"""
import os
import logging
from pandas import read_csv, DataFrame
from ksrates.fc_ks_database_utilities import (TableSpec, connection, existing_columns, read_columns,
                                              upsert, list_keys, compress, decompress)


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
	return decompress(row[0])


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
				rows.append([row[0], row[1], row[2], decompress(row[3])])
	df = DataFrame(rows, columns=["pair_key"] + KS_LIST_TSV_COLUMNS).set_index("pair_key")
	df.index.name = None # same as the unnamed index column of the legacy TSV
	return df


def write_to_ortholog_db(db_path, latin_sp1, latin_sp2, peak=None, ks_list=None, reset_peak=False):
	"""
	Write the peak and/or the Ks list of a species pair in one atomic statement. Only what is given
	is written: whatever the row already holds (e.g. written by another concurrent run) is left untouched.

	:param db_path: path to the sqld server address file (see connection)
	:param latin_sp1, latin_sp2: latin names of the two species (in any order)
	:param peak: (mode, mode_sd) tuple, or None
	:param ks_list: list of Ks values, or None
	:param reset_peak: erase the stored peak in the same statement. Used when a new Ks list is stored without a new
	                   peak: the old peak was estimated from a different list, and a stored peak must always come
	                   from the stored Ks list
	:return: True if something was written, False if there was nothing to write
	:raises KsDatabaseError: if the write fails
	"""
	if peak is None and ks_list is None and not reset_peak:
		return False
	name1, name2 = _sorted_latin_names(latin_sp1, latin_sp2)
	values = {"latin_sp1": name1, "latin_sp2": name2}
	if peak is not None:
		values["mode"], values["mode_sd"] = float(peak[0]), float(peak[1])
	if ks_list is not None:
		values["ks_values"] = compress([float(ks) for ks in ks_list])
	return upsert(db_path, ORTHOLOG_TABLE, pair_key(name1, name2), values,
	              clear=("mode", "mode_sd") if (reset_peak and peak is None) else ())


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
