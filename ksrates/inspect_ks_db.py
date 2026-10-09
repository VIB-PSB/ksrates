import csv
import os
import pickle
import zlib
from datetime import datetime
import ksrates.fc_ks_database_utilities as fcKsDbUtil
import ksrates.fc_ks_database_paralog as fcKsDbPara
import ksrates.fc_ks_database_ortholog as fcKsDbOrtho

_TABLE = "paralog_ks"
_PAIR_COLUMNS = ['Paralog1', 'Paralog2', 'Family', 'Node', 'Ks', 'AlignmentCoverage', 'AlignmentIdentity', 'AlignmentLength']

# i-ADHoRe output files stored as raw text columns (see fc_ks_database_paralog.IADHORE_FILES):
# these have their own file-specific layout, unrelated to _PAIR_COLUMNS, so they are dumped back out
# as individual text files rather than folded into the flat TSV.
_IADHORE_FILES = ['anchorpoints', 'multiplicons', 'segments', 'list_elements', 'multiplicon_pairs']
_IADHORE_COLUMNS = [f"{name}_txt" for name in _IADHORE_FILES]


def export_paralog_tsv(db_path, species_filter=None):
	"""
	Dump the full paralog content of the Ks database to disk, next to the address file itself. Two
	kinds of output are produced, both named from the address file's basename with the current
	timestamp appended (e.g. "ks_database_server_address.txt" -> "ks_database_server_address_YYYYMMDD_HHMMSS.tsv"):
	- a single flat TSV file, one row per gene pair (not per species): columns are latin_name,
	  analysis_type, Paralog1, Paralog2, Family, Node, Ks, AlignmentCoverage, AlignmentIdentity,
	  AlignmentLength.
	- the raw i-ADHoRe output files (anchorpoints.txt, multiplicons.txt, segments.txt,
	  list_elements.txt, multiplicon_pairs.txt) stored per species, written back out one subdirectory
	  per species (named after its latin name) under a matching "..._iadhore_files" directory.
	Useful for inspecting or analyzing the underlying data outside of ksrates (e.g. in Excel, pandas,
	awk), since the database's blobs and text columns themselves aren't directly readable.

	:param db_path: path to the sqld server address file (see fcKsDbUtil.connection)
	:param species_filter: if given, only export species whose latin name contains this substring
	                        (case-insensitive)
	"""
	# Ensures the table exists (e.g. a fresh server that no paralogs-ks run has written to yet)
	# before querying it, rather than letting the SELECT below fail with "no such table".
	fcKsDbUtil.initialize_database(db_path)

	db_base = os.path.splitext(os.path.basename(db_path))[0]
	timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
	output_tsv_path = os.path.join(os.path.dirname(db_path), f"{db_base}_{timestamp}.tsv")
	iadhore_dir = os.path.join(os.path.dirname(db_path), f"{db_base}_{timestamp}_iadhore_files")

	# Fetch species names first to avoid fetching all large BLOBs at once (exceeds websocket message limit).
	with fcKsDbUtil.connection(db_path) as client:
		result = client.execute(f"SELECT latin_name FROM {_TABLE} ORDER BY latin_name")
		species_list = [row[0] for row in result.rows]

	n_pairs = 0
	n_species = 0
	n_iadhore_files = 0
	with open(output_tsv_path, "w", newline="") as outfile:
		writer = csv.writer(outfile, delimiter="\t")
		writer.writerow(['latin_name', 'analysis_type'] + _PAIR_COLUMNS)

		# Fetch each species' data one at a time using read functions (may handle large BLOBs better)
		for latin_name in species_list:
			if species_filter and species_filter.lower() not in latin_name.lower():
				continue

			n_species += 1

			# Fetch paranome, anchors, recret separately to avoid oversized messages
			for analysis_key, display_name in [('paranome', 'paranome'), ('anchors', 'anchors'), ('recret', 'reciprocally_retained')]:
				data = fcKsDbPara.read_analysis_data(db_path, latin_name, analysis_key)
				if data is None:
					continue
				n = len(data.get('Ks', []))
				for i in range(n):
					writer.writerow([latin_name, display_name] + [data[col][i] for col in _PAIR_COLUMNS])
					n_pairs += 1

			# Fetch i-ADHoRe files separately
			iadhore_texts = fcKsDbPara.read_anchor_iadhore_files(db_path, latin_name)
			for file_name in _IADHORE_FILES:
				text = iadhore_texts.get(file_name)
				if text is None:
					continue
				species_dir = os.path.join(iadhore_dir, latin_name.replace(" ", "_"))
				os.makedirs(species_dir, exist_ok=True)
				with open(os.path.join(species_dir, f"{file_name}.txt"), "w") as iadhore_file:
					iadhore_file.write(text)
				n_iadhore_files += 1

	print(f"Exported {n_pairs} gene pairs from {n_species} species to [{output_tsv_path}]")
	if n_iadhore_files:
		print(f"Exported {n_iadhore_files} i-ADHoRe output files to [{iadhore_dir}]")


def delete_paralog_species(db_path, species_filter):
	"""
	Delete species matching species_filter (case-insensitive substring of their latin name) from
	the Ks database, after listing the matches and asking for confirmation. Used to discard
	a species' stored data, e.g. after it was written incompletely or with now-outdated data.

	:param db_path: path to the sqld server address file (see fcKsDbUtil.connection)
	:param species_filter: substring to match against latin names (case-insensitive); required, to
	                        avoid accidentally deleting the entire database's content
	"""
	fcKsDbUtil.initialize_database(db_path)

	with fcKsDbUtil.connection(db_path) as client:
		result = client.execute(f"SELECT latin_name FROM {_TABLE} ORDER BY latin_name")
		all_species = [row[0] for row in result.rows]

		matches = [name for name in all_species if species_filter.lower() in name.lower()]
		if not matches:
			print(f"No species matching '{species_filter}' found in database.")
			return

		print(f"The following {len(matches)} species will be deleted from the database:")
		for name in matches:
			print(f"  {name}")
		text = input("Confirm deleting (y/N)? ").lower()
		if text not in ("y", "yes"):
			print("Cancelled")
			return

		placeholders = ', '.join(['?'] * len(matches))
		client.execute(f"DELETE FROM {_TABLE} WHERE latin_name IN ({placeholders})", tuple(matches))
		print(f"Deleted {len(matches)} species from database.")


def list_paralog_species(db_path, species_filter=None):
	"""
	Export a TSV listing which species are in the Ks database and which analysis
	types each has data for (paranome, anchors, reciprocally retained, i-ADHoRe files). Unlike
	export_paralog_tsv, this only checks presence/absence and is meant for a quick overview.

	:param db_path: path to the sqld server address file (see fcKsDbUtil.connection)
	:param species_filter: if given, only list species whose latin name contains this substring
	                        (case-insensitive)
	"""
	fcKsDbUtil.initialize_database(db_path)

	db_base = os.path.splitext(os.path.basename(db_path))[0]
	timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
	output_tsv_path = os.path.join(os.path.dirname(db_path), f"{db_base}_species_list_{timestamp}.tsv")

	with fcKsDbUtil.connection(db_path) as client:
		iadhore_presence = " AND ".join(f"{col} IS NOT NULL" for col in _IADHORE_COLUMNS)
		result = client.execute(f"""
			SELECT latin_name, paranome IS NOT NULL, anchors IS NOT NULL, reciprocally_retained IS NOT NULL,
			       ({iadhore_presence})
			FROM {_TABLE} ORDER BY latin_name
		""")

	n_species = 0
	lines = ['\t'.join(['latin_name', 'paranome', 'anchors', 'reciprocally_retained', 'iadhore_files'])]
	for latin_name, has_paranome, has_anchors, has_recret, has_iadhore in result.rows:
		if species_filter and species_filter.lower() not in latin_name.lower():
			continue
		n_species += 1
		lines.append('\t'.join([latin_name, str(bool(has_paranome)), str(bool(has_anchors)), str(bool(has_recret)), str(bool(has_iadhore))]))

	# No trailing newline after the last row: a naive line-count on this file (e.g. wc -l, or
	# splitting on "\n") should equal exactly 1 header + n_species data lines, not one more.
	with open(output_tsv_path, "w") as outfile:
		outfile.write('\n'.join(lines))

	print(f"Listed {n_species} species to [{output_tsv_path}]")



def _matching_pair_keys(db_path, pair_filter):
	"""
	All the ortholog pair keys of the database whose key contains pair_filter
	(case-insensitive), or all of them if no filter.
	"""
	keys = fcKsDbUtil.list_keys(db_path, fcKsDbOrtho.ORTHOLOG_TABLE)
	if pair_filter:
		keys = [key for key in keys if pair_filter.lower() in key.lower()]
	return keys


def export_ortholog_tsvs(db_path, pair_filter=None):
	"""
	Dump the ortholog content of the Ks database to disk, next to the address file itself, as two TSV files in
	the exact format of the legacy ortholog TSV databases (peak_database_path and ks_list_database_path), so that
	a dump can be used as a drop-in replacement of those files. Both are named from the address file's basename
	with the current timestamp appended (e.g. "ks_database_server_address.txt" ->
	"ks_database_server_address_ortholog_peaks_YYYYMMDD_HHMMSS.tsv" and "..._ortholog_ks_lists_YYYYMMDD_HHMMSS.tsv").

	:param db_path: path to the sqld server address file (see fcKsDbUtil.connection)
	:param pair_filter: if given, only export the species pairs whose key (latin names joined by "_")
	                    contains this substring (case-insensitive)
	"""
	fcKsDbUtil.initialize_database(db_path)

	db_base = os.path.splitext(os.path.basename(db_path))[0]
	timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
	peaks_path = os.path.join(os.path.dirname(db_path), f"{db_base}_ortholog_peaks_{timestamp}.tsv")
	ks_lists_path = os.path.join(os.path.dirname(db_path), f"{db_base}_ortholog_ks_lists_{timestamp}.tsv")

	peaks = fcKsDbOrtho.read_all_peaks(db_path)
	if pair_filter:
		peaks = peaks[[pair_filter.lower() in key.lower() for key in peaks.index]]
	peaks.index.name = None # legacy TSVs have an unnamed index column
	peaks.to_csv(peaks_path, sep="\t")

	# Ks lists are read in chunks of pairs, since a single message holding too many of them can exceed the websocket size limit
	pair_keys_with_ks_list = fcKsDbOrtho.pair_keys_in_database(db_path, "ks_list")
	keys_with_ks_list = [key for key in _matching_pair_keys(db_path, pair_filter) if key in pair_keys_with_ks_list]
	with open(ks_lists_path, "w") as outfile:
		outfile.write("\tSpecies1\tSpecies2\tKs_Values\n")
		for key in keys_with_ks_list:
			row = fcKsDbOrtho.read_ks_lists(db_path, [key]).loc[key]
			outfile.write(f"{key}\t{row['Species1']}\t{row['Species2']}\t{row['Ks_Values']}\n")

	print(f"Exported {len(peaks)} ortholog peaks to [{peaks_path}]")
	print(f"Exported {len(keys_with_ks_list)} ortholog Ks lists to [{ks_lists_path}]")


def list_ortholog_pairs(db_path, pair_filter=None):
	"""
	Export a TSV listing which species pairs are in the ortholog part of the Ks database and whether each has
	a peak and a Ks list. Only checks presence/absence, as a quick overview.

	:param db_path: path to the sqld server address file (see fcKsDbUtil.connection)
	:param pair_filter: if given, only list the pairs whose key contains this substring (case-insensitive)
	"""
	fcKsDbUtil.initialize_database(db_path)

	db_base = os.path.splitext(os.path.basename(db_path))[0]
	timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
	output_tsv_path = os.path.join(os.path.dirname(db_path), f"{db_base}_ortholog_pairs_list_{timestamp}.tsv")

	with fcKsDbUtil.connection(db_path) as client:
		result = client.execute("SELECT pair_key, mode IS NOT NULL AND mode_sd IS NOT NULL, ks_values IS NOT NULL FROM ortholog_ks ORDER BY pair_key")

	lines = ['\t'.join(['pair_key', 'peak', 'ks_list'])]
	n_pairs = 0
	for key, has_peak, has_ks_list in result.rows:
		if pair_filter and pair_filter.lower() not in key.lower():
			continue
		n_pairs += 1
		lines.append('\t'.join([key, str(bool(has_peak)), str(bool(has_ks_list))]))
	with open(output_tsv_path, "w") as outfile:
		outfile.write('\n'.join(lines))

	print(f"Listed {n_pairs} species pairs to [{output_tsv_path}]")


def delete_ortholog_pairs(db_path, pair_filter):
	"""
	Delete the species pairs matching pair_filter (case-insensitive substring of their key) from the ortholog part
	of the Ks database, after listing the matches and asking for confirmation. Used to discard the stored data of a
	pair, e.g. to have its peak estimated again.

	:param db_path: path to the sqld server address file (see fcKsDbUtil.connection)
	:param pair_filter: substring to match against pair keys; required, to avoid deleting the entire content
	"""
	fcKsDbUtil.initialize_database(db_path)

	matches = _matching_pair_keys(db_path, pair_filter)
	if not matches:
		print(f"No species pair matching '{pair_filter}' found in database.")
		return

	print(f"The following {len(matches)} species pairs will be deleted from the database:")
	for key in matches:
		print(f"  {key}")
	text = input("Confirm deleting (y/N)? ").lower()
	if text not in ("y", "yes"):
		print("Cancelled")
		return

	fcKsDbUtil.delete_keys(db_path, fcKsDbOrtho.ORTHOLOG_TABLE, matches)
	print(f"Deleted {len(matches)} species pairs from database.")
