import csv
import os
import pickle
import zlib
from datetime import datetime
import ksrates.fc_consolidate_paralog_ks as fc_consolidate_paralog_ks

_TABLE = "paralog_ks"
_PAIR_COLUMNS = ['Paralog1', 'Paralog2', 'Family', 'Node', 'Ks', 'AlignmentCoverage', 'AlignmentIdentity', 'AlignmentLength']

# i-ADHoRe output files stored as raw text columns (see fc_consolidate_paralog_ks._IADHORE_FILES):
# these have their own file-specific layout, unrelated to _PAIR_COLUMNS, so they are dumped back out
# as individual text files rather than folded into the flat TSV.
_IADHORE_FILES = ['anchorpoints', 'multiplicons', 'segments', 'list_elements', 'multiplicon_pairs']
_IADHORE_COLUMNS = [f"{name}_txt" for name in _IADHORE_FILES]


def export_full_tsv(db_path, species_filter=None):
	"""
	Dump the full content of the paralog Ks database to disk, next to the address file itself. Two
	kinds of output are produced, both named from the address file's basename with the current
	timestamp appended (e.g. "paralog_ks_server_address.txt" -> "paralog_ks_server_address_YYYYMMDD_HHMMSS.tsv"):
	- a single flat TSV file, one row per gene pair (not per species): columns are latin_name,
	  analysis_type, Paralog1, Paralog2, Family, Node, Ks, AlignmentCoverage, AlignmentIdentity,
	  AlignmentLength.
	- the raw i-ADHoRe output files (anchorpoints.txt, multiplicons.txt, segments.txt,
	  list_elements.txt, multiplicon_pairs.txt) stored per species, written back out one subdirectory
	  per species (named after its latin name) under a matching "..._iadhore_files" directory.
	Useful for inspecting or analyzing the underlying data outside of ksrates (e.g. in Excel, pandas,
	awk), since the database's blobs and text columns themselves aren't directly readable.

	:param db_path: path to the sqld server address file (see fc_consolidate_paralog_ks._connect)
	:param species_filter: if given, only export species whose latin name contains this substring
	                        (case-insensitive)
	"""
	# Ensures the table exists (e.g. a fresh server that no paralogs-ks run has written to yet)
	# before querying it, rather than letting the SELECT below fail with "no such table".
	if not fc_consolidate_paralog_ks.initialize_paralog_db(db_path):
		print(f"Could not use paralog Ks database [{db_path}].")
		return

	db_base = os.path.splitext(os.path.basename(db_path))[0]
	timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
	output_tsv_path = os.path.join(os.path.dirname(db_path), f"{db_base}_{timestamp}.tsv")
	iadhore_dir = os.path.join(os.path.dirname(db_path), f"{db_base}_{timestamp}_iadhore_files")

	# Fetch species names first to avoid fetching all large BLOBs at once (exceeds websocket message limit).
	client = fc_consolidate_paralog_ks._connect(db_path)
	result = client.execute(f"SELECT latin_name FROM {_TABLE} ORDER BY latin_name")
	species_list = [row[0] for row in result.rows]
	client.close()

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
				data = fc_consolidate_paralog_ks.read_analysis_data(db_path, latin_name, analysis_key)
				if data is None:
					continue
				n = len(data.get('Ks', []))
				for i in range(n):
					writer.writerow([latin_name, display_name] + [data[col][i] for col in _PAIR_COLUMNS])
					n_pairs += 1

			# Fetch i-ADHoRe files separately
			iadhore_texts = fc_consolidate_paralog_ks.read_anchor_iadhore_files(db_path, latin_name)
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


def delete_species(db_path, species_filter):
	"""
	Delete species matching species_filter (case-insensitive substring of their latin name) from
	the paralog Ks database, after listing the matches and asking for confirmation. Used to discard
	a species' stored data, e.g. after it was written incompletely or with now-outdated data.

	:param db_path: path to the sqld server address file (see fc_consolidate_paralog_ks._connect)
	:param species_filter: substring to match against latin names (case-insensitive); required, to
	                        avoid accidentally deleting the entire database's content
	"""
	if not fc_consolidate_paralog_ks.initialize_paralog_db(db_path):
		print(f"Could not use paralog Ks database [{db_path}].")
		return

	client = fc_consolidate_paralog_ks._connect(db_path)
	result = client.execute(f"SELECT latin_name FROM {_TABLE} ORDER BY latin_name")
	all_species = [row[0] for row in result.rows]

	matches = [name for name in all_species if species_filter.lower() in name.lower()]
	if not matches:
		client.close()
		print(f"No species matching '{species_filter}' found in database.")
		return

	print(f"The following {len(matches)} species will be deleted from the database:")
	for name in matches:
		print(f"  {name}")
	text = input("Confirm deleting (y/N)? ").lower()
	if text not in ("y", "yes"):
		client.close()
		print("Cancelled")
		return

	placeholders = ', '.join(['?'] * len(matches))
	client.execute(f"DELETE FROM {_TABLE} WHERE latin_name IN ({placeholders})", tuple(matches))
	client.close()
	print(f"Deleted {len(matches)} species from database.")


def list_species(db_path, species_filter=None):
	"""
	Export a TSV listing which species are in the paralog Ks database and which analysis
	types each has data for (paranome, anchors, reciprocally retained, i-ADHoRe files). Unlike
	export_full_tsv, this only checks presence/absence and is meant for a quick overview.

	:param db_path: path to the sqld server address file (see fc_consolidate_paralog_ks._connect)
	:param species_filter: if given, only list species whose latin name contains this substring
	                        (case-insensitive)
	"""
	if not fc_consolidate_paralog_ks.initialize_paralog_db(db_path):
		print(f"Could not use paralog Ks database [{db_path}].")
		return

	db_base = os.path.splitext(os.path.basename(db_path))[0]
	timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
	output_tsv_path = os.path.join(os.path.dirname(db_path), f"{db_base}_species_list_{timestamp}.tsv")

	client = fc_consolidate_paralog_ks._connect(db_path)
	iadhore_presence = " OR ".join(f"{col} IS NOT NULL" for col in _IADHORE_COLUMNS)
	result = client.execute(f"""
		SELECT latin_name, paranome IS NOT NULL, anchors IS NOT NULL, reciprocally_retained IS NOT NULL,
		       ({iadhore_presence})
		FROM {_TABLE} ORDER BY latin_name
	""")
	client.close()

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
