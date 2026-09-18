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
	db_base = os.path.splitext(os.path.basename(db_path))[0]
	timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
	output_tsv_path = os.path.join(os.path.dirname(db_path), f"{db_base}_{timestamp}.tsv")
	iadhore_dir = os.path.join(os.path.dirname(db_path), f"{db_base}_{timestamp}_iadhore_files")

	client = fc_consolidate_paralog_ks._connect(db_path)
	columns = ['latin_name', 'paranome', 'anchors', 'reciprocally_retained'] + _IADHORE_COLUMNS
	result = client.execute(f"SELECT {', '.join(columns)} FROM {_TABLE} ORDER BY latin_name")
	rows = result.rows
	client.close()

	n_pairs = 0
	n_species = 0
	n_iadhore_files = 0
	with open(output_tsv_path, "w", newline="") as outfile:
		writer = csv.writer(outfile, delimiter="\t")
		writer.writerow(['latin_name', 'analysis_type'] + _PAIR_COLUMNS)

		for latin_name, paranome_blob, anchors_blob, recret_blob, *iadhore_texts in rows:
			if species_filter and species_filter.lower() not in latin_name.lower():
				continue
			n_species += 1
			for analysis_type, blob in [("paranome", paranome_blob), ("anchors", anchors_blob), ("reciprocally_retained", recret_blob)]:
				if blob is None:
					continue
				data = pickle.loads(zlib.decompress(blob))
				n = len(data.get('Ks', []))
				for i in range(n):
					writer.writerow([latin_name, analysis_type] + [data[col][i] for col in _PAIR_COLUMNS])
					n_pairs += 1

			for file_name, text in zip(_IADHORE_FILES, iadhore_texts):
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
