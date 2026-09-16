import csv
import pickle
import sqlite3

_TABLE = "paralog_ks"
_PAIR_COLUMNS = ['Paralog1', 'Paralog2', 'Family', 'Node', 'Ks', 'AlignmentCoverage', 'AlignmentIdentity', 'AlignmentLength']


def export_full_tsv(db_path, output_tsv_path, species_filter=None):
	"""
	Dump the full content of the paralog Ks database into a single flat TSV file, one row per
	gene pair (not per species): columns are latin_name, analysis_type, Paralog1, Paralog2, Family, Node, Ks,
	AlignmentCoverage, AlignmentIdentity, AlignmentLength. Useful for inspecting or analyzing the
	underlying data outside of ksrates (e.g. in Excel, pandas, awk), since the SQLite blobs
	themselves aren't human-readable.

	:param db_path: path to the paralog Ks SQLite database file
	:param output_tsv_path: path where the flat TSV dump will be written
	:param species_filter: if given, only export species whose latin name contains this substring
	                        (case-insensitive)
	"""
	conn = sqlite3.connect(db_path)
	rows = conn.execute(f"SELECT latin_name, paranome, anchors, reciprocally_retained FROM {_TABLE} ORDER BY latin_name").fetchall()
	conn.close()

	n_pairs = 0
	n_species = 0
	with open(output_tsv_path, "w", newline="") as outfile:
		writer = csv.writer(outfile, delimiter="\t")
		writer.writerow(['latin_name', 'analysis_type'] + _PAIR_COLUMNS)

		for latin_name, paranome_blob, anchors_blob, recret_blob in rows:
			if species_filter and species_filter.lower() not in latin_name.lower():
				continue
			n_species += 1
			for analysis_type, blob in [("paranome", paranome_blob), ("anchors", anchors_blob), ("reciprocally_retained", recret_blob)]:
				if blob is None:
					continue
				data = pickle.loads(blob)
				n = len(data.get('Ks', []))
				for i in range(n):
					writer.writerow([latin_name, analysis_type] + [data[col][i] for col in _PAIR_COLUMNS])
					n_pairs += 1

	print(f"Exported {n_pairs} gene pairs from {n_species} species to [{output_tsv_path}]")
