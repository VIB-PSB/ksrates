import os
import glob
import logging
from pandas import DataFrame, read_csv
import ksrates.fc_check_input as fcCheck
from ksrates.fc_wgd import _OUTPUT_KS_FILE_PATTERN_PARA, _OUTPUT_KS_FILE_PATTERN_ANCHORS, _OUTPUT_KS_FILE_PATTERN_RR_OMCL

# Columns needed to recompute per-family/per-node weights later (via filter_compute_weights)
# for an arbitrary Ks/alignment-quality cutoff, without having to re-read the original TSV file.
_REQUIRED_COLUMNS = ['Family', 'Node', 'Ks', 'AlignmentCoverage', 'AlignmentIdentity', 'AlignmentLength']


def initialize_paralog_db(db_path):
	"""
	Initialize paralog Ks database file if it doesn't exist.
	Creates file with proper column headers.

	:param db_path: path to the paralog Ks database file
	:return: True if the database file exists or was successfully created, False if it could not be created
	         (e.g. parent directory missing, permission denied). Callers are responsible for acting on failure.
	"""
	if os.path.isfile(db_path):
		return True

	logging.info(f"Paralog Ks list database [{db_path}] not found: creating a new one.")
	try:
		with open(db_path, "w+") as outfile:
			outfile.write('\tParanome\tAnchors\tReciprocally_retained\n')
		return True
	except Exception as e:
		logging.error(f"Could not create paralog Ks database at [{db_path}]: {str(e)}")
		return False


def _extract_columns_from_tsv(tsv_path):
	"""
	Read a wgd Ks TSV file and extract the columns needed to recompute weights later
	(Family, Node, Ks, AlignmentCoverage, AlignmentIdentity, AlignmentLength), as a dict of lists.
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
	Write consolidated paralog Ks data to database. Weights are not stored: they are recomputed
	at analysis time from the stored Family/Node/alignment columns for whichever Ks/alignment
	cutoff is requested.

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

	if ks_paranome is not None or ks_anchors is not None or ks_recret is not None:
		logging.info("  - Writing consolidated Ks lists to database")

		# Remove any existing row for this species to avoid duplicates (e.g. from a re-run)
		try:
			with open(db_path, "r") as f:
				db_df = read_csv(f, sep="\t", index_col=0)
			if latin_name in db_df.index:
				db_df = db_df.drop(latin_name)
				with open(db_path, "w") as fw:
					fw.write(db_df.to_csv(sep="\t"))
		except Exception:
			pass

		ks_list_new_row = DataFrame([[ks_paranome, ks_anchors, ks_recret]],
		                             columns=['Paranome', 'Anchors', 'Reciprocally_retained'],
		                             index=[latin_name])
		with open(db_path, "a+") as outfile_ks_list:
			outfile_ks_list.write(ks_list_new_row.to_csv(sep="\t", header=None))
		return True
	else:
		logging.warning(f"  No paralog Ks data could be extracted.")
		return False


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