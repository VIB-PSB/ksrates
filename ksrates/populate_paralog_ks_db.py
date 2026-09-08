import os
import sys
import logging
import glob
import ksrates.fc_configfile as fcConf
import ksrates.fc_consolidate_paralog_ks as fc_consolidate_paralog_ks
from ksrates.utils import init_logging


def populate_paralog_ks_db(config_file, expert_config_file):
	"""
	Populates the paralog Ks database from existing TSV files.
	This utility allows users to consolidate previously generated paralog Ks TSV files
	into the new unified database format.

	:param config_file: configuration file
	:param expert_config_file: expert configuration file (can be empty string)
	"""
	config = fcConf.Configuration(config_file, expert_config_file)
	species = config.get_species()
	init_logging("Populating paralog Ks database from TSV files", config.get_logging_level())
	logging.info("Loading parameters and input files")

	paranome = config.get_paranome()
	colinearity = config.get_colinearity()
	reciprocal_retention = config.get_reciprocal_retention()
	ks_list_paralog_db_path = config.get_paralog_ks_db()

	num_gfs = config.get_num_reciprocal_retention_gfs(reciprocal_retention)
	bottom = config.use_bottom_gfs_instead_of_top(reciprocal_retention)
	rank_type = config.get_reciprocal_retention_rank_type(reciprocal_retention)

	if not paranome and not colinearity and not reciprocal_retention:
		logging.error('At least one of the "paranome", "collinearity" or "reciprocal_retention" parameters in the configuration file needs to be set to "yes".')
		logging.error("Exiting.")
		sys.exit(1)

	logging.info("")
	logging.info(f"Populating paralog Ks database for species [{species}]")

	# Initialize database file if it doesn't exist
	fc_consolidate_paralog_ks.initialize_paralog_db(ks_list_paralog_db_path)

	# Extract Ks data from TSV files
	logging.info("- Extracting paralog Ks data from TSV files")
	ks_data = fc_consolidate_paralog_ks.extract_paralog_ks_from_tsv(species, paranome, colinearity, reciprocal_retention,
	                                                                   num_gfs, rank_type, bottom)

	# Write to database
	logging.info("")
	success = fc_consolidate_paralog_ks.write_to_paralog_db(species, ks_data, ks_list_paralog_db_path)

	if success:
		logging.info(f"Successfully populated database for [{species}]")
		logging.info("")
		logging.info("Done")
	else:
		logging.error(f"No paralog Ks data could be extracted for {species}.")
		logging.error("Exiting")
		sys.exit(1)


def populate_paralog_ks_db_batch(paralog_distributions_path, db_path, species_list=None, force_overwrite=False, num_gfs=2000, paranome=True, anchors=True, reciprocal_retention=True):
	"""
	Batch populate the paralog Ks database from species found in a paralog_distributions directory.
	Discovers all wgd_<species> subdirectories and processes their Ks TSV files.
	If species_list is provided, only processes those species; otherwise processes all found.

	:param paralog_distributions_path: path to the paralog_distributions directory (or similar directory containing wgd_* subdirs)
	:param db_path: path where the consolidated paralog Ks database will be saved
	:param species_list: optional list of species names to process (if None, processes all found in directory)
	:param force_overwrite: whether to overwrite species that already exist in the database (default: False)
	:param num_gfs: number of gene families for reciprocally retained (default: 2000)
	:param paranome: whether to extract paranome Ks values (default: True)
	:param anchors: whether to extract anchor pair Ks values (default: True)
	:param reciprocal_retention: whether to extract reciprocally retained Ks values (default: True)
	"""
	init_logging("Batch populating paralog Ks database from TSV files", logging.INFO)
	logging.info("Loading paralog distributions from directory")
	logging.info("")

	# Discover all wgd_* subdirectories
	paralog_distributions_path = os.path.abspath(paralog_distributions_path)
	if not os.path.isdir(paralog_distributions_path):
		logging.error(f"Directory not found: {paralog_distributions_path}")
		logging.error("Exiting")
		sys.exit(1)

	wgd_dirs = glob.glob(os.path.join(paralog_distributions_path, "wgd_*"))
	if not wgd_dirs:
		logging.error(f"No wgd_* subdirectories found in {paralog_distributions_path}")
		logging.error("Exiting")
		sys.exit(1)

	# Extract species names from wgd_* directory names
	all_species_found = []
	for wgd_dir in sorted(wgd_dirs):
		species_name = os.path.basename(wgd_dir).replace("wgd_", "")
		all_species_found.append(species_name)

	# Use provided species list or all found species
	if species_list:
		# Verify that requested species exist in the directory
		missing_species = [sp for sp in species_list if sp not in all_species_found]
		if missing_species:
			logging.error(f"The following requested species were not found in {paralog_distributions_path}:")
			for sp in missing_species:
				logging.error(f"  - {sp}")
			logging.error("Exiting")
			sys.exit(1)
		species_to_process = species_list
	else:
		species_to_process = all_species_found

	logging.info(f"Found {len(all_species_found)} species total in {paralog_distributions_path}")
	logging.info(f"Processing {len(species_to_process)} species:")
	for sp in species_to_process:
		logging.info(f"  - {sp}")
	logging.info("")

	# Initialize database
	fc_consolidate_paralog_ks.initialize_paralog_db(db_path)

	# Read existing database to check which species are already present
	import pandas
	from ast import literal_eval
	existing_species = set()
	try:
		with open(db_path, "r") as f:
			db_df = pandas.read_csv(f, sep="\t", index_col=0)
			existing_species = set(db_df.index)
			if existing_species:
				logging.info(f"Database already contains {len(existing_species)} species")
	except Exception:
		pass

	# Process each species
	successful = []
	failed = []
	skipped = []
	overwritten = []
	data_types_found = {}
	alternative_recret_files = {}  # Track species with alternative recret files

	for species in species_to_process:
		# Check if species already exists in database
		will_overwrite = species in existing_species
		if will_overwrite and not force_overwrite:
			logging.warning(f"Skipping [{species}]: already exists in database (use --force to overwrite)")
			skipped.append(species)
			continue

		if will_overwrite and force_overwrite:
			logging.info(f"Processing [{species}]: will overwrite existing data")
		else:
			logging.info(f"Processing [{species}]:")

		try:
			# Temporarily change paralog_distributions directory for extraction
			original_cwd = os.getcwd()
			parent_dir = os.path.dirname(paralog_distributions_path)
			os.chdir(parent_dir)

			# Extract Ks data from TSV files
			ks_data = fc_consolidate_paralog_ks.extract_paralog_ks_from_tsv(
				species, paranome, anchors, reciprocal_retention,
				num_gfs=num_gfs, rank_type="lambda", bottom=False
			)

			# Check for alternative recret files
			wgd_dir = os.path.join(os.path.dirname(paralog_distributions_path), os.path.basename(paralog_distributions_path), f"wgd_{species}")
			if reciprocal_retention and os.path.isdir(wgd_dir):
				all_recret_files = glob.glob(os.path.join(wgd_dir, f"{species}.ks_recret_top*.tsv"))
				if all_recret_files:
					used_file = os.path.basename([f for f in all_recret_files if f"{num_gfs}" in f][0]) if any(f"{num_gfs}" in f for f in all_recret_files) else None
					alternative_files = [os.path.basename(f) for f in all_recret_files if f"{num_gfs}" not in f]
					if alternative_files:
						alternative_recret_files[species] = alternative_files

			# Track which data types were found
			found_types = []
			if ks_data['paranome'] is not None:
				found_types.append('paranome')
			if ks_data['anchors'] is not None:
				found_types.append('anchors')
			if ks_data['recret'] is not None:
				found_types.append('recret')

			data_types_found[species] = found_types

			# Write to database
			success = fc_consolidate_paralog_ks.write_to_paralog_db(species, ks_data, db_path)

			if success:
				if will_overwrite:
					logging.info(f"  Successfully overwrote database entry for [{species}]")
					overwritten.append(species)
				else:
					logging.info(f"  Successfully populated database for [{species}]")
					successful.append(species)
			else:
				logging.warning(f"  No data extracted for [{species}]")
				failed.append(species)

			os.chdir(original_cwd)
			logging.info("")

		except Exception as e:
			logging.error(f"  Error processing [{species}]: {str(e)}")
			failed.append(species)
			os.chdir(original_cwd)
			logging.info("")

	# Summary
	logging.info("=" * 70)
	logging.info("BATCH PROCESSING SUMMARY")
	logging.info("=" * 70)
	logging.info(f"Processed: {len(successful)} new, {len(overwritten)} overwritten, {len(skipped)} skipped, {len(failed)} failed out of {len(species_to_process)} species")

	if successful:
		logging.info("")
		logging.info("Successfully added new species:")
		logging.info("  Species                Paranome  Anchors   RecRet")
		logging.info("  " + "-" * 64)
		for sp in successful:
			types_found = data_types_found.get(sp, [])
			paranome_status = "[YES]" if 'paranome' in types_found else "[NO]"
			anchors_status = "[YES]" if 'anchors' in types_found else "[NO]"
			recret_status = "[YES]" if 'recret' in types_found else "[NO]"
			logging.info(f"  {sp:20}  {paranome_status:8}  {anchors_status:8}  {recret_status:8}")

	if overwritten:
		logging.info("")
		logging.info("Successfully overwritten existing species:")
		logging.info("  Species                Paranome  Anchors   RecRet")
		logging.info("  " + "-" * 64)
		for sp in overwritten:
			types_found = data_types_found.get(sp, [])
			paranome_status = "[YES]" if 'paranome' in types_found else "[NO]"
			anchors_status = "[YES]" if 'anchors' in types_found else "[NO]"
			recret_status = "[YES]" if 'recret' in types_found else "[NO]"
			logging.info(f"  {sp:20}  {paranome_status:8}  {anchors_status:8}  {recret_status:8}")

	if skipped:
		logging.info("")
		logging.warning(f"Skipped (already in database, use --force to overwrite): {len(skipped)} species")
		for sp in skipped:
			logging.warning(f"  - {sp}")

	if failed:
		logging.info("")
		logging.error(f"Failed to process: {len(failed)} species")
		for sp in failed:
			logging.error(f"  - {sp}")

	if alternative_recret_files:
		logging.info("")
		logging.warning("=" * 70)
		logging.warning("Alternative reciprocally retained (recret) files available:")
		logging.warning("=" * 70)
		for species, alt_files in alternative_recret_files.items():
			logging.warning(f"  [{species}]:")
			for alt_file in alt_files:
				# Extract the number from filename (e.g., "top4000" from "species.ks_recret_top4000.tsv")
				gf_num = alt_file.split("_")[-1].replace(".tsv", "").replace("top", "")
				logging.warning(f"    - {alt_file}")
				logging.warning(f"      To consolidate this, run: ksrates populate-paralog-ks-db <dir> {species} --database {db_path} --num-gfs {gf_num} --force")

	logging.info("")
	logging.info(f"Database saved to: {db_path}")
	logging.info("=" * 70)
	logging.info("")
	logging.info("Done")
