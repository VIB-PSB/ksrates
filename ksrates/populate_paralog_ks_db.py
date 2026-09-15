import os
import sys
import logging
import glob
from pandas import read_csv
import ksrates.fc_configfile as fcConf
import ksrates.fc_consolidate_paralog_ks as fc_consolidate_paralog_ks
from ksrates.utils import init_logging


def populate_paralog_ks_db_batch(config_dir_path, paralog_distributions_path, db_path, force_overwrite=False, num_gfs=2000):
	"""
	Batch populate the paralog Ks database from config files and a paralog distribution directory.

	Reads all config files in config_dir_path to extract species information (informal and latin names),
	then consolidates paralog Ks data from matching directories in paralog_distributions_path.
	Uses latin names as database keys for consistency across analyses.

	:param config_dir_path: path to directory containing config files for species to consolidate
	:param paralog_distributions_path: path to the paralog_distributions directory (containing wgd_* subdirs)
	:param db_path: path where the consolidated paralog Ks database will be saved
	:param force_overwrite: whether to overwrite species that already exist in the database (default: False)
	:param num_gfs: number of gene families for reciprocally retained (default: 2000)
	"""
	init_logging("Batch populating paralog Ks database from config files and TSV files", logging.INFO)

	trigger_exit = False

	logging.info("Checking configuration files from input directory")

	# Check if config directory exists
	config_dir_path = os.path.abspath(config_dir_path)
	if not os.path.isdir(config_dir_path):
		logging.error(f"Config directory not found: {config_dir_path}")
		trigger_exit = True

	# Check if there are config files in it
	config_files = glob.glob(os.path.join(config_dir_path, "*.txt"))
	if not config_files:
		logging.error(f"No config files found in {config_dir_path}")
		trigger_exit = True
	logging.info(f"Found {len(config_files)} files in {config_dir_path}")

	# Extract species info from config files
	species_info = {}  # informal_name -> latin_name
	species_configs = {}  # informal_name -> config_file_path

	for config_file in sorted(config_files):
		try:
			config = fcConf.Configuration(config_file, "")
			informal_name = config.get_species()
			latin_names = config.get_latin_names()
			latin_name = latin_names.get(informal_name)
			species_info[informal_name] = latin_name
			species_configs[informal_name] = config_file
		except Exception as e:
			logging.error(f"Could not read file [{os.path.basename(config_file)}]: {str(e)}")
			continue
	logging.info("")

	if not species_info:
		logging.error("No valid species found in config files")
		logging.error("Exiting")
		sys.exit(1)
	

	logging.info("Checking paralog distributions from input directory")

	# Check paralog distributions directory exists
	paralog_distributions_path = os.path.abspath(paralog_distributions_path)
	if not os.path.isdir(paralog_distributions_path):
		logging.error(f"Directory not found: {paralog_distributions_path}")
		trigger_exit = True

	# Check if there are any wgd directories in it
	wgd_dirs = glob.glob(os.path.join(paralog_distributions_path, "wgd_*"))
	if not wgd_dirs:
		logging.error(f"No wgd_* subdirectories found in {paralog_distributions_path}")
		trigger_exit = True
	logging.info("")

	logging.info("Checking paralog Ks database from input path")
	# Initialize database
	fc_consolidate_paralog_ks.initialize_paralog_db(db_path)

	# Read existing database and remove any entries for species we're about to process
	existing_species = set()
	species_to_remove = set()
	try:
		with open(db_path, "r") as f:
			db_df = read_csv(f, sep="\t", index_col=0)
			existing_species = set(db_df.index)
			if existing_species:
				logging.info(f"Database already contains {len(existing_species)} species")
				# Find species we're about to process that already exist
				species_to_remove = existing_species.intersection(set(species_info.values()))
				if species_to_remove and force_overwrite:
					logging.info(f"Removing {len(species_to_remove)} existing entries that will be overwritten...")
					db_df = db_df.drop(list(species_to_remove))
					with open(db_path, "w") as fw:
						fw.write(db_df.to_csv(sep="\t"))
					logging.info("Cleaned database, ready to append new entries")
	except Exception:
		pass
	logging.info("")

	logging.info(f"Processing {len(species_info)} species:")

	# Process each species
	successful = []
	failed = []
	skipped = []
	overwritten = []
	data_types_found = {}
	alternative_recret_files = {}

	for informal_name, latin_name in species_info.items():
		# Check if species already exists in database (by latin name)
		will_overwrite = latin_name in existing_species
		if will_overwrite and not force_overwrite:
			logging.warning(f"Skipping [{informal_name}] ({latin_name}): already exists in database (use --force to overwrite)")
			skipped.append(informal_name)
			continue

		if will_overwrite and force_overwrite:
			logging.info(f"Overwriting [{informal_name}] ({latin_name}):")
		else:
			logging.info(f"Adding [{informal_name}] ({latin_name}):")

		try:
			# Get config parameters
			config_file = species_configs[informal_name]
			config = fcConf.Configuration(config_file, "")
			paranome = config.get_paranome()
			colinearity = config.get_colinearity()
			reciprocal_retention = config.get_reciprocal_retention()
			bottom = config.use_bottom_gfs_instead_of_top(reciprocal_retention)
			rank_type = config.get_reciprocal_retention_rank_type(reciprocal_retention)

			# Temporarily change to parent directory for extraction
			original_cwd = os.getcwd()
			parent_dir = os.path.dirname(paralog_distributions_path)
			os.chdir(parent_dir)

			# Extract Ks data from TSV files
			ks_data = fc_consolidate_paralog_ks.extract_paralog_ks_from_tsv(
				informal_name, paranome, colinearity, reciprocal_retention,
				num_gfs=num_gfs, rank_type=rank_type, bottom=bottom
			)

			# Check for alternative recret files
			wgd_dir = os.path.join(os.path.dirname(paralog_distributions_path), os.path.basename(paralog_distributions_path), f"wgd_{informal_name}")
			if reciprocal_retention and os.path.isdir(wgd_dir):
				all_recret_files = glob.glob(os.path.join(wgd_dir, f"{informal_name}.ks_recret_top*.tsv"))
				if all_recret_files:
					used_file = os.path.basename([f for f in all_recret_files if f"{num_gfs}" in f][0]) if any(f"{num_gfs}" in f for f in all_recret_files) else None
					alternative_files = [os.path.basename(f) for f in all_recret_files if f"{num_gfs}" not in f]
					if alternative_files:
						alternative_recret_files[informal_name] = alternative_files

			# Track which data types were found
			found_types = []
			if ks_data['paranome'] is not None:
				found_types.append('paranome')
			if ks_data['anchors'] is not None:
				found_types.append('anchors')
			if ks_data['recret'] is not None:
				found_types.append('recret')

			data_types_found[informal_name] = found_types

			# Write to database using latin name
			success = fc_consolidate_paralog_ks.write_to_paralog_db(latin_name, ks_data, db_path)

			if success:
				if will_overwrite:
					overwritten.append(informal_name)
				else:
					successful.append(informal_name)
			else:
				logging.warning(f"  No data extracted for [{informal_name}] ({latin_name})")
				failed.append(informal_name)

			os.chdir(original_cwd)

		except Exception as e:
			logging.error(f"  Error processing [{informal_name}] ({latin_name}): {str(e)}")
			failed.append(informal_name)
			try:
				os.chdir(original_cwd)
			except:
				pass
	logging.info("")

	# Summary
	logging.info("=" * 70)
	logging.info("BATCH PROCESSING SUMMARY")
	logging.info("=" * 70)
	logging.info(f"Processed: {len(successful)} new, {len(overwritten)} overwritten, {len(skipped)} skipped, {len(failed)} failed, out of {len(species_info)} species")

	if successful:
		logging.info("")
		logging.info("Successfully added new species:")
		logging.info("  Species                Paranome  Anchors   RecRet")
		logging.info("  " + "-" * 64)
		for informal_name in successful:
			types_found = data_types_found.get(informal_name, [])
			paranome_status = "[YES]" if 'paranome' in types_found else "[NO]"
			anchors_status = "[YES]" if 'anchors' in types_found else "[NO]"
			recret_status = "[YES]" if 'recret' in types_found else "[NO]"
			logging.info(f"  {informal_name:20}  {paranome_status:8}  {anchors_status:8}  {recret_status:8}")

	if overwritten:
		logging.info("")
		logging.info("Successfully overwritten existing species:")
		logging.info("  Species                Paranome  Anchors   RecRet")
		logging.info("  " + "-" * 64)
		for informal_name in overwritten:
			types_found = data_types_found.get(informal_name, [])
			paranome_status = "[YES]" if 'paranome' in types_found else "[NO]"
			anchors_status = "[YES]" if 'anchors' in types_found else "[NO]"
			recret_status = "[YES]" if 'recret' in types_found else "[NO]"
			logging.info(f"  {informal_name:20}  {paranome_status:8}  {anchors_status:8}  {recret_status:8}")

	if skipped:
		logging.info("")
		logging.warning(f"Skipped (already in database, use --force to overwrite): {len(skipped)} species")
		for informal_name in skipped:
			logging.warning(f"  - {informal_name}")

	if failed:
		logging.info("")
		logging.error(f"Failed to process: {len(failed)} species")
		for informal_name in failed:
			logging.error(f"  - {informal_name}")

	if alternative_recret_files:
		logging.info("")
		logging.warning("=" * 70)
		logging.warning("Alternative reciprocally retained (recret) files available:")
		logging.warning("=" * 70)
		for informal_name, alt_files in alternative_recret_files.items():
			logging.warning(f"  [{informal_name}]:")
			for alt_file in alt_files:
				gf_num = alt_file.split("_")[-1].replace(".tsv", "").replace("top", "")
				logging.warning(f"    - {alt_file}")
				logging.warning(f"      To consolidate this, add its config to the config directory and use --num-gfs {gf_num} --force")

	logging.info("")
	logging.info(f"Database saved to: {db_path}")
	logging.info("=" * 70)
	logging.info("")
	logging.info("Done")
