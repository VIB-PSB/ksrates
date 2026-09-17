import os
import sys
import logging
import glob
import ksrates.fc_configfile as fcConf
import ksrates.fc_consolidate_paralog_ks as fc_consolidate_paralog_ks
from ksrates.utils import init_logging


def populate_paralog_ks_db_batch(config_dir_path, paralog_distributions_path, db_path, force_overwrite=False, num_gfs=2000):
	"""
	Batch populate the paralog Ks database from config files and a paralog distribution directory.

	Reads all config files in config_dir_path to extract species information (informal and latin names),
	then consolidates paralog Ks data from matching directories in paralog_distributions_path.
	Attempts to consolidate all types of paralog Ks data (paranome, anchors and reciprocally retained),
	throwing warnings if any of these files are not found (NOTE: possibly they were just never generated).
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
	logging.info("")

	logging.info(f"Processing {len(species_info)} species:")

	# Process each species
	successful = []
	updated = []
	failed = []
	skipped = []
	overwritten = []
	data_types_found = {}
	iadhore_found = {}
	alternative_recret_files = {}

	for informal_name, latin_name in species_info.items():
		# Check if species already exists in database (by latin name), and if so which analysis
		# types/i-ADHoRe bundle it already has data for (e.g. paranome and anchors were stored on
		# a first pass, and reciprocal retention only finished days later): without --force, only
		# what's actually missing gets extracted and written, so an already-populated type is
		# never touched, but a species isn't skipped wholesale just because it exists.
		exists = fc_consolidate_paralog_ks.species_exists(db_path, latin_name)
		if exists and not force_overwrite:
			existing_types = fc_consolidate_paralog_ks.existing_data_types(db_path, latin_name)
		else:
			existing_types = {"paranome": False, "anchors": False, "recret": False, "iadhore": False}

		paranome_enabled = force_overwrite or not existing_types["paranome"]
		anchors_enabled = force_overwrite or not existing_types["anchors"]
		recret_enabled = force_overwrite or not existing_types["recret"]
		iadhore_enabled = force_overwrite or not existing_types["iadhore"]

		if exists and not force_overwrite and not (paranome_enabled or anchors_enabled or recret_enabled or iadhore_enabled):
			logging.warning(f"Skipping [{informal_name}] ({latin_name}): already fully populated (use --force to overwrite)")
			skipped.append(informal_name)
			continue

		if force_overwrite and exists:
			logging.info(f"Overwriting [{informal_name}] ({latin_name}):")
		elif exists:
			missing = [t for t, enabled in [('paranome', paranome_enabled), ('anchors', anchors_enabled), ('recret', recret_enabled)] if enabled] + (['i-ADHoRe files'] if iadhore_enabled else [])
			logging.info(f"Filling in missing data for [{informal_name}] ({latin_name}): {', '.join(missing)}")
		else:
			logging.info(f"Adding [{informal_name}] ({latin_name})")

		try:
			# Get config parameters. Note: we don't trust the config's current paranome/colinearity/
			# reciprocal_retention yes/no toggles, since they may have changed since the TSV files were
			# generated (e.g. a species run with recret=yes in the past, but config now set to recret=no).
			# Instead, we always search among the enabled types for whichever actually exist.
			config_file = species_configs[informal_name]
			config = fcConf.Configuration(config_file, "")
			bottom = config.use_bottom_gfs_instead_of_top(True)
			rank_type = config.get_reciprocal_retention_rank_type(True)

			# Temporarily change to parent directory for extraction
			original_cwd = os.getcwd()
			parent_dir = os.path.dirname(paralog_distributions_path)
			os.chdir(parent_dir)

			# Extract Ks data from TSV files, but only for the types that are actually enabled
			# above (missing, or forced); a type already present in the database is left untouched.
			ks_data = fc_consolidate_paralog_ks.extract_paralog_ks_from_tsv(
				informal_name, paranome_enabled=paranome_enabled, anchors_enabled=anchors_enabled,
				reciprocal_retention_enabled=recret_enabled, num_gfs=num_gfs, rank_type=rank_type, bottom=bottom
			)

			# Check for alternative recret files
			wgd_dir = os.path.join(os.path.dirname(paralog_distributions_path), os.path.basename(paralog_distributions_path), f"wgd_{informal_name}")
			if os.path.isdir(wgd_dir):
				all_recret_files = glob.glob(os.path.join(wgd_dir, f"{informal_name}.ks_recret_top*.tsv"))
				if all_recret_files:
					used_file = os.path.basename([f for f in all_recret_files if f"{num_gfs}" in f][0]) if any(f"{num_gfs}" in f for f in all_recret_files) else None
					alternative_files = [os.path.basename(f) for f in all_recret_files if f"{num_gfs}" not in f]
					if alternative_files:
						alternative_recret_files[informal_name] = alternative_files

			# Write to database using latin name (a no-op, returning False, if nothing was enabled
			# above or nothing was actually found on disk for the enabled types)
			wrote_ks_data = fc_consolidate_paralog_ks.write_to_paralog_db(latin_name, ks_data, db_path)

			# Also pick up the i-ADHoRe output files if present on disk and enabled above, same as
			# wgd_paralogs.py does right after a fresh colinearity run: extract_paralog_ks_from_tsv()
			# above only ever reads the three Ks TSVs, so without this the accessory files a past
			# colinearity run already produced would never make it into the database via this batch
			# command.
			if iadhore_enabled:
				iadhore_texts = fc_consolidate_paralog_ks.extract_anchor_iadhore_files(informal_name)
				wrote_iadhore = any(text is not None for text in iadhore_texts.values())
				if wrote_iadhore:
					fc_consolidate_paralog_ks.write_anchor_iadhore_files(latin_name, iadhore_texts, db_path)
			else:
				wrote_iadhore = False

			# Final status shown in the summary table below: the union of what the species already
			# had before this run and whatever was freshly written just now (rather than only what
			# this run itself touched), so an already-populated type isn't shown as missing.
			data_types_found[informal_name] = [
				t for t, already, fresh in [
					('paranome', existing_types["paranome"], ks_data['paranome'] is not None),
					('anchors', existing_types["anchors"], ks_data['anchors'] is not None),
					('recret', existing_types["recret"], ks_data['recret'] is not None),
				] if already or fresh
			]
			iadhore_found[informal_name] = existing_types["iadhore"] or wrote_iadhore

			if wrote_ks_data or wrote_iadhore:
				if force_overwrite and exists:
					overwritten.append(informal_name)
				elif exists:
					updated.append(informal_name)
				else:
					successful.append(informal_name)
			elif exists:
				# Something was missing, but the TSV/i-ADHoRe files for it still aren't on disk
				# yet (e.g. reciprocal retention hasn't finished) - nothing new to add this run.
				logging.warning(f"  Still missing data for [{informal_name}] ({latin_name}); nothing new found on disk")
				skipped.append(informal_name)
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
	logging.info(f"Processed: {len(successful)} new, {len(updated)} updated, {len(overwritten)} overwritten, "
	             f"{len(skipped)} skipped, {len(failed)} failed, out of {len(species_info)} species")

	spacer_informal_names = max(len(s) for s in species_info.keys()) + 5
	spacer_latin_names = max(len(s) for s in species_info.values()) + 3

	def print_species_table(title, informal_names):
		logging.info("")
		logging.info(title)
		logging.info(f"  {'Species':{spacer_informal_names}} {'Latin name':{spacer_latin_names}} {'Paranome':9} {'Anchors':9} {'RecRet':9} {'iADHoRe':9}")
		logging.info("  " + "-" * 64)
		for informal_name in informal_names:
			types_found = data_types_found.get(informal_name, [])
			paranome_status = "[YES]" if 'paranome' in types_found else "[NO]"
			anchors_status = "[YES]" if 'anchors' in types_found else "[NO]"
			recret_status = "[YES]" if 'recret' in types_found else "[NO]"
			iadhore_status = "[YES]" if iadhore_found.get(informal_name) else "[NO]"
			logging.info(f"  {informal_name:{spacer_informal_names}} {species_info[informal_name]:{spacer_latin_names}} {paranome_status:9}  {anchors_status:9}  {recret_status:9}  {iadhore_status:9}")

	if successful:
		print_species_table("Successfully added new species:", successful)

	if updated:
		print_species_table("Filled in missing data for existing species:", updated)

	if overwritten:
		print_species_table("Successfully overwritten existing species:", overwritten)

	if skipped:
		logging.info("")
		logging.warning(f"Skipped: {len(skipped)} species")
		for informal_name in skipped:
			logging.warning(f"  - {informal_name} ({species_info.get(informal_name)})")

	if failed:
		logging.info("")
		logging.error(f"Failed to process: {len(failed)} species")
		for informal_name in failed:
			logging.error(f"  - {informal_name} ({species_info.get(informal_name)})")

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
