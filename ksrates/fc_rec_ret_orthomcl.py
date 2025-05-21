import os
import sys
import re
import io
import logging
import tarfile
from pandas import DataFrame, concat, Series, read_csv
from numpy import nan
import subprocess
from Bio import SeqIO
from ksrates.utils import translate_cds
from wgd_ksrates.utils import read_fasta
import shutil

original_37_names =     {'alyr37':'Arabidopsis lyrata',   'atha37':'Arabidopsis thaliana',  'atri37':'Amborella trichopoda',  'bdis37':'Brachypodium distachyon', \
                         'brap37':'Brassica rapa',        'cari37':'Cicer arietinum',       'ccaj37':'Cajanus cajan',         'clan37':'Citrullus lanatus', \
                         'cmel37':'Cucumis melo',         'cpap37':'Carica papaya',         'crub37':'Capsella rubella',      'csat37':'Cucumis sativus', \
                         'fves37':'Fragaria vesca',       'gmax37':'Glycine maxima',        'grai37':'Gossypium raimondii',   'hvul37':'Hordeum vulgare', \
                         'jcur37':'Jatropha curcas',      'ljap37':'Lotus japonicus',       'lusi37':'Linum usitatissiumum',  'macu37':'Musa acuminata', \
                         'mesc37':'Manihot esculenta',    'mtru37':'Medicago truncatula',   'osat37':'Oryza sativa',          'pbre37':'Pyrus bretschneideri', \
                         'pdac37':'Phoenix dactylifera',  'pmum37':'Prunus mume',           'pper37':'Prunus persica',        'ptri37':'Populus trichocarpa', \
                         'rcom37':'Ricinus communis',     'sbic37':'Sorghum bicolor',       'sita37':'Setaria italica',       'slyc37':'Solanum lycopersicum', \
                         'stub37':'Solanum tuberosum',    'tcac37':'Theobroma cacao',       'tpar37':'Thellungiella parvula', 'vvin37': 'Vitis vinifera', \
                         'zmay37':'Zea mays'}


def check_recret_name_compatibility(informal_species_name_list, reciprocal_retention):
    """
    Checks if the informal species' name(s) provided in the list are equal to any ID of the original 37 angiosperms species
    used by Li (2016) and Tasdighian (2017) to generate the reciprocally retained GF ranking.
    If recret has been asked and there is at least an overlapping name, then the pipeline exists due
    to incompatibility with the implemented OrthoMCL workflow to get recret GFs for focal species. 
    Note: This is very unlikely to happen, since the original species IDs are e.g. "atha37", "zmay37"...

    :param informal_species_name_list: list of informal species names (can be all species in Newick tree or just focal species)
    :param reciprocal_retention: boolean stating if reciprocal retention analysis has been requested (True)
    """
    trigger_exit = False
    species_same_name_with_37 = []
    for informal_name in informal_species_name_list:
        if informal_name in original_37_names.keys():
            species_same_name_with_37.append(informal_name)
                
    if len(species_same_name_with_37) != 0:
        logging.info("Checking species name compatibility with reciprocal retention pipeline")
        # First state if recret analysis has been requested or not from configuration
        if reciprocal_retention:
            trigger_exit = True
            logging.error("Since reciprocal retention analysis has been requested in the analysis configuration,")
        else: 
            logging.warning("Currently reciprocal retention analysis has not been requested;")
            logging.warning("however, should you want to execute it in the future,")
        # Then explain that having the same name is a problem or might be so in the future
        common_text_part = []
        common_text_part.append("the informal species' names provided in the Newick tree must be different from")
        common_text_part.append("those internally assigned to the original 37 angiosperms in the OrthoMCL workflow.")
        common_text_part.append("These latter names are 4-letter words containing the first genus letter and")
        common_text_part.append("the first three species letters, plus number '37': e.g. 'atha37' for Arabidopsis thaliana.")
        common_text_part.append("For more details, please refer to Methods.")
        # Then either require or suggest to make the name changes
        if reciprocal_retention:
            for line in common_text_part:
                logging.error(line)
            logging.error("The following informal species names in the configuration file need to be changed:")
        else:
            for line in common_text_part:
                logging.warning(line)
            logging.warning("Consider changing the following conflicting names in the configuration file to prevent future compatibility issues:")
        # Finally list the species names whose name is problematic
        list_of_problematic_names = []
        for informal_name in species_same_name_with_37:
            list_of_problematic_names.append(f"- {informal_name}")
        if reciprocal_retention:
            for line in list_of_problematic_names:
                logging.error(line)
        else:
            for line in list_of_problematic_names:
                logging.warning(line)
                logging.warning("")
    return trigger_exit

####################  MAKE ORTHOGROUPS FROM SCRATCH ####################

def merge_original_and_new_fastas(all_fasta_filepaths, merged_fasta_outpath):
    """
    Merges the FASTA sequences of the focal species together with the original FASTAs 
    of the 37 angiosperms used to generate the ranking in Li et al (2016) and
    Tasdighian et al (2017). The output will be used by diamond.
    
    :param all_fasta_filepaths: path to directory containing original FASTAs of 37 angiosperms
    :param merged_fasta_outpath: path to the file with all merged sequences
    """        
    command = ["cat"]
    command.extend(all_fasta_filepaths)
    logging.debug(' '.join(command))
    logging.debug("")
    
    try:
        with open(merged_fasta_outpath, 'w') as merged_out:
            sp = subprocess.run(command, stdout=merged_out, stderr=subprocess.PIPE, encoding='utf-8', check=True)
    
    except subprocess.CalledProcessError as e:
        logging.error(f"Merging focal species FASTA with original 37 species FASTAs failed with return code: {e.returncode}")
        if e.stderr:
            logging.error("Merging FASTA standard error output:\n" + e.stderr.strip())
    
        logging.error("Exiting.")
        sys.exit(e.returncode)
    logging.debug("")
    return


def run_diamond_makedb(merged_fasta_file, database_filepath, n_threads):
    """
    Generates the diamond database with FASTAs from 37 species and focal species.
    
    :param merged_fasta_file: merged file from 37-species original fastas plus focal species fasta
    :param database_filepath: output path for diamond database
    :param n_threads: number of threads to use
    """   
    # Define makedatabase command
    command = ["diamond", "makedb", "--in", merged_fasta_file, "-d", database_filepath, "--threads", str(n_threads)]
    logging.info(' '.join(command))
    logging.info("")
    
    try:
        sp = subprocess.run(command, stdout=subprocess.PIPE, stderr=subprocess.PIPE, encoding='utf-8', check=True)
        out = sp.stdout.strip()
        err = sp.stderr.strip()
        if out:
            logging.info("diamond makedb output:\n" + out)
        if err:
            logging.error("diamond makedb standard error output:\n" + err)
    
    except subprocess.CalledProcessError as e:
        logging.error(f"diamond makedb execution failed with return code: {e.returncode}")
        if e.stderr:
            logging.error("diamond makedb standard error output:\n" + e.stderr.strip())
            
        if os.path.exists(database_filepath):
            logging.error(f"Removing faulty database file [{database_filepath}]:")
            os.remove(database_filepath)
        
        logging.error("Exiting.")
        sys.exit(e.returncode)
    return
    
def run_diamond_search(database_filepath, merged_fasta_file, dmd_homology_table_filepath, 
                       n_threads, max_target_seqs=750):
    """
    Performs a diamond blastp all-vs.-all search within a file containing all sequences coming
    from the original FASTAs of the 37 angiosperms used to generate the ranking, plus the sequences
    of the FASTA file of the focal species.
    NOTE: Ideally, I should be able to do that for all species in the user dataset, as a separate command?
    NOTE: add param evalue cutoff?
    
    :param database_filepath: diamond database of the merged FASTAs
    :param merged_fasta_file: merged file from 37-species original fastas plus focal species fasta
    :param dmd_homology_table_filepath: output homology table filepath
    :param n_threads: number of threads to use for blastp
    :pram max_target_seqs: max number of target sequences per query to report alignments for [default: 750]
    """
    # Define homology search command    
    command = ["diamond", "blastp", '-d', database_filepath, '-q', merged_fasta_file, 
               "-o", dmd_homology_table_filepath,
               "--threads", str(n_threads), "--max-target-seqs", str(max_target_seqs)]
    logging.info(' '.join(command))
    logging.info("")
    
    try:
        sp = subprocess.run(command, stdout=subprocess.PIPE, stderr=subprocess.PIPE, encoding='utf-8', check=True)
        out = sp.stdout.strip()
        err = sp.stderr.strip()
        if out:
            logging.info("diamond blastp output:\n" + out)
        if err: # NOTE: stdErr captures the actual log output
            logging.error("diamond blastp standard error output:\n" + err)
            
    except subprocess.CalledProcessError as e:
        logging.error(f"diamond blastp execution failed with return code: {e.returncode}")
        if e.stderr: # NOTE: stdErr captures the actual log output
            logging.error("diamond blastp standard error output:\n" + e.stderr.strip())
        if e.returncode == -11:
            logging.error("This looks like a segmentation fault, you may want to look for and "
                        "remove any core dump files (core.*).")
            logging.error("Try to increase the memory and/or lower the number of threads.")
        
        # Remove output file as it is likely incomplete or broken
        if os.path.exists(dmd_homology_table_filepath):
            os.remove(dmd_homology_table_filepath)
        logging.error("Exiting.")
        
        sys.exit(e.returncode)
    return


def make_all_gg_file_for_orthomcl(species_name, focal_species_aa_filepath, 
                                  original_37_fasta_dir_targz, output_all_gg_file):
    """
    Generates the all.gg file needed by OrthoMCL.
    Its format is one line per FASTA file showing the basename of the file, a column and 
    then the space-separated list of all sequence identifiers.
    E.g. one of the merged file was "tpar.fasta":
        tpar37: Tp1g33280 Tp3g00890 ...

    :param focal_species_aa_filepath: path to the focal species FASTA file
    :param original_37_fasta_dir_targz: path to the original 37 species FASTA file (tar.gz)
    :param output_all_gg_file: output path for the all.gg file
    """
    # Open 37 angiosperm tarfile 
    tar = tarfile.open(original_37_fasta_dir_targz)
    
    # Open the output all.gg file
    with open(output_all_gg_file, 'w') as all_gg_out:
        
        # First read the gene IDs of focal species' FASTA file
        logging.debug("Adding focal species' gene IDs")
        with open(focal_species_aa_filepath, "r") as focal_species_aa:
            # Read the FASTA file and extract gene IDs
            gene_ids = [record.id for record in SeqIO.parse(focal_species_aa, 'fasta')]
            # Write the species_name, followed by the gene IDs, to the output all.gg file
            # Use species_name as focal species' identifier
            all_gg_out.write(species_name + ': ' + ' '.join(gene_ids) + '\n')
        
        # Then do the same for each of the 37 angiosperms FASTAs
        logging.debug("Adding original 37 species' gene IDs")
        for original_fasta in tar.getmembers():
            content_original_fasta = tar.extractfile(original_fasta)
            # Get the basename of the FASTA file as species identifier
            basename_with_suffix = original_fasta.name
            basename = os.path.splitext(basename_with_suffix)[0]
            # Decode the binary content to text
            content_original_fasta_textobj = io.TextIOWrapper(content_original_fasta, encoding="utf-8")
            # Read the FASTA file and extract gene IDs
            gene_ids = [record.id for record in SeqIO.parse(content_original_fasta_textobj, 'fasta')]
            # Write the basename, followed by the gene IDs, to the output all.gg file
            all_gg_out.write(basename + ': ' + ' '.join(gene_ids) + '\n')
    return output_all_gg_file


def run_orthomcl_command(homology_search_table, parsed_homology_table, 
                         output_all_gg_file, orthomcl_outdir_focal, inflation=1.5,
                         use_original_orthomcl_version=False):
    """
    Run either OrthoMCL mode 3 or 4. Mode 3 accepts the homology table and the gg file, then parses the table
    (generating the BPO file) and finally runs MCL. Mode 4 accepts instead an already existing parsed homology table 
    (BPO) and runs directly MCL. Mode 4 is useful when BPO has already been generated (by a previous run, or 
    independentently through stand-alone OrthoMCL).
    
    :param homology_search_table: BLAST- or diamond-like table with homolog pairs (default format?)
    :param parsed_homology_table: parsed version of homology table ("BPO file") generated by a previous OrthoMCL run
    :param output_all_gg_file: path to all.gg, listing all gene IDs of each FASTA files (one FASTA per row)
    :param inflation: OrthoMCL parameter for orthogroup granulosity [default: 1.5]
    :return: filepath to OrthoMCL output file with orthogroups
    """
    # Set the following environmental variables to guarantee that output ends up in dedicated dir
    os.environ["ORTHOMCL_WORKING_DIR"] = os.path.join(orthomcl_outdir_focal, "")
    os.environ["ORTHOMCL_DATA_DIR"] = os.path.join(orthomcl_outdir_focal, "")
    
    # Define OrthoMCL executable
    if not use_original_orthomcl_version: # Edited version OrthoMCLight (much faster)
        orthomcl_executable = "orthomclight.pl"
    elif use_original_orthomcl_version:   # Original version OrthoMCL v1.4 (much slower)
        orthomcl_executable = "orthomcl.pl"
    
    # If homology table provided, then parse the table and run mcl (OrthoMCL mode: 3)
    if homology_search_table != "":
        command = [orthomcl_executable, "--mode", "3", \
                                "--blast_file", homology_search_table, \
                                "--gg_file", output_all_gg_file, \
                                "--inflation", str(inflation)]
    # Else if user provided existing parsed homology table, skip parsing and run directly MCL (OrthoMCL mode: 4)
    elif parsed_homology_table != "":
        command = [orthomcl_executable, "--mode", "4", \
                                "--bpo_file", parsed_homology_table, \
                                "--gg_file", output_all_gg_file, \
                                "--inflation", str(inflation)]
    logging.info(' '.join(command))
    logging.info("")
    
    try:
        p = subprocess.Popen(command, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, encoding='utf-8', bufsize=1)
        stdout = []
        
        # Print the output lines until the name of OrthoMCL outpur dir is printed; then save it.
        # Find OrthoMCL output dir name (depends on month, date and whether there are duplicates)
        # E.g. Sep_17 or Sep_18_2
        found_outdir_name = False
        while found_outdir_name == False:
            line = p.stdout.readline().strip()
            stdout.append(line)
            logging.info(line)
            if line == "### WORKING DIRECTORY: ###": # First line
                next_line_with_outdir = next(p.stdout, None).strip() # Second line
                if next_line_with_outdir:
                    stdout.append(next_line_with_outdir)
                    logging.info(next_line_with_outdir)
                    # Get clean outdir-path and name and exit this first while loop
                    orthomcl_dirpath = next_line_with_outdir.rstrip("/")
                    match = re.search(r'/(\w+_\d+(?:_\d+)?)\/', next_line_with_outdir)
                    if match:
                        orthomcl_dirname = match.group(1)
                        found_outdir_name = True
                    else:
                        logging.error("The name of OrthoMCL output directory can't be found in the OrthoMCL output.")
                else:
                    logging.error("Unexpected OrthoMCL output log: output directory name not found at the second line.")

        # Continue printing the output log lines until the OrthoMCL is done
        while True:
            line = p.stdout.readline().strip()
            stdout.append(line)
            logging.info(line)
            # If the line is empty and the program finished running, then break
            if line == '' and p.poll() != None:
                # poll is None when program still running; it's the exit code (0) when finished
                break
        
    except: # subprocess.CalledProcessError as e:
        logging.error("An exception has occured. OrthoMCL run has been stopped.")
        logging.error("Exiting.")
        sys.exit(1)

    return orthomcl_dirpath


def translate_focal_fasta_for_recret(species_name, focal_species_fasta, focal_species_aa_filepath):
    """
    Translates focal species FASTA file.
    
    :param species_name: focal species' name
    :param focal_species_fasta: filepath to the focal species' FASTA
    :param focal_species_aa_filepath: filepath for the outputed translated focal species' FASTA
    """
    # - Translate focal species FASTA CDS to AA and print it to file
    # Filename is the same basename of original FASTA file, plus suffix ".faa"; located in OrthoMCL dir
    if not os.path.exists(focal_species_aa_filepath) or os.path.getsize(focal_species_aa_filepath) == 0:
        logging.debug(f'Translating CDS file of [{species_name}] for diamond run...')
        cds_sequences, protein_sequences = None, None
        cds_sequences = read_fasta(focal_species_fasta) # TODO: read_fasta appears both in wgd_ksrates and ksrates... same function?
        protein_sequences = translate_cds(cds_sequences)
        with open(focal_species_aa_filepath, "w+") as f:
            for id in protein_sequences:
                f.write(f">{id}\n")
                f.write(f"{protein_sequences[id]}\n")
        logging.debug("")
    return

def merge_all_amino_fasta_for_recret(species_name, focal_species_aa_filepath,
                                   orthomcl_outdir_focal, original_37_fasta_dir_targz):
    """
    Merged translated focal species FASTA with the translated original FASTAs
    of the 37 angiosperms used for the ranking in Li et al. (2016).
    
    :param species_name: focal species' name
    :param focal_species_aa_filepath: filepath for the outputed translated focal species' FASTA
    :param orthomcl_outdir_focal: directory of OrthoMCL output for focal species
    :param original_37_fasta_dir_targz: path to tar containing the original 37 angiosperms' FASTAs (or a mock version for the test run)
    :return merged_fasta_outpath: file of merged 38 aminoacidic FASTAs 
    """
    # - Merge the translated FASTA of the focal species with the original FASTA files of the 37 angiosperms
    # The 37 angiosperms FASTAs are read from a tar file
    tar = tarfile.open(original_37_fasta_dir_targz)
    # Define the merged output filepath
    merged_fasta_outfile = f"{species_name}_original_37.faa"
    merged_fasta_outpath = os.path.join(orthomcl_outdir_focal, merged_fasta_outfile)
    # Merge to output file
    logging.debug(f"Merging FASTA files of original 37 angiosperms plus focal species into [{species_name}].")
    with open(merged_fasta_outpath, "wb") as merged_fasta_out:
        # First add focal species FASTA
        with open(focal_species_aa_filepath, "rb") as focal_species_aa:
            shutil.copyfileobj(focal_species_aa, merged_fasta_out)
        # Then add each of the 37 angiosperms FASTAs
        for original_fasta in tar.getmembers():
            content_original_fasta = tar.extractfile(original_fasta)
            shutil.copyfileobj(content_original_fasta, merged_fasta_out)
    logging.debug("")
    return merged_fasta_outpath


def make_homology_table_for_orthomcl(species_name, orthomcl_outdir_focal, 
                                     n_threads, merged_fasta_outpath):
    """
    Performs a diamond homology search on the merged FASTA to indentify homologous pairs.
    
    :param species_name: focal species' name
    :param orthomcl_outdir_focal: directory of OrthoMCL output for focal species
    :param n_threads: number of threads used to parallelize diamond
    """
    # - Make diamond database file
    database_filename = f"{species_name}_original_37.dmnd"
    database_filepath = os.path.join(orthomcl_outdir_focal, database_filename)    
    if os.path.exists(database_filepath):
        logging.info(f"diamond database file [{database_filename}] already exists.")
    else:
        logging.info(f"diamond database file [{database_filename}] not found, will generate it")
        run_diamond_makedb(merged_fasta_outpath, database_filepath, n_threads)
        logging.info("")

    # - Generate homology table output
    dmd_homology_table_filename = f"{species_name}_original_37.dmd.tsv"
    dmd_homology_table_filepath = os.path.join(orthomcl_outdir_focal, dmd_homology_table_filename)
    if not os.path.exists(dmd_homology_table_filepath):
        logging.info(f"diamond homology table [{dmd_homology_table_filename}] not found, will generate it")
        run_diamond_search(database_filepath, merged_fasta_outpath, dmd_homology_table_filepath, n_threads)
    else:
        logging.info(f"diamond homology table [{dmd_homology_table_filename}] already exists.")
    return dmd_homology_table_filepath


def generate_gfs_orthomcl(species_name, focal_species_aa_filepath, original_37_fasta_dir_targz, orthomcl_output_renamed_path, 
                          orthomcl_outdir_focal, n_threads, inflation, dmd_homology_table_filepath,
                          parsed_homology_table_filepath, use_original_orthomcl_version):
    """
    Generates the "all.gg" file needed by OrthoMCL, which contains all gene IDs per FASTA.
    Runs OrthoMCL: either mode 3 accepting as input the homology search table, or mode 4 
    accepting as input the already parsed version of the homology table (BPO file) generated
    in a previous OrthoMCL run, if user has provided it.
    
    :param species_name: focal species' name
    :param focal_species_aa_filepath: filepath for the outputed translated focal species' FASTA
    :param original_37_fasta_dir_targz: path to tar containing the original 37 angiosperms' FASTAs (or a mock version for the test run)
    :param orthomcl_output_renamed_path: path to the renamed OrthoMCL orthogroups output file (.out)
    :param orthomcl_outdir_focal: directory of OrthoMCL output for focal species
    :param n_threads: number of threads used to parallelize diamond
    :param inflation: OrthoMCL parameter for cluster granularity [default: 1.5]
    :param dmd_homology_table_filepath: path to homology search table
    :param parsed_homology_table_filepath: user-provided path to existing parse homology table (BPO file)
    """    
    # - Make all.gg file for OrthoMCL
    output_all_gg_file = os.path.join(orthomcl_outdir_focal, "all.gg")
    if os.path.exists(output_all_gg_file) and os.path.getsize(output_all_gg_file) != 0:
        logging.debug(f"File [{output_all_gg_file}] for OrthoMCL run already exists.")
    else:
        logging.debug(f"Generating [{output_all_gg_file}] for OrthoMCL run.")
        make_all_gg_file_for_orthomcl(species_name, focal_species_aa_filepath, original_37_fasta_dir_targz, output_all_gg_file)
    logging.debug("")

    # - Run OrthoMCL for orthogroup clustering
    logging.info("Starting OrthoMCL run... Will take some time!")
    orthomcl_dirpath = run_orthomcl_command(dmd_homology_table_filepath, parsed_homology_table_filepath, 
                                            output_all_gg_file, orthomcl_outdir_focal, inflation=inflation,
                                            use_original_orthomcl_version=use_original_orthomcl_version)
    # - Copy and rename output to final recret directory
    orthomcl_output_orig_path = os.path.join(orthomcl_dirpath, "all_orthomcl.out")
    shutil.copy(orthomcl_output_orig_path, orthomcl_output_renamed_path)
    logging.info("")
    return


def clean_orthomcl_output(output_orthomcl_path):
    """
    Convert the orthomcl output (.out) into a dataframe showing orthogroup ID,
    gene count, number of taxa and gene list.
    
    :param output_orthomcl_path: path to raw OrthoMCL output
    """
    output_orthomcl_dirty = read_csv(output_orthomcl_path, sep="\t", index_col=0, header=None)

    names_genecount_mix_taxa = DataFrame.from_records(output_orthomcl_dirty.index.str.split("(")) #      ORTHOMCL4  137 genes,1 taxa):
    names = names_genecount_mix_taxa.iloc[:,0] # ORTHOMCL1

    genecount_taxa_dirty = DataFrame.from_records(names_genecount_mix_taxa.iloc[:,1].str.split(",")) # 2 genes  1 taxa):
    genecount = DataFrame.from_records(genecount_taxa_dirty.iloc[:,0].str.split(" "))[0].astype(int) # 2

    taxa_with_string = DataFrame.from_records(genecount_taxa_dirty.iloc[:,1].str.split(")"))[0] # 1 taxa
    taxa = DataFrame.from_records(taxa_with_string.str.split(" "))[0].astype(int) # 1

    genelist_dirty = output_orthomcl_dirty.copy() # ORTHOMCL0(46 genes,2 taxa):    PDK_30s1005061L001(pdac) PDK_30s1009311L001(p
    genelist_dirty.index = range(0, len(genelist_dirty)) # 0     PDK_30s1005061L001(pdac) PDK_30s1009311L001(p.
    genelist = genelist_dirty.iloc[:,0].str.strip().str.split(" ") # 0 [PDK_30s1005061L001(pdac), PDK_30s1009311L001(pdac) ...

    orthomcl_output = concat([names, genecount, taxa, genelist], axis=1)
    orthomcl_output.columns = ["Orthogroup", "Gene count", "Number of taxa", "Gene list"]
    
    return orthomcl_output


def count_original_and_focal_genes_in_orthomcl_gfs(orthomcl_output_renamed_path):
    """
    Counts how many original 37 species's genes there are and how many
    focal species' genes there are in every new OrthoMCL GF.
    
    :param orthomcl_output_renamed_path: OrthoMCL orthogroup output, has been renamed and moved
    :return: OrthoMCL output updated with number of original and new genes per orthogroup
    """
    # - Convert OrthoMCL orthogroup table into a clean pandas DataFrame
    #   Example of cleaned-up layout:
    #           Orthogroup   Gene count   Number of taxa                Gene list
    #     0      ORTHOMCL0           46                2   [Gene0001(species), ... ]
    clean_orthomcl_gfs = clean_orthomcl_output(orthomcl_output_renamed_path)

    # - Get list and set of genes per orthogroup
    gene_list_with_species_names = clean_orthomcl_gfs["Gene list"].copy()
    clean_orthomcl_gfs["Gene set"] = clean_orthomcl_gfs["Gene list"].apply(lambda genes: {gene.split("(")[0] for gene in genes})

    # - Get number of original (old) genes belonging to the 37 angiosperms (based on whether they are in the list of original species names)
    # NOTE: focal species name is required by configuration to be different from any of the original 37species: e.g. "atha37" is forbidden
    # Replace genes belonging to focal species with a "nan"
    original_genes_and_nan_focal = gene_list_with_species_names.apply(lambda genes: {gene.split("(")[0] if gene.split("(")[1].split(")")[0] in original_37_names else nan for gene in genes})
    # Remove the focal species' nan values. How to find them? "nan" is the only value for which holds that "value != value"
    clean_orthomcl_gfs["IDs original genes"] = original_genes_and_nan_focal.apply(lambda genes: {gene for gene in genes if gene == gene})
    clean_orthomcl_gfs["Num original genes"] = clean_orthomcl_gfs["IDs original genes"].str.len()

    # - Get the number of genes of the FOCAL (or NEW) species
    # Replace genes from the original 37 species with a nan
    focal_genes_and_nan_original = gene_list_with_species_names.apply(lambda genes: {gene.split("(")[0] if gene.split("(")[1].split(")")[0] not in original_37_names else nan for gene in genes})
    # Remove the original species' nan values. How to find them? "nan" is the only value for which holds that "value != value"
    clean_orthomcl_gfs["IDs focal species genes"] = focal_genes_and_nan_original.apply(lambda genes: {gene for gene in genes if gene == gene})
    clean_orthomcl_gfs["Num focal species genes"] = clean_orthomcl_gfs["IDs focal species genes"].str.len()

    return clean_orthomcl_gfs


def match_and_reconstruct_gfs(top, gf_ids_ranked_to_top, original_and_focal_genes_in_orthomcl_gfs, original_top_gf_genes_raw,
                        matched_orthomcl_gfs_path, max_extra_original_genes_in_new_gfs, min_common_old_genes_in_clade_gfs):
    """
    Matches the original GFs to the new GFs and select the good matches for GF reconstruction.
    
    This is done first by macthing an original top GF to any newly generated OrthoMCL GFs
    sharing at least one of its original genes (original = coming from the 37 angiosperms' genes);
    then, only the matched new GFs containing few non-shared (extra) original genes are accepted,
    in order to not inflate the reconstructed GF with enexpected original genes.
    
    Note that sharing genes is possible because the new OrthoMCL GFs are generated with
    the same genome versions of the the 37 angiosperms used to make the lambda ranking (Li et al., 2016).

    Detailed steps and output column layout:
    1) Match original GF to new GF(s) based on gene sharing:
    - ID and number of genes of original GF. Col_1 [Original GF ID] and Col_2 [Original GF genes].
    - ID and number of genes of each matched new GF. Col_3 [Matched GF IDs] and Col_4 [Matched GF genes].
    - Number of focal species' genes in each matched new GF. Col_7 [Focal species genes].

    2) Compare original GF with each matched new GF:
    - Number of original genes shared between the original GF and each matched new GF. Col_5 [Shared original genes].
    - Number of original genes shared between the original GF and and all matched GFs. Col_8 [Merged shared original genes].
    - Number of original genes not shared (extra) between the original GF and each matched GF. Col_6 [Extra original genes].
        
    3) Select and merge the well-matched new GF(s) into a reconstructed GF:
    - Selects each matched new GF based on max_extra_original_genes_in_new_gfs. Col_9 [Accept/Reject].
    - Merges the accepted matches into the reconstructed GF.
    - Number of original genes shared between original GF and reconstricted GF. Col_10 [Merged accepted shared original genes].
    
    4) Compare original GF with reconstructed GF:
    - Number of original genes not shared (extra) between original GF and reconstructed GF. Col_11 [Merged accepted extra original genes].
    - # TODO: NOT USED, REMOVE? Percentage of original genes not shared (extra), wrt the number of shared original genes. Col_12 [Perc merged accepted extra original genes].
    - IDs of non-shared (extra) original genes in reconstructed GF. Col_15 [IDs extra original genes].
    
    - Number of missing original genes in the reconstructed GF. Col_13 [Missing original genes].
    - Percentage of missing original genes in the reconstructed GF, wrt number of genes in the original GF. Col_14 [Perc missing original genes].
    - IDs of missing original genes in reconstructed GF. Col_16 [IDs missing original genes].
    
    5) Evaluate GF reconstruction quality:
    - Assigns quality score to reconstructed GF: good if there are only few missing original genes. Col_17 [Coverage score].

    :param top: number of original top GFs requested, e.g. 2000
    :param gf_ids_ranked_to_top: IDs of the original top GFs
    :param original_and_focal_genes_in_orthomcl_gfs: OrthoMCL table plus analysis of 37 species' and focal species' genes. 
    :param original_top_gfs_raw_path: DataFrame with member genes of all 9178 core-angiosperm GFs used in from Tasdighian (2017)
    :param matched_orthomcl_gfs_path: path to table containing the matching between original top GFs and new OrthoMCL GFs
    :param max_extra_original_genes_in_new_gfs: max number of original genes not shared between the original 
    top GF and the matched GFs for this match to be accepted
    :param min_common_old_genes_in_clade_gfs: NOTE [not used] minimum required number of original genes shared
    :return: table listing the original GFs sharing genes with the new OrthoMCL GFs, plus...
    """
    # Read and clean-up orthogroup gene members of the original 9178 GFs
    original_top_gf_genes_raw.columns = ["Gene family ID", "Number of species", "Genes"]
    original_top_gfs_gene_column = original_top_gf_genes_raw["Genes"].str.split(",")
    original_top_gfs_gene_set = original_top_gfs_gene_column.apply(lambda genes: {gene.split('|')[1] for gene in genes})
    original_top_gfs = concat([original_top_gf_genes_raw["Gene family ID"], original_top_gf_genes_raw["Number of species"], original_top_gfs_gene_set], axis=1)
    # Example of cleaned-up layout:
    #   Gene family ID    Number of species    Genes
    #   ORTHO002222                      36    {Cotton_D_gene_10033278, ... }
    
    # Select only the top X GFs
    original_top_gfs_genes = original_top_gfs.loc[original_top_gfs["Gene family ID"].isin(gf_ids_ranked_to_top)]

    # Initialize DataFrame of reconstructed gene families
    top_gf_matched_to_orthomcl_gfs = DataFrame()
    
    # Loop through the original top X GFs
    for gf_index, gf in gf_ids_ranked_to_top.iteritems():
        
        # 1) Match original GF to any new GF sharing at least one original gene
        # ---------------------------------------------------------------------

        # Get the genes of original GF
        original_genes = original_top_gfs_genes.loc[original_top_gfs_genes["Gene family ID"] == gf]["Genes"].iloc[0]
        
        # Get the ID of the new GF(s) that share some genes with the orginal GF; they are the "matched new GFs"
        gf_intersection_dirty = original_and_focal_genes_in_orthomcl_gfs["Gene set"].apply(lambda clade_genes: clade_genes.intersection(original_genes))
        gf_intersection = gf_intersection_dirty.loc[gf_intersection_dirty.str.len() != 0]
        gf_intersection_sorted = gf_intersection.sort_index()
        # If there are no new GF sharing any genes, skip current original GF
        if gf_intersection_sorted.empty:
            continue
        
        # Get data rows of the matched new GF(s)
        clade_orthomcl_with_intersection = original_and_focal_genes_in_orthomcl_gfs.iloc[list((gf_intersection_sorted.index))].sort_index()
        
        # Column 1: ID of original GF
        # Note that it contains same ID repeated by the number of matched new GFs
        original_gf_id_repeated = Series([gf]*len(clade_orthomcl_with_intersection), index=list(clade_orthomcl_with_intersection.index))
        original_gf_id_repeated.name = "Original GF ID"
        # Column 2: Number of genes in original GF
        # Note that it contains same number repeated by the number of matched new GFs
        original_gf_num_genes = Series([len(original_genes)]*len(clade_orthomcl_with_intersection), index=list(clade_orthomcl_with_intersection.index))
        original_gf_num_genes.name = "Original GF genes"
        original_gf_num_genes_int = original_gf_num_genes.iloc[0]
        
        # Column 3: ID of each matched new GF
        clade_gf_ids = clade_orthomcl_with_intersection["Orthogroup"]
        clade_gf_ids.name = "Matched GF IDs"
        # Column 4: Number of genes in each matched new GF
        num_genes_clade_gfs = clade_orthomcl_with_intersection["Gene count"]
        num_genes_clade_gfs.name = "Matched GF genes"

        # Column 7: Number of focal species' genes in each matched new GF
        number_new_genes = clade_orthomcl_with_intersection["Num focal species genes"]
        number_new_genes.name = "Focal species genes"


        # 2) Compare original GF and each matched new GF
        # ----------------------------------------------

        # Column 5: Number of shared original genes between the original GF and each matched new GF
        num_common_genes = gf_intersection_sorted.str.len()
        num_common_genes.name = "Shared original genes"
        # List of all shared genes
        list_common_genes = gf_intersection_sorted.copy()
        # Merge the shared genes from all the matched GFs.
        merged_common_old_genes_ids = set()
        # I use Python sets() to be sure I don't retain duplicates, but actually OrthoMCL doesn't assign a gene to more than one GF...
        for row in list_common_genes:
            if merged_common_old_genes_ids.intersection(row) != set():
                logging.error("Gene redundancy during top GF reconstruction. Shouldn't happen! Needs debugging.")
                logging.error("Exiting.")
                sys.exit(1)
            merged_common_old_genes_ids = merged_common_old_genes_ids.union(row)
        # Column 8: Number of shared original genes between the original GF and all matched GFs
        # Same number repeated by the number of matched GFs
        merged_common_genes_int = len(merged_common_old_genes_ids)
        merged_common_genes = Series([merged_common_genes_int]*len(clade_orthomcl_with_intersection), index=list(clade_orthomcl_with_intersection.index))
        merged_common_genes.name = "Merged shared original genes"
                
        # Column 6: Number of original genes in macthed new GFs that are not appearing in the original GF (therefore, "extra original genes")
        # Number original genes in new GFs (i.e. no matter whether they are shared or not with the original GF)
        num_old_genes_clade_gfs = clade_orthomcl_with_intersection["Num original genes"]
        # Number of extra original genes in matched new GFs
        num_extra_old_genes_clade_gfs = num_old_genes_clade_gfs - num_common_genes
        num_extra_old_genes_clade_gfs.name = "Extra original genes"


        # 3) Reconstruct the GF using the well-matched new GF(s)
        # ------------------------------------------------------

        # Column 8: Accept (True) or reject (False) matched new GFs based on number of extra original genes
        # The "True" matches have only few extra original genes (<= max_extra_old)
        accept_reject_matched_clade_gf = (num_extra_old_genes_clade_gfs <= max_extra_original_genes_in_new_gfs)
        accept_reject_matched_clade_gf.name = "Accept/Reject"
        # TODO: perhaps also check for "at least one shared gene" or "at least some shared genes (>= min_common_old)"
        # accept_reject_matched_clade_gf = (num_extra_old_genes_clade_gfs <= max_extra_original_genes_in_new_gfs) & (num_common_genes >= 1)


        # TODO: CHECK THIS OUT. Automatically set a True match to False because of default criteria
        # | ((num_extra_old_genes_clade_gfs == 14) & (num_common_genes >= conf.min_common_old_genes_in_clade_gfs))
        # Manually set to false when there are 14 extra original genes and only 1 shared original gene
        # Cause 1 is too little to justify the integration of 14 extra genes...
        set_to_false_nonetheless = (num_extra_old_genes_clade_gfs == 14) & (num_common_genes <= min_common_old_genes_in_clade_gfs)
        set_to_false_nonetheless_index = set_to_false_nonetheless.loc[set_to_false_nonetheless == True].index
        if len(set_to_false_nonetheless_index) != 0:
            logging.debug(f"Automatically setting a True match to False: marker GF {gf} matches OrthoMCL GF {clade_orthomcl_with_intersection.at[set_to_false_nonetheless_index[0], 'Orthogroup']} but has actually too many extra genes compared to shared genes.")
            accept_reject_matched_clade_gf.loc[set_to_false_nonetheless_index] = False
        # else:
        #     logging.debug("No need to automatically set a macthed GF to False.")
        
        # TODO: These additional checks are not run, and must be thoroughly re-evaluated before being accepted
        # run_additional_checks(gf, gf_index, accept_reject_matched_clade_gf, num_common_genes, num_extra_old_genes_clade_gfs,
        #                   original_gf_num_genes_int, clade_orthomcl_with_intersection,
        #                   gf_to_set_false_by_asking_user=False, retain_even_if_many_extra_old_genes=False)    
        
        # Column 9: Merged count of shared original genes that belong to the accepted GFs ("Accept/Reject" = True)
        num_common_genes_accepted = num_common_genes.loc[accept_reject_matched_clade_gf == True]
        # Merge all accepted elements (we already know there are no duplicates!)
        merged_common_genes_retained_int = num_common_genes_accepted.sum()
        # Same number repeated by the number of matched GFs
        merged_common_genes_retained = Series([merged_common_genes_retained_int]*len(clade_orthomcl_with_intersection), index=list(clade_orthomcl_with_intersection.index))
        merged_common_genes_retained.name = "Merged accepted shared original genes"


        # 4) Compare original GF with reconstructed GF
        # --------------------------------------------

        # Column 10: Merged accepted extra original genes
        # Should be as least as possible
        num_extra_old_genes_retained = num_extra_old_genes_clade_gfs.loc[accept_reject_matched_clade_gf == True]
        merged_extra_old_genes_retained_int = num_extra_old_genes_retained.sum()
        # Same number repeated by the number of matched GFs
        merged_extra_old_genes_retained = Series([merged_extra_old_genes_retained_int]*len(clade_orthomcl_with_intersection), 
                                                index=list(clade_orthomcl_with_intersection.index))
        merged_extra_old_genes_retained.name = "Merged accepted extra original genes"
        merged_extra_old_genes_retained_int = merged_extra_old_genes_retained.iloc[0]
        
        # Column 11: Percentage of merged accepted extra original genes
        # Calculated as percentage of extra original genes out of the complete set of genes in the matched new GFs.
        # Should be as least as possible
        perc_merged_extra_old_genes_retained = round((merged_extra_old_genes_retained/merged_common_genes_retained)*100, 2) # TODO: MEANING? WHY USEFUL?
        perc_merged_extra_old_genes_retained.name = "Perc merged accepted extra original genes"
        perc_merged_extra_old_genes_retained_float = perc_merged_extra_old_genes_retained.iloc[0]

        # Column 12: Difference between expected original genes and collected original genes after filtering for the True matches
        # Column X: Percentage of missing original genes out of total number of original genes in original GF
        if original_gf_num_genes_int >= merged_common_genes_retained_int:
            missing_old_genes = original_gf_num_genes - merged_common_genes_retained_int
            missing_old_genes_int = missing_old_genes.iloc[0]
            perc_missing_old_genes = round((missing_old_genes/original_gf_num_genes)*100, 2)
            perc_missing_old_genes_float = perc_missing_old_genes.iloc[0]
        else:
            missing_old_genes = merged_common_genes_retained_int - original_gf_num_genes
            logging.debug(f"Clade GFs has more original genes than the original GF: {missing_old_genes}")
        missing_old_genes.name = "Missing original genes"
        perc_missing_old_genes.name = "Perc missing original genes"
        
        # Get the IDs of the missing original gene
        missing_old_genes_ids = original_genes - gf_intersection
        missing_from_any_matched_gf = original_genes - merged_common_old_genes_ids
        # Make a Series to be added as a column
        ids_missing_old_genes_in_clade_orthomcl_retained_matches = missing_from_any_matched_gf
        # Same ID list repeated by the number of matched new GFs
        ids_missing_old_genes_in_clade_orthomcl_retained_matches = Series([ids_missing_old_genes_in_clade_orthomcl_retained_matches]*len(clade_orthomcl_with_intersection), index=list(clade_orthomcl_with_intersection.index))
        ids_missing_old_genes_in_clade_orthomcl_retained_matches.name = "IDs missing original genes"
        
        # Column 13: List of IDs of extra original genes (only belonging to the accepted macthed new GFs)
        # Get the index of the accepted (True) matches
        retained_gfs_index = num_extra_old_genes_retained.index
        # Get the list of all genes in the new matched GFs
        ids_tot_genes_in_clade_orthomcl_retained_matches = original_and_focal_genes_in_orthomcl_gfs["Gene set"].loc[retained_gfs_index]
        # Get the list of shared original genes in the new matched GFs
        ids_old_common_genes_in_clade_orthomcl_retained_matches = list_common_genes.loc[retained_gfs_index]
        # Get the list of genes belonging to the new species in the new matched GFs
        ids_new_genes_in_clade_orthomcl_retained_matches = original_and_focal_genes_in_orthomcl_gfs["IDs focal species genes"].loc[retained_gfs_index]
        # Get the extra original genes, by removing from the total genes the shared original genes and the genes coming from the new species
        ids_extra_old_genes_in_clade_orthomcl_retained_matches = ids_tot_genes_in_clade_orthomcl_retained_matches - ids_old_common_genes_in_clade_orthomcl_retained_matches - ids_new_genes_in_clade_orthomcl_retained_matches
        ids_extra_old_genes_in_clade_orthomcl_retained_matches.name = "IDs extra original genes"
        
        
        # 5) Evaluate GF reconstruction quality
        # -------------------------------------

        # Column 14: Assign score based on percentage of missing original genes and percentage of extra original genes
        # TODO: CHECK THIS OUT. Based on the percentage of extra genes (the smaller the better)
        if perc_merged_extra_old_genes_retained_float <= 5:
            automatic_score_str = "Good"
        else:
            automatic_score_str = "Bad"
        # Based on the percentage of missing genes (the smaller the better); it overwrites the previous one!
        if perc_missing_old_genes_float < 0:
            automatic_score_str = "Bad" # TODO: check if it ever happened and decide how to handle it
        elif perc_missing_old_genes_float <= 5:
            automatic_score_str = "Good"
        elif perc_missing_old_genes_float > 5 and perc_missing_old_genes_float <= 10:
            automatic_score_str = "Fair"
        elif perc_missing_old_genes_float > 10:
            automatic_score_str = "Bad"
        automatic_score = Series([automatic_score_str]*len(clade_orthomcl_with_intersection), index=list(clade_orthomcl_with_intersection.index))
        automatic_score.name = "Coverage score"

        # Print statement
        number_accepted_and_rejected_matched_new_gfs = len(gf_intersection_sorted)
        number_accepted_matched_new_gfs = len(accept_reject_matched_clade_gf.loc[accept_reject_matched_clade_gf == True])
        logging.debug(f"Original GF {gf} (ranked {gf_index}) was matched to {len(gf_intersection_sorted)} new GF(s) sharing some original genes")
        if number_accepted_matched_new_gfs == number_accepted_and_rejected_matched_new_gfs:
            logging.debug(f"All matched GFs were used to reconstruct the original GF:")
        elif number_accepted_matched_new_gfs != 0:
            logging.debug(f"However, only {number_accepted_matched_new_gfs} matches could be accepted for GF reconstruction:")
        else:
            logging.debug(f"However, no matches could be accepted for GF reconstruction:")
        logging.debug(f"- Original genes included in the original GF:            {original_gf_num_genes_int}")
        logging.debug(f"- Original genes retrieved in the reconstructed GF:      {merged_common_genes_retained_int}")
        logging.debug(f"- Missing original genes in the reconstructed GF:        {missing_old_genes_int} ({perc_missing_old_genes_float}%)")
        logging.debug(f"- Extra original genes acquired in the reconstructed GF: {merged_extra_old_genes_retained_int}")
        logging.debug("")
        
        # All columns together:
        comparison_overlapping_gfs = concat([original_gf_id_repeated, original_gf_num_genes, clade_gf_ids, num_genes_clade_gfs,
                                                num_common_genes, num_extra_old_genes_clade_gfs, number_new_genes,
                                                merged_common_genes, accept_reject_matched_clade_gf,
                                                merged_common_genes_retained,
                                                merged_extra_old_genes_retained, perc_merged_extra_old_genes_retained,
                                                missing_old_genes, perc_missing_old_genes,
                                                ids_extra_old_genes_in_clade_orthomcl_retained_matches,
                                                ids_missing_old_genes_in_clade_orthomcl_retained_matches,
                                                automatic_score], axis=1) # list_common_genes as last col
        
        top_gf_matched_to_orthomcl_gfs = concat([top_gf_matched_to_orthomcl_gfs, comparison_overlapping_gfs])
                    
    # Print to file
    top_gf_matched_to_orthomcl_gfs.to_csv(matched_orthomcl_gfs_path, sep="\t")
    return


def run_additional_checks(gf, gf_index, accept_reject_matched_clade_gf, num_common_genes, num_extra_old_genes_clade_gfs,
                          original_gf_num_genes_int, clade_orthomcl_with_intersection,
                          gf_to_set_false_by_asking_user=False, retain_even_if_many_extra_old_genes=False):
    """
    """
    # TODO: CHECK THIS OUT. Manually set a True match to False because of user's criteria (by asking them with the input() function)
    if gf_to_set_false_by_asking_user == True:
        logging.warning(f"Do you want to set to False some other matched GFs, based on non-default criteria?")
        logging.warning(f"E.g. Lamiales: marker ORTHO000675 matches the small GF ORTHOMCL20267, which has only")
        logging.warning(f"1 shared original gene, but also 5 extra original genes coming from the same species (Ccaj).")
        while gf_to_set_false_by_asking_user != "":
            gf_to_set_false_by_asking_user = input("Write here the matched GF NUMBER (e.g. 20267, not the full ID), or just press enter to skip: ")
            if gf_to_set_false_by_asking_user != "":
                gf_to_set_false_by_asking_user = int(gf_to_set_false_by_asking_user)
                accept_reject_matched_clade_gf.loc[gf_to_set_false_by_asking_user] = False
    # else:
    #     logging.debug("No need to manually set a macthed GF to false.")

    # if conf.clade == "lamiales" and gf == "ORTHO000675":
    #     print(f"The following matches were set from True to False by user:")
    #     matched_gf_index = 20267
    #     print(f"Clade: {conf.clade}, marker GF: {gf}, match GF index: {matched_gf_index}")
    #     accept_reject_matched_clade_gf.loc[matched_gf_index] = False
        
    # TODO: CHECK THIS OUT. Perhaps it is instead needed to set a False GF back to True: some extra original, but many good original genes.
    if retain_even_if_many_extra_old_genes:
        matches_currently_set_to_false = num_common_genes.loc[accept_reject_matched_clade_gf == False]
        set_to_true_nonetheless = (matches_currently_set_to_false >= 10)
        set_to_true_nonetheless_index_list = set_to_true_nonetheless.loc[set_to_true_nonetheless == True].index
        if len(set_to_true_nonetheless_index_list) != 0:
            logging.warning(f'Marker GF {gf} matches the following OrthoMCL GF(s), which have many "extra" original genes but also many good original genes (>= 10).')
            for set_to_true_nonetheless_index in set_to_true_nonetheless_index_list:
                logging.warning(f'{clade_orthomcl_with_intersection.loc[set_to_true_nonetheless_index, "Orthogroup"]} (ranked {gf_index})')
                num_common = num_common_genes.loc[set_to_true_nonetheless_index]
                num_extra = num_extra_old_genes_clade_gfs.loc[set_to_true_nonetheless_index]
                logging.warning(f'Despite having {num_extra} extra original genes, it also shares {num_common}/{original_gf_num_genes_int} original genes.')
                logging.warning(f"This match is currently set to False; do you want to convert it to True and include it in the marker GF reconstruction?")
                logging.warning(f"If yes, remember to check whether you can get rid of the extra original genes based on the matched GF gene tree.")
                logging.warning(f"Otherwise, if you eventually want to exclude this match, rerun the code and rejectthis option.")
                set_to_true_nonetheless_choice = input(f"Do you want to include the match(es) listed above for this GF reconstruction? [y/N]: ").lower()
                if set_to_true_nonetheless_choice == "y":
                    logging.info(f"Manually setting the above False match(es) to True.")
                    accept_reject_matched_clade_gf.loc[set_to_true_nonetheless_index] = True
        # else:
        #     logging.debug("No need to manually set a macthed GF to true.")
    
    # OLD VERSION ONLY COUNTING ON max_extra_old: accept_reject_matched_clade_gf = num_extra_old_genes_clade_gfs <= conf.max_extra_original_genes_in_new_gfs


def log_match_quality(top, matched_orthomcl_gfs):
    """
    Reads the number of original top GFs sharing some genes with the new OrthoMCL GFs
    (i.e. have at least one matched new GF).
    --> The number is the number of unique original GF IDs in matched_orthomcl_gfs.
        E.g. 1500/2000 original top GFs share some genes
    These are the original GFs with the potential of being reconstructed.
    
    Reads also the quality of the gene overlap (coverage) between the  
    filtered original GFs and their matched new GFs.
    Quality is based on percentage of missing genes (the more missing, the worse the quality).
    --> The quality is read from the "Coverage score" column in matched_orthomcl_gfs.
        E.g.
        - "Good" gene coverage : 1000/1500
        - "Fair" gene coverahe : 300/1500
        - "Bad" gene coverage  : 200/1500
    Note: all matched new GFs are taken into account, namely both the accepted matches
    (Accept/Reject = True) and rejected matches (False); they were rejected
    when containing too many extra original genes.
    
    :param top: top requested GFs from the reciprocal retention ranking
    :param matched_orthomcl_gfs: table with matches between original top GF and new OrthoMCL GFs
    """
    # Get reconstruction table with an "absolute index" (i.e. [0, 1, ..., len(table)] instead of the OrthoMCL GF indeces (e.g. [566, 1407, ...])
    matched_orthomcl_gfs_absolute_index = matched_orthomcl_gfs.copy()
    matched_orthomcl_gfs_absolute_index.index = range(len(matched_orthomcl_gfs))
    # Get the non-redundant list of top X GF IDs from the full table (length is equal to "top X", e.g. 2000)
    all_covered_top_gfs_ids = matched_orthomcl_gfs_absolute_index["Original GF ID"].drop_duplicates()
    num_covered_top_gfs = len(all_covered_top_gfs_ids)

    # State how many original GFs are sharing some genes with the new OrthoMCL GFs
    # This means that they have the potential to be reconstructed!
    if num_covered_top_gfs == top:
        logging.info(f"All top {top} GFs share some genes with the new OrthoMCL GFs")
    else:
        logging.warning(f"Only {num_covered_top_gfs}/{top} top GFs share some genes with the new OrthoMCL GFs (i.e. they have matches); all the remaining top GFs will be discarded")
        logging.warning(f"Note: a lower number than the requested {top} can be expected to a certain extent;")
        logging.warning(f"however keep in mind that too low numbers can generate poor Ks distributions!")
    # - Use the indexes of the non-redundant GF ID list to filter a non-redundant set of lines from the full table
    matched_orthomcl_gfs_absolute_index_top = matched_orthomcl_gfs_absolute_index.iloc[all_covered_top_gfs_ids.index]
    # - Get the score from this non-redundant set of lines
    automatic_score_series = matched_orthomcl_gfs_absolute_index_top.loc[matched_orthomcl_gfs_absolute_index_top.index, "Coverage score"]
    # - Count how many Good, Fair and Bad
    automatic_score_counts = automatic_score_series.value_counts()
    try:
        num_good = automatic_score_counts["Good"]
    except KeyError:
        num_good = 0
    ratio_good = round(num_good/num_covered_top_gfs*100, 2)
    try:
        num_fair = automatic_score_counts["Fairly"]
    except KeyError:
        num_fair = 0
    ratio_fair = round(num_fair/num_covered_top_gfs*100, 2)
    try:
        num_bad = automatic_score_counts["Bad"]
    except KeyError:
        num_bad = 0
    ratio_bad = round(num_bad/num_covered_top_gfs*100, 2)
    logging.info(f"Evaluating gene coverage of the {num_covered_top_gfs} original GFs, based on missing original genes:")
    logging.info(f'- Good gene coverage : {num_good}/{num_covered_top_gfs} ({ratio_good}%)')
    logging.info(f'- Fair gene coverage : {num_fair}/{num_covered_top_gfs} ({ratio_fair}%)')
    logging.info(f'- Bad gene coverage  : {num_bad}/{num_covered_top_gfs} ({ratio_bad}%)')
    
    return


def reconstruct_top_gfs(species_name, output_mcl_path, output_mcl_path_with_gf_ids, 
                        original_and_focal_genes_in_orthomcl_gfs, matched_orthomcl_gfs):
    """
    Reconstructs the original top GFs by merging their accepted matched new GFs.
    Then it kicks out the 37 angiosperms' genes in each reconstructed GF 
    to finally get the focal species' reconstructed top GFs.
    Outputs as a MCL-like gene family file, with one GF per row.
    
    Note: the number of reconstructed GFs can be lower then the requested one,
    e.g. 1460/2000, because of some filtering in previous steps:
    1) Not all original top 2000 GFs might match the new OrthoMCL GFs (i.e. shared some genes)
    2) Some of the GFs that did share some genes were however introducing too many 
       extra original genes, unexpectadly inflating the counts (Accept/Reject = False)

    :param species_name: focal species' informal name
    :param output_mcl_path: path to output file with reconstructed GF (MCL-like format)
    :param output_mcl_path_with_gf_ids: path to output file with reconstructed GF (MCL-like format with original top GF ID)
    :param original_and_focal_genes_in_orthomcl_gfs: OrthoMCL table plus analysis of 37 species' and focal species' genes. 
    :param matched_orthomcl_gfs: table with matches between original top GF and new OrthoMCL GFs
    :return: MCL-like table of reconstructed top X GFs for focal species
    """
    # Merge genes form all accepted matched new GFs to generate the reconstructed GF
    # ----------------------------------------------------------------------------------
    logging.info("Merge genes form all accepted matched new GFs to generate the reconstructed GFs")
    
    # - Map each gene ID to its species ID
    # Make a clean gene ID list by removing species id in brackets: e.g. from [gene1(species1), ...] to [gene1, ...]
    original_and_focal_genes_in_orthomcl_gfs["Gene IDs"] = original_and_focal_genes_in_orthomcl_gfs["Gene list"].apply(lambda genes: [gene.split("(")[0] for gene in genes])
    # Make a clean species ID list by removing the gene id: e.g. from [gene1(species1), ...] to [species1, ...]
    original_and_focal_genes_in_orthomcl_gfs["Species IDs"] = original_and_focal_genes_in_orthomcl_gfs["Gene list"].apply(lambda genes: [gene.split("(")[1].split(")")[0] for gene in genes])

    map_gene_id_to_species_id = {}
    for gene_id_list, species_id_list in zip(original_and_focal_genes_in_orthomcl_gfs["Gene IDs"], original_and_focal_genes_in_orthomcl_gfs["Species IDs"]):
        for gene_id, species_id in zip(gene_id_list, species_id_list):
            map_gene_id_to_species_id[gene_id] = species_id

    # - Merge together all genes from the accepted matched new GFs (Accept/Reject = True) to obtain the reconstructed GF
    num_matched_orthomcl_gfs = len(matched_orthomcl_gfs)
    matched_orthomcl_gfs_accepted = matched_orthomcl_gfs.loc[matched_orthomcl_gfs["Accept/Reject"] == True]
    num_matched_orthomcl_gfs_accepted = len(matched_orthomcl_gfs_accepted)
    # Log if some original GFs had exclusively poor matches with the new OrthoMCL GFs 
    # (i.e. too many extra original genes!) and therefore will be ignored 
    if num_matched_orthomcl_gfs_accepted != num_matched_orthomcl_gfs:
        logging.warning(f"Only {num_matched_orthomcl_gfs_accepted}/{num_matched_orthomcl_gfs} top GFs have at least one accepted match within the new OrthoMCL GFs and can thus be reconstructed; the remaining ones will be discarded")
        logging.warning(f"Note: a lower number than {num_matched_orthomcl_gfs} can be expected to a certain extent;")
        logging.warning(f"however, keep in mind that too low numbers may generate poor Ks distributions!")
    logging.info("")
    
    reconstructed_top_gfs = {} # key: original GF, value: gene IDs from all accepted matched new GFs
    # Loop through top X original GFs that had at least one accepted matched new GFs (Accept/Reject = True)
    for gf in matched_orthomcl_gfs_accepted["Original GF ID"]:
        reconstructed_top_gfs[gf] = set()
        
        # Get matched GF IDs to the original GF
        all_matches = matched_orthomcl_gfs_accepted.loc[matched_orthomcl_gfs_accepted["Original GF ID"] == gf]
        
        # Loop through all matched new GFs of the current original GF
        for index, row in all_matches.iterrows():
            match = row["Matched GF IDs"]
            
            # if match["Coverage score"] == "Bad":
            #     continue
            
            # Get the genes of accepted GFs from OrthoMCL output file
            match_orthomcl_out = original_and_focal_genes_in_orthomcl_gfs.loc[original_and_focal_genes_in_orthomcl_gfs["Orthogroup"] == match]
            match_orthomcl_out_genes = match_orthomcl_out.loc[index, "Gene set"]
            
            # # For some GFs there has been a manual removal of certain gene ourgoups because of belonging to non-core GFs)
            # if gf in markers_to_be_cleaned and match in markers_to_be_cleaned[gf]:
            #     match_orthomcl_out_genes_checked = match_orthomcl_out_genes - markers_to_be_cleaned[gf][match]
            # else: # These other markers didn't need any manual cleanup
            match_orthomcl_out_genes_checked = match_orthomcl_out_genes
            
            # Keep the merged genes in memory
            reconstructed_top_gfs[gf] = reconstructed_top_gfs[gf].union(match_orthomcl_out_genes_checked)
                
    reconstructed_top_gfs_genes = DataFrame.from_dict(reconstructed_top_gfs, orient="index")
    reconstructed_top_gfs_df_species = reconstructed_top_gfs_genes.applymap(lambda x: map_gene_id_to_species_id[x], na_action="ignore")
    
    # Kick out the original 37 angiosperms' genes to obtain the reconstructed rec.ret. GF for the focal species
    # ---------------------------------------------------------------------------------------------------------

    focal_rec_ret_gfs = {}
    # Keep count of the GFs for which the focal species has:
    # - no genes at all (will be ignored from now on)
    # - only one gene (singleton GFs, also will be ignored because Ks estimation needs at least two paralog genes)
    # - two or more genes (these GFs will be used to build the Ks plot)
    singletons_gfs_to_be_ignored = 0
    gfs_with_no_gene_of_focal = 0
    gfs_with_two_or_more_genes_focal = 0

    # Removing the original 37 angiosperms' genes
    for index_reconsturcted_gf, reconstructed_top_gf in reconstructed_top_gfs_df_species.iterrows():
        species_id_dropna = reconstructed_top_gf.dropna()
        # If the focal species is not in the reconstructed GF, we'll ignore this GF
        if species_name not in species_id_dropna.values:
            gfs_with_no_gene_of_focal += 1
        # If the focal species appears at least once...
        else:
            # Retain only the focal species' genes
            focal_genes_index = species_id_dropna.loc[species_id_dropna == species_name].index
            focal_genes_ids = reconstructed_top_gfs_genes.loc[index_reconsturcted_gf, focal_genes_index]
            focal_rec_ret_gfs[index_reconsturcted_gf] = sorted(list(focal_genes_ids))
            # Count number of singletons and number of GFs with at least two members
            # - 2 focal genes is the minimum needed to be able to compute Ks value(s)
            # - 1 focal gene can't be used and will be ignored for Ks purposes
            if species_id_dropna.value_counts().loc[species_name] >= 2:
                gfs_with_two_or_more_genes_focal += 1
            else:
                singletons_gfs_to_be_ignored += 1
    
    # Log number of GFs with 0, 1 or at least 2 focal species' genes 
    ratio_no_genes = round(gfs_with_no_gene_of_focal/num_matched_orthomcl_gfs_accepted*100, 2)
    ratio_one_gene = round(singletons_gfs_to_be_ignored/num_matched_orthomcl_gfs_accepted*100, 2)
    ratio_two_genes = round(gfs_with_two_or_more_genes_focal/num_matched_orthomcl_gfs_accepted*100, 2)
    logging.info(f"Filtering for reconstructed GFs including at least 2 focal species' paralogs (minimum required to calculate Ks values):")
    logging.info(f"- {gfs_with_two_or_more_genes_focal}/{num_matched_orthomcl_gfs_accepted} ({ratio_two_genes}%) top GFs have 2 or more genes for focal species and can be used to calculate Ks values")
    logging.info(f"- {singletons_gfs_to_be_ignored}/{num_matched_orthomcl_gfs_accepted} ({ratio_one_gene}%) top GFs have only 1 gene for focal species; will be ignored")
    logging.info(f"- {gfs_with_no_gene_of_focal}/{num_matched_orthomcl_gfs_accepted} ({ratio_no_genes}%) top GFs have 0 genes for focal species; will be ignored") 
    if gfs_with_two_or_more_genes_focal == 0:
        logging.warning("No gene families were found with at least two focal genes")
        logging.warning("Will skip Ks estimate and construction of reciprocally retained Ks distribution!")
        sys.exit(0) # exit code 0 because no actual errors were thrown
    elif gfs_with_two_or_more_genes_focal != num_matched_orthomcl_gfs_accepted:
        logging.warning(f"Note: a lower number of usable GFs than {num_matched_orthomcl_gfs_accepted} can be expected to a certain extent;")
        logging.warning(f"however, keep in mind that too low numbers may generate poor Ks distributions!")
    logging.info("")

    # Convert to DataFrame
    focal_rec_ret_gfs_df = DataFrame.from_dict(focal_rec_ret_gfs, orient="index")
    
    # Sort rows by decreasing number of focal genes
    num_focal_genes = focal_rec_ret_gfs_df.apply(lambda x: len(x.dropna()), axis=1)
    num_focal_genes.index = focal_rec_ret_gfs_df.index
    num_focal_genes_sorted = num_focal_genes.sort_values(ascending=False)
    focal_rec_ret_gfs_df_sorted = focal_rec_ret_gfs_df.reindex(num_focal_genes_sorted.index)
    
    # Retain only gene families with at least two focal genes, because singletons are not suitable for Ks estimate
    focal_rec_ret_gfs_df_sorted_for_ks = focal_rec_ret_gfs_df_sorted.dropna(thresh=2)
    
    # Generate file of reconstructed GFs with at least two 2 focal genes (MCL-like format)
    logging.info(f"Generating MCL-like file [{os.path.basename(output_mcl_path)}]")
    focal_rec_ret_gfs_df_sorted_for_ks.to_csv(output_mcl_path, sep="\t", index=None, header=None)
    
    # For transparency, generate analogous file of reconstructed GFs (MCL-like format) including GF IDs and listing also the singletons
    focal_rec_ret_gfs_df_sorted.to_csv(output_mcl_path_with_gf_ids, sep="\t", header=None)
