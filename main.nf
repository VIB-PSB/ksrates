#!/usr/bin/env nextflow

// should probably use our own Logger...
LOG_OUTPUT = true  // log process and other output via log.info
LOG_PROCESS_INFO = false  // log process details before submission

version_file = file("${workflow.projectDir}/ksrates/_version.py").readLines()[0] =~ /__version__ = "(.*)"/
version = version_file[0][1]

/*
 * Pipeline input parameters
 */

// Parameter to automatically delete or not leftover folders at the end of the pipeline
params.preserve = false

// Giving the configuration file through the "input" process section
// NOTE:
// $params.config is what the user has entered after --config in the command line, for example ./config_files/config_elaeis.txt 
// $configfile is the config file absolute path, e.g. /home/config_files/config_elaeis.txt
// $config is the basename of the config file, e.g. config_elaeis.txt
configfile = file(params.config)
if (configfile.isEmpty()) {
     newConfigFile = true
} else {
     newConfigFile = false
}

// Giving the expert configuration file as input
// NOTE:
// $params.expert is what the user has entered after --expert in the command line, for example ./config_files/config_expert.txt 
// $expert_configfile is the expert config file absolute path, e.g. /home/config_files/config_expert.txt
// $expert_config is the basename of the expert config file, e.g. config_expert.txt

// Set parameter expert to false by default (can be overwritten to something else by using --expert in command line)
params.expert = false
user_defined_expert_config_missing = false
default_expert_config_file_available = true
// If user specified the expert config file through command line (e.g. to specify a certain path), overwrite "false" and read file
if (params.expert) {
    expert_configfile = file(params.expert)
    // If user has defined a file/path, but it doesn't exists, flag it and interrupt pipeline at later checkpoint
    if (expert_configfile.exists() == false) {
        user_defined_expert_config_missing = true
    }
}
else {
    // If parameter --expert not used in command line but file config_expert.txt is actually present in the launching folder, use it
    if (file("config_expert.txt").exists()) {
        expert_configfile = file("config_expert.txt")
        // Note that default_expert_config_file_available is already set to true
    }
    // If parameter --expert not used in command line and default filename config_expert.txt doesn't exist,
    // set the parameter to an empty string so that default parameters will be used
    else {
        expert_configfile = ""
        default_expert_config_file_available = false
    }
}

// Collect input configuration files in a string that will be given to each ksrates command.
// The string always contains the standard config file, plus additionally the expert config file
// if this latter was provided by the user through the --expert option or if it was found 
// with default name "config_expert.txt".
config_args = "${configfile}"
// If expert_configfile is NOT a string, it means that it has been provided and will be used
if (expert_configfile !instanceof String) {
    config_args = config_args + " --expert ${expert_configfile}"
}

log.info ""
log.info ""
log.info """\
         K S R A T E S   -   N E X T F L O W   P I P E L I N E   (v${version})
         ----------------------------------------------------------------
         
         Configuration file:                    $params.config
         Logs folder:                           logs_${workflow.sessionId.toString().substring(0,8)}
         Preserve leftover files:               $params.preserve
         """
         .stripIndent()

log.info "Command line:               $workflow.commandLine"
log.info "Launch directory:           $workflow.launchDir"
//log.info "PWD:                        $PWD"
log.info "Work directory:             $workflow.workDir"
//log.info "Base dir:                   $baseDir"
log.info "ksrates directory:          $workflow.projectDir"
//log.info "ksrates Git repo:           $workflow.repository"
//log.info "ksrates Git revision:       $workflow.revision"
//log.info "ksrates container:          $workflow.container"
//log.info "ksrates container engine:   $workflow.containerEngine"
//log.info "Configuration files:        $workflow.configFiles"
//log.info "Session ID:                 $workflow.sessionId"
//log.info "Stats:                      $workflow.stats"
log.info ""
log.info "Start time:                 $workflow.start"
log.info ""
log.info ""

/*
 * Check if Nextflow runtime environment variable $NXF_ANSI_LOG is set,
 * if not set it to true (the default) so it is accessible within script.
 */
try {
    NXF_ANSI_LOG
} catch (Exception e) {
    NXF_ANSI_LOG = "true"
}

//log.info "NXF_ANSI_LOG: $NXF_ANSI_LOG"

// command-line option -ansi-log has preference over NXF_ANSI_LOG
if ( workflow.commandLine.containsIgnoreCase("-ansi-log true") ) {
    ansiLog = true
} else if ( workflow.commandLine.containsIgnoreCase("-ansi-log false") ) {
    ansiLog = false
} else if ( workflow.commandLine.containsIgnoreCase("-bg") ) {
    // Nextflow -bg option turns ANSI log off (unless explicit -ansi-log true)
    ansiLog = false
} else if ( NXF_ANSI_LOG.equalsIgnoreCase("false") ) {
    ansiLog = false
} else {
    // default should be true
    ansiLog = true
}
//log.info "ansiLog: $ansiLog"


// Dictionary mapping each process with the associated log file
logs_names = [
"checkConfig" : null,
"setupAdjustment" : "setup_adjustment.log",
"setParalogAnalysis" : "wgd_paralogs.log",
"setOrthologAnalysis": "set_orthologs.log",
"estimatePeaks": "estimate_peaks.log",
"wgdParalogs": "wgd_paralogs.log", 
"wgdOrthologs": "wgd_orthologs_",
"plotOrthologDistrib": "plot_ortholog_distributions.log", 
"doRateAdjustment": "rate_adjustment.log",
"paralogsAnalyses": "paralogs_analyses.log",
"drawTree": "rate_adjustment.log" ]


/*
 * Function to log various process information before submission/execution
 * if print is set to true or if any process setting deviates from defaults
 * (executor = 'local', cpus = 1, and no memory setting).
 */
def logProcessInfo(task, speciesLabel="", printIndex=false, print=LOG_PROCESS_INFO) {
    if (print || task.executor != 'local' || task.cpus != 1 || task.memory) {
        log.info "[${task.process}]${ speciesLabel ? ' [' + speciesLabel + '] ' : ' ' }" +\
                 "Process${ printIndex ? ' (' + task.index + ') ' : ' ' }" +\
                 "queued on executor '${task.executor}' with " +\
                 "${task.cpus} ${ task.cpus == 1 ? 'CPU' : 'CPUs' }" +\
                 "${ task.memory ? ' and ' + task.memory : '' }"
    }
}

/*
 * Function to log process output to stdout/console (if print is set to true).
 * Each non-empty log line will have the given processName prepended in square brackets.
 * Note that the output will be logged all together only after the process terminates.
 */
def logProcessOutput(out, processName, print=LOG_OUTPUT) {
    if (print) {
        out.splitText() {
            if ( it.equals("\n") ) {
                // remove newline since subsequent log.info will append another one
                ""
            } else if ( it.endsWith("\n") ) {
                // remove newline since subsequent log.info will append another one
                "[$processName] " + it.replace('\n', '')
            } else {
                "[$processName] $it"
            }
        }
        .subscribe {
            if ( "".equals(it) && ansiLog ) {
                // ignore, i.e. do not print empty lines if ansiLog is true
                // to make output more compact
            } else {
                log.info it.trim()
            }
        }
    }
}


////////////////////////////////////////////////////////////


/*
 * Calling ksrates workflow
 */

workflow {

    checkConfig(
        configfile,
        expert_configfile
    )
    // Log output to stdout/console (always, i.e. no matter if LOG_OUTPUT).
    logProcessOutput(checkConfig.out.outCheckConfig, "checkConfig", true)


    setupAdjustment(
        configfile,
        checkConfig.out.trigger_setupAdjustment_channel
    )
    // Log output to stdout/console.
    logProcessOutput(setupAdjustment.out.outSetupAdjustment, "setupAdjustment")


    setParalogAnalysis(
        setupAdjustment.out.species_channel,

        setupAdjustment.out.logs_folder_channel,

        configfile,
        
        expert_configfile
    )
    // Log output to stdout/console.
    logProcessOutput(setParalogAnalysis.out.outSetParalogAnalysis, "setParalogAnalysis")


    setOrthologAnalysis(
        setupAdjustment.out.check_ortholog_pairs_channel,
        
        setupAdjustment.out.logs_folder_channel,
        
        configfile
    )
    // Log output to stdout/console.
    logProcessOutput(setOrthologAnalysis.out.outSetOrthologAnalysis, "setOrthologAnalysis")


    estimatePeaks(
        setOrthologAnalysis.out.file_for_estimatePeak_channel,

        setupAdjustment.out.logs_folder_channel,

        configfile
    )
    // Log output to stdout/console.
    logProcessOutput(estimatePeaks.out.outEstimatePeaks, "estimatePeaks")


    wgdParalogs(
        setupAdjustment.out.species_channel,

        setParalogAnalysis.out.trigger_wgdPara_channel,

        setupAdjustment.out.logs_folder_channel,

        configfile
    )
    // Log output to stdout/console.
    logProcessOutput(wgdParalogs.out.outParalogs, "wgdParalogs")


    wgdOrthologs(
        setupAdjustment.out.species_channel,

        setOrthologAnalysis.out.species_pairs_for_wgd_Orthologs_channel.splitCsv(sep:'\t'),

        setupAdjustment.out.logs_folder_channel, 
        
        configfile
    )
    // Log output to stdout/console.
    logProcessOutput(wgdOrthologs.out.outOrthologs, "wgdOrthologs")


    plotOrthologDistrib(
        setupAdjustment.out.species_channel,

        setupAdjustment.out.trigger_plotOrtholog_from_setupAdjustment_channel.mix(
            estimatePeaks.out.trigger_plotOrtholog_from_estimatePeak_together_with_wgdOrthologs_channel.merge(
                wgdOrthologs.out.trigger_plotOrtholog_from_wgdOrtholog_together_with_estimatePeak_channel.collect()
                ),
            setOrthologAnalysis.out.trigger_plotOrthologs_together_with_wgdOrtholog_channel.merge(
                wgdOrthologs.out.trigger_plotOrtholog_from_wgdOrtholog_channel.collect()
                ), 
            setOrthologAnalysis.out.trigger_plotOrthologs_together_with_estimatePeak_channel.merge(
                estimatePeaks.out.trigger_plotOrtholog_from_estimatePeak_channel
                )
        ),

        setupAdjustment.out.logs_folder_channel, 

        configfile
    )
    // Log output to stdout/console.
    logProcessOutput(plotOrthologDistrib.out.outPlotOrthologDistrib, "plotOrthologDistrib")


    doRateAdjustment(
        setupAdjustment.out.species_channel,

        setParalogAnalysis.out.trigger_doRateAdjustment_from_setParalog_channel.mix(
            estimatePeaks.out.trigger_doRateAdjustment_from_estimatePeak_channel,
            wgdOrthologs.out.trigger_doRateAdjustment_from_wgdOrtholog_channel,
            wgdParalogs.out.trigger_doRateAdjustment_from_wgdParalog_channel
            ),
        
        setupAdjustment.out.logs_folder_channel, 
        
        configfile,

        expert_configfile
    )
    // Log output to stdout/console.
    logProcessOutput(doRateAdjustment.out.outDoRateAdjustment, "doRateAdjustment")
    

    paralogsAnalyses(
        setupAdjustment.out.species_channel,

        setupAdjustment.out.logs_folder_channel, 

        configfile,

        doRateAdjustment.out.trigger_paralogsAnalyses_from_doRateAdjustment_channel.collect()
    )
    // Log output to stdout/console.
    logProcessOutput(paralogsAnalyses.out.outParalogsAnalyses, "paralogsAnalyses")


    drawTree(
        setupAdjustment.out.species_channel,

        setupAdjustment.out.logs_folder_channel, 

        configfile,

        doRateAdjustment.out.trigger_drawTree_from_doRateAdjustment_channel.collect()
    )
    // Log output to stdout/console.
    logProcessOutput(drawTree.out.outDrawTree, "drawTree")

}


////////////////////////////////////////////////////////////


process checkConfig {

    executor 'local'
    maxForks 1

    input:
        file config /// from configfile DONE
        file expert_config /// from expert_configfile DONE

    output:
        stdout emit: outCheckConfig /// DONE
        env trigger_pipeline, emit: trigger_setupAdjustment_channel /// into trigger_setupAdjustment_channel DONE

    script:
    logProcessInfo(task)
    """
    processDir=\$PWD
    cd $PWD

    if [ ! -f ${config} ]; then
        echo ""
        echo -n "Configuration file [${config}] not found: it will now be generated... "

        ksrates generate-config ${config} >> \$processDir/generate_config.txt
        RET_CODE=\$?
        echo "done [\${RET_CODE}]"

        echo "Please fill in the required parameters:"
        echo " focal_species, newick_tree, latin_names, fasta_filenames, and"
        echo " gff_filename, gff_feature and gff_attribute (if applicable)"
        echo "Then rerun the Nextflow pipeline"

        trigger_pipeline=false
    else
        echo ""
        echo -n "Configuration file [${config}] found"
        trigger_pipeline=true
    fi

    # Expert config file
    # If the file at the user-defined path is not found, warn and exit
    if [ ${user_defined_expert_config_missing} = "true" ]; then
        echo ""
        echo "User-defined expert configuration file [${params.expert}] not found:"
        echo "please check the input path after the '--expert' parameter in the command line and rerun the analysis"
        exit 1
    # If the file with default name can't be found, will use default parameters 
    elif [ ${default_expert_config_file_available} = "false" ]; then
        echo ""
        echo -n "Expert configuration file [config_expert.txt] not available: will use default parameters"
    # If the file with default name is found, will use parameter therein listed
    else
        echo ""
        echo -n "Expert configuration file [${expert_config}] found"
    fi

    cd \$processDir
    """
}



/* 
 * Process that extracts from the input tree the ortholog trios and
 * the ortholog species pairs used for the rate-adjustment.
 */
process setupAdjustment {

    executor 'local'
    maxForks 1

    input:
        file config /// from configfile DONE
        val trigger_pipeline /// from trigger_setupAdjustment_channel DONE

    output:
        stdout emit: outSetupAdjustment /// DONE
        env species, emit: species_channel /// into species_channel DONE
        env logs_folder, emit: logs_folder_channel /// into logs_folder_channel DONE
        path "ortholog_pairs_*.tsv", emit: check_ortholog_pairs_channel /// into check_ortholog_pairs_channel DONE
        env trigger_plot_orthologs, emit: trigger_plotOrtholog_from_setupAdjustment_channel /// into trigger_plotOrtholog_from_setupAdjustment_channel DONE

    when:
        trigger_pipeline == "true"

    script:
    logProcessInfo(task)
    """
    echo ""
    trigger_plot_orthologs=false

    processDir=\$PWD
    cd $PWD

    species=`grep "^[[:space:]]*focal_species[[:space:]]*=" ${configfile} | cut -d "=" -f 2 | xargs`
    echo "Focal species: \$species"

    # Generating folders for output files, especially to have the log_folder since the very beginning
    if [ ! -d rate_adjustment ]; then
        mkdir rate_adjustment
    fi
    if [ ! -d rate_adjustment/\$species ]; then
        mkdir rate_adjustment/\$species
    fi    
    workID=`echo ${workflow.sessionId} | cut -c1-8`
    if [ ! -d rate_adjustment/\$species/logs_\$workID ]; then
        mkdir rate_adjustment/\$species/logs_\$workID
    fi
    logs_folder="rate_adjustment/\$species/logs_\$workID"


    echo -n "Extracting ortholog species pairs and trios from Newick tree... "
    echo "NF internal work directory for [setupAdjustment] process:\n\$processDir\n" > \${logs_folder}/${logs_names["setupAdjustment"]}

    ksrates init ${config_args} --nextflow >> \${logs_folder}/${logs_names["setupAdjustment"]} 2>&1

    RET_CODE=\$?
    echo "done [\${RET_CODE}]"

    cat rate_adjustment/\$species/ortholog_pairs_\${species}.tsv > \${processDir}/ortholog_pairs_\${species}.tsv

    # If all the ortholog data are present in the databases, already trigger plotOrthologs to plot the ortholog distributions
    if [ -z "`tail -n +2 rate_adjustment/\$species/ortholog_pairs_\${species}.tsv`" ]; then    # if file size NOT greater than 0, so if the file is empty and there aren't unknown pairs
        echo "All ortholog species pair data already present in ortholog databases"
        trigger_plot_orthologs=true
    else
        echo "The list of required ortholog species pairs can be found in:"
        echo " rate_adjustment/\$species/ortholog_pairs_\${species}.tsv"
        echo "The list of ortholog species trios can be found in:"
        echo " rate_adjustment/\$species/ortholog_trios_\${species}.tsv"
    fi

    echo "Log can be found in:"
    echo " \${logs_folder}/${logs_names["setupAdjustment"]}"

    cd \$processDir
    """
}



/*
 * Process that checks if the .ks.tsv file(s) containing the paralog Ks values
 * for the focal species are already present. If not, triggers
 * wgdParalogs for their estimate.
 */ 
process setParalogAnalysis {

    executor 'local'
    maxForks 1

    input:
        val species /// from species_channel
        val logs_folder /// from logs_folder_channel
        file config /// from configfile
        file expert_config /// from expert_configfile

    output:
        stdout emit: outSetParalogAnalysis /// DONE
        env trigger_wgdPara, emit: trigger_wgdPara_channel /// into trigger_wgdPara_channel DONE
        env trigger_doRateAdjustment_from_para, emit: trigger_doRateAdjustment_from_setParalog_channel /// into trigger_doRateAdjustment_from_setParalog_channel DONE

    script:
    logProcessInfo(task)
    """
    echo ""
    echo "[$species] Setting up paralog wgd analysis for focal species"

    trigger_doRateAdjustment_from_para=false
    trigger_wgdPara=false

    paranome_status="not_required"
    colinearity_status="not_required"

    paranome=`grep "^[[:space:]]*paranome[[:space:]]*=" ${config} | cut -d "=" -f 2 | xargs | tr '[:upper:]' '[:lower:]'`
    colinearity=`grep "^[[:space:]]*collinearity[[:space:]]*=" ${config} | cut -d "=" -f 2 | xargs | tr '[:upper:]' '[:lower:]'`

    processDir=\$PWD
    cd $PWD

    echo "NF internal work directory for [setParalogAnalysis] process:\n\$processDir\n" > $logs_folder/${logs_names["setParalogAnalysis"]}

    # Triggering wgdParalog process only if ".ks.tsv" (and ".ks_anchors.tsv") files are missing

    if [ \${paranome} = "yes" ]; then
        if [ ! -f paralog_distributions/wgd_${species}/${species}.ks.tsv ]; then
            echo "[$species] Will run whole-paranome wgd analysis"
            echo "[$species] Paralog TSV file not found [${species}.ks.tsv]" >> $logs_folder/${logs_names["setParalogAnalysis"]}
            echo "[$species] Whole-paranome wgd pipeline will be started\n" >> $logs_folder/${logs_names["setParalogAnalysis"]}
            paranome_status="todo"
        else
            paranome_status="already_done"
        fi
    fi
    if [ \${colinearity} = "yes" ]; then
        if [ ! -f paralog_distributions/wgd_${species}/${species}.ks_anchors.tsv ]; then
            echo "[$species] Will run anchor-pair (collinearity) wgd analysis"
            echo "[$species] Anchor pairs TSV file not found [${species}.ks_anchors.tsv]" >> $logs_folder/${logs_names["setParalogAnalysis"]}
            echo "[$species] Collinearity wgd pipeline will be started\n" >> $logs_folder/${logs_names["setParalogAnalysis"]}
            colinearity_status="todo"
        else
            colinearity_status="already_done"
        fi
    fi

    if [ \$paranome_status = "todo" ] || [ \$colinearity_status = "todo" ]; then
        # Trigger wgdParalog process to get the missing tsv file(s)
        trigger_wgdPara=true
    else
        echo "[$species] Paralog data already present; skipping paralog wgd analysis"
        echo "[$species] Paralog TSV file(s) already present; skipping paralog wgd pipeline\n" >> $logs_folder/${logs_names["setParalogAnalysis"]}
        
        # Trigger doRateAdjustment process to plot (at least) the paralog distribution in the mixed plot
        trigger_doRateAdjustment_from_para=true
    fi

    echo "----------------------------------------------------------------\n" >> $logs_folder/${logs_names["setParalogAnalysis"]}
    echo "[$species] Log can be found in:"
    echo "[$species]  $logs_folder/${logs_names["setParalogAnalysis"]}"

    cd \$processDir
    """
}



/*
 * Process that receives the list of all species pairs missing in database and that 
 * 1) triggers wgdOrthologs for pairs without the ks.tsv file and/or 2) triggers estimatePeaks
 * for pairs without ortholog peak data in database but that do have the .ks.tsv file.
 */
process setOrthologAnalysis {

    executor 'local'
    maxForks 1

    input:
        file ortholog_pairs /// from check_ortholog_pairs_channel DONE
        val logs_folder /// from logs_folder_channel DONE
        file config /// from configfile DONE

    output:
        stdout emit: outSetOrthologAnalysis /// DONE
        path "tmp_species_pairs_for_wgdOrtholog.txt", optional: true, emit: species_pairs_for_wgd_Orthologs_channel /// into species_pairs_for_wgd_Orthologs_channel DONE
        path "tmp_species_pairs_for_estimatePeak.txt", optional: true, emit: file_for_estimatePeak_channel /// into file_for_estimatePeak_channel DONE
        env estimatePeak_not_needed, emit: trigger_plotOrthologs_together_with_wgdOrtholog_channel /// into trigger_plotOrthologs_together_with_wgdOrtholog_channel DONE
        env wgdOrtholog_not_needed, emit: trigger_plotOrthologs_together_with_estimatePeak_channel /// into trigger_plotOrthologs_together_with_estimatePeak_channel DONE

    script:
    logProcessInfo(task)
    """
    echo ""
    echo "Setting up required ortholog analyses"

    processDir=\$PWD
    cd $PWD

    echo "NF internal work directory for [setOrthologAnalysis] process:\n\$processDir\n" >> $logs_folder/${logs_names["setOrthologAnalysis"]}

    wgdOrtholog_not_needed=false
    estimatePeak_not_needed=false

    # If the file with the ortholog species pairs is empty, then Ks estimation and peak estimation steps are skipped.
    if [ -z "`tail -n +2 \${processDir}/$ortholog_pairs`" ]; then
        echo "No ortholog species pairs listed; skipping ortholog wgd and peak analyses"
        echo "No species pairs are listed for wgd ortholog pipeline or ortholog peak estimate." >> $logs_folder/${logs_names["setOrthologAnalysis"]}
    fi 

    while read -r species1 species2 || [ -n "\$species1" ]; do
        if [ ! -f ortholog_distributions/wgd_\${species1}_\${species2}/\${species1}_\${species2}.ks.tsv ]; then
            echo "\$species1\t\$species2" >> \${processDir}/tmp_species_pairs_for_wgdOrtholog.txt
            echo "[\$species1 – \$species2] Will run ortholog wgd analysis"
            echo "[\$species1 – \$species2] Ortholog TSV file not present [\${species1}_\${species2}.ks.tsv]" >> $logs_folder/${logs_names["setOrthologAnalysis"]}
            echo "[\$species1 – \$species2] Ortholog wgd pipeline will be started" >> $logs_folder/${logs_names["setOrthologAnalysis"]}
        else
            echo "[\$species1 – \$species2] Will run ortholog peak analysis"
            echo "[\$species1 – \$species2] Ortholog TSV file already present; skipping ortholog wgd pipeline" >> $logs_folder/${logs_names["setOrthologAnalysis"]}
            echo "[\$species1 – \$species2] Update of ortholog database(s) will be started." >> $logs_folder/${logs_names["setOrthologAnalysis"]}

            # if the tmp file does not exist yet, add headers
            if [ ! -f \${processDir}/tmp_species_pairs_for_estimatePeak.txt ]; then
                echo "Species1\tSpecies2" >> \${processDir}/tmp_species_pairs_for_estimatePeak.txt
            fi
            echo "\$species1\t\$species2" >> \${processDir}/tmp_species_pairs_for_estimatePeak.txt   # add a species pair
        fi
        echo " \n----------------------------------------------------------------\n" >> $logs_folder/${logs_names["setOrthologAnalysis"]}
    done < <(tail -n +2 \${processDir}/$ortholog_pairs) # skipping the headers

    # if the file for wgdOrthologs is empty (size not > 0)
    if [ ! -s \${processDir}/tmp_species_pairs_for_wgdOrtholog.txt ]; then
        wgdOrtholog_not_needed=true
    else
        if [ ! -d ortholog_distributions ]; then
            mkdir ortholog_distributions
        fi
    fi

    # if the file content for estimatePeaks (skipping headers) is empty (empty string check)
    if [ -z "`tail -n +2 \${processDir}/tmp_species_pairs_for_estimatePeak.txt`" ]; then
        estimatePeak_not_needed=true
    fi

    echo "Log can be found in:"
    echo " $logs_folder/${logs_names["setOrthologAnalysis"]}"

    cd \$processDir
    """
}



/*
 * Process that estimates the ortholog Ks distribution peak for species pairs
 * that have already their .ks.tsv files, but that for some reasons
 * are not present anymore in the ortholog database(s) (i.e. have been deleted).
 */
process estimatePeaks {

    input:
        file species_pairs_for_peak /// from file_for_estimatePeak_channel
        val logs_folder /// from logs_folder_channel
        file config /// from configfile
        
    output:
        stdout emit: outEstimatePeaks /// DONE
        val true, emit: trigger_doRateAdjustment_from_estimatePeak_channel /// into trigger_doRateAdjustment_from_estimatePeak_channel DONE
        val true, emit: trigger_plotOrtholog_from_estimatePeak_channel /// into trigger_plotOrtholog_from_estimatePeak_channel DONE
        val true, emit: trigger_plotOrtholog_from_estimatePeak_together_with_wgdOrthologs_channel /// into trigger_plotOrtholog_from_estimatePeak_together_with_wgdOrthologs_channel DONE

    script:
    logProcessInfo(task)
    """
    echo ""
    echo -n "`date "+%T"` Starting missing ortholog peak analyses... "

    processDir=\$PWD
    cd $PWD
    echo "NF internal work directory for [estimatePeaks (${task.index})] process:\n\$processDir\n" >> $logs_folder/${logs_names["estimatePeaks"]}

    echo "Updating ortholog peak database" >> $logs_folder/${logs_names["estimatePeaks"]}

    ksrates orthologs-analysis ${config_args} --ortholog-pairs=\$processDir/$species_pairs_for_peak >> $logs_folder/${logs_names["estimatePeaks"]} 2>&1

    RET_CODE=\$?
    echo "done [\${RET_CODE}] `date "+%T"`"

    echo "Log can be found in:"
    echo " $logs_folder/${logs_names["estimatePeaks"]}"

    cd \$processDir
    """
}



/*
 * Process that estimates the paranome and/or anchor pair Ks values
 * for the focal species.
 */
process wgdParalogs {

    input:
        val species /// from species_channel
        val trigger_wgdPara /// from trigger_wgdPara_channel
        val logs_folder /// from logs_folder_channel
        file config /// from configfile
        
    output:
        stdout emit: outParalogs /// DONE
        val true, emit: trigger_doRateAdjustment_from_wgdParalog_channel /// into trigger_doRateAdjustment_from_wgdParalog_channel DONE

    when:
        trigger_wgdPara == "true"

    script:
    logProcessInfo(task, species, printIndex=false, print=true)
    """
    echo ""
    echo -n "[$species] `date "+%T"` Starting paralog wgd analysis... "

    processDir=\$PWD
    cd $PWD
    echo "NF internal work directory for [wgdParalogs (${task.index})] process:\n\$processDir\n" >> $logs_folder/${logs_names["wgdParalogs"]}

    echo "Using ${task.cpus} thread(s)\n">> $logs_folder/${logs_names["wgdParalogs"]}

    ksrates paralogs-ks ${config_args} --n-threads=${task.cpus} >> $logs_folder/${logs_names["wgdParalogs"]} 2>&1

    RET_CODE=\$?
    echo "done [\${RET_CODE}] `date "+%T"`"

    echo "[$species] wgd log can be found in:"
    echo "[$species]  $logs_folder/${logs_names["wgdParalogs"]}"

    cd \$processDir
    """
}



/*
 * Process that estimates the ortholog Ks values for a species pair and
 * estimates also the ortholog Ks distribution peak.
 */
process wgdOrthologs {

    input:
        val species /// from species_channel DONE
        tuple val(species1), val(species2) /// from species_pairs_for_wgd_Orthologs_channel.splitCsv(sep:'\t') DONE
        val logs_folder /// from logs_folder_channel DONE
        file config /// from configfile DONE
        
    output:
        stdout emit: outOrthologs /// DONE
        val true , emit: trigger_doRateAdjustment_from_wgdOrtholog_channel /// into trigger_doRateAdjustment_from_wgdOrtholog_channel DONE
        val true , emit: trigger_plotOrtholog_from_wgdOrtholog_channel /// into trigger_plotOrtholog_from_wgdOrtholog_channel DONE
        val true , emit: trigger_plotOrtholog_from_wgdOrtholog_together_with_estimatePeak_channel /// into trigger_plotOrtholog_from_wgdOrtholog_together_with_estimatePeak_channel DONE

    script:
    logProcessInfo(task, "$species1 – $species2", printIndex=true, print=true)
    """
    echo ""
    echo -n "[$species1 – $species2] `date "+%T"` Starting ortholog wgd analysis... "

    processDir=\${PWD}
    cd $PWD
    echo "NF internal work directory for [wgdOrthologs (${task.index})] process:\n\$processDir\n" > $logs_folder/${logs_names["wgdOrthologs"]}${species1}_${species2}.log

    echo "Using ${task.cpus} thread(s)\n">> $logs_folder/${logs_names["wgdOrthologs"]}${species1}_${species2}.log

    ksrates orthologs-ks ${config_args} $species1 $species2 --n-threads=${task.cpus} >> $logs_folder/${logs_names["wgdOrthologs"]}${species1}_${species2}.log 2>&1

    RET_CODE=\$?
    echo "done [\${RET_CODE}] `date "+%T"`"
    
    echo -n "[$species1 – $species2] `date "+%T"` Starting ortholog peak analysis... "
    echo "\n" >> $logs_folder/${logs_names["wgdOrthologs"]}${species1}_${species2}.log
    echo "Species1\tSpecies2\n$species1\t$species2" > \${processDir}/tmp_${species1}_${species2}.txt
    
    ksrates orthologs-analysis ${config_args} --ortholog-pairs=\${processDir}/tmp_${species1}_${species2}.txt >> $logs_folder/${logs_names["wgdOrthologs"]}${species1}_${species2}.log 2>&1

    RET_CODE=\$?
    echo "done [\${RET_CODE}] `date "+%T"`"

    echo "[$species1 – $species2] wgd log and peak analysis log can be found in:"
    echo "[$species1 – $species2]  $logs_folder/${logs_names["wgdOrthologs"]}${species1}_${species2}.log"

    cd \$processDir
    """
}



/*
 * Process that generates figures collecting all the ortholog distributions 
 * used for the adjustment and that highlights their peaks.
 */
process plotOrthologDistrib {

    maxForks 1

    input:
        val species /// from species_channel
        val trigger /// from trigger_plotOrtholog_from_setupAdjustment_channel.mix(trigger_plotOrtholog_from_estimatePeak_together_with_wgdOrthologs_channel.merge(trigger_plotOrtholog_from_wgdOrtholog_together_with_estimatePeak_channel.collect()), trigger_plotOrthologs_together_with_wgdOrtholog_channel.merge(trigger_plotOrtholog_from_wgdOrtholog_channel.collect()), trigger_plotOrthologs_together_with_estimatePeak_channel.merge(trigger_plotOrtholog_from_estimatePeak_channel))
        val logs_folder /// from logs_folder_channel
        file config /// from configfile
        
    output:
        stdout emit: outPlotOrthologDistrib /// DONE

    when:
        trigger == "true" || !trigger.contains("false")
        /*
         * String "false" or string "true" comes from setupAdjustment; to be an accepted trigger must be string "true"
         * [true, many true] comes from trigger_plotOrtholog_from_estimatePeak_together_with_wgdOrthologs_channel.merge(trigger_plotOrtholog_from_wgdOrtholog_together_with_estimatePeak_channel.collect()
         *    it is always a trigger because it always doesn't contain the string "false"
         * ["false"/"true", many true] from trigger_plotOrthologs_together_with_wgdOrtholog_channel.merge(trigger_plotOrtholog_from_wgdOrtholog_channel.collect())
         *    to be a trigger must be ["true", many true], so it must not contain string "false"
         * ["false"/"true", many true] from trigger_plotOrthologs_together_with_estimatePeak_channel.merge(trigger_plotOrtholog_from_estimatePeak_channel)
         *    to be a trigger must be ["true", many true], so it must not contain string "false"
         *
         * An accepted trigger comes from setupAdjustment if all ortholog data are already present and ready to be plotted; from estimatePeaks together with wgdOrthologs
         * when they have both finished with all the peaks; from setOrthologs and wgdOrthologs after this latter has done with all ortholog distributions;
         * from setOrtholog and estimatePeaks when all missing peaks are computed.
         */

    script:
    logProcessInfo(task, speciesLabel="", printIndex=true)
    """
    echo ""
    echo -n "`date "+%T"` Plotting ortholog distributions... "

    processDir=\$PWD
    cd $PWD
    echo "NF internal work directory for [plotOrthologDistrib] process:\n\$processDir\n" > $logs_folder/${logs_names["plotOrthologDistrib"]}

    ksrates plot-orthologs ${config_args} >> $logs_folder/${logs_names["plotOrthologDistrib"]} 2>&1
    
    RET_CODE=\$?
    echo "done [\${RET_CODE}] `date "+%T"`"

    echo "Log can be found in:"
    echo " $logs_folder/${logs_names["plotOrthologDistrib"]}"

    cd \$processDir
    """
}



/*
 * Process that performs the rate-adjustment and generates the
 * mixed distribution plots.
 */
process doRateAdjustment {

    executor 'local'
    maxForks 1

    input:
        val species /// from species_channel
        val trigger /// from trigger_doRateAdjustment_from_setParalog_channel.mix(trigger_doRateAdjustment_from_estimatePeak_channel, trigger_doRateAdjustment_from_wgdOrtholog_channel, trigger_doRateAdjustment_from_wgdParalog_channel)
        val logs_folder /// from logs_folder_channel
        file config /// from configfile
        file expert_config /// from expert_configfile
        
    output:
        stdout emit: outDoRateAdjustment
        val true, emit: trigger_paralogsAnalyses_from_doRateAdjustment_channel /// into trigger_paralogsAnalyses_from_doRateAdjustment_channel
        val true, emit: trigger_drawTree_from_doRateAdjustment_channel /// into trigger_drawTree_from_doRateAdjustment_channel

    when:
        trigger == "true" || trigger == true
        /*
         * String "true" comes from setParalog process, which returns either a string "true" or a string "false" (it uses "env" variable, not a boolean);
         * Boolean true come from the other triggers, which return a boolean "val" true (trigger_doRateAdjustment_from_estimatePeak_channel, trigger_doRateAdjustment_from_wgdOrtholog_channel, trigger_doRateAdjustment_from_wgdParalog_channel)
         *
         * An accepted trigger can come from setParalog (if paralog data is already present it returns string "true", so that at least the paralog distribution will be plotted),
         *     or from estimatePeaks (updates the databases and sends a boolean true) or from wgdOrthologs (when finishes a single run returns boolean true) or finally
         *     from wgdParalogs (when finishes returns boolean true).
         */

    script:
    logProcessInfo(task, speciesLabel="", printIndex=true)
    """
    echo ""
    paranome=`grep "^[[:space:]]*paranome[[:space:]]*=" ${config} | cut -d "=" -f 2 | xargs | tr '[:upper:]' '[:lower:]'`
    colinearity=`grep "^[[:space:]]*collinearity[[:space:]]*=" ${config} | cut -d "=" -f 2 | xargs | tr '[:upper:]' '[:lower:]'`

    missing_paranome=false
    if [ \${paranome} = "yes" ] && [ ! -s $PWD/paralog_distributions/wgd_${species}/${species}.ks.tsv ]; then
        missing_paranome=true
    fi
    missing_anchorpairs=false
    if [ \${colinearity} = "yes" ] && [ ! -s $PWD/paralog_distributions/wgd_${species}/${species}.ks_anchors.tsv ]; then
        missing_anchorpairs=true
    fi

    if [ \${missing_paranome} = "true" ]; then
        echo "Whole-paranome data not yet available; skipping rate-adjustment"
        exit 0
    fi
    if [ \${missing_anchorpairs} = "true" ]; then
        echo "Anchor-pair data not yet available; skipping rate-adjustment"
        exit 0
    fi

    processDir=\$PWD
    cd $PWD
    echo "NF internal work directory for [doRateAdjustment (${task.index})] process:\n\$processDir\n" >> $logs_folder/${logs_names["doRateAdjustment"]}

    echo -n "`date "+%T"` Starting rate-adjustment analysis... "

    ksrates orthologs-adjustment ${config_args} >> $logs_folder/${logs_names["doRateAdjustment"]} 2>&1

    RET_CODE=\$?
    echo "done [\${RET_CODE}] `date "+%T"`"

    echo "\n" >> $logs_folder/${logs_names["doRateAdjustment"]}

    echo -n "`date "+%T"` Plotting mixed distributions... "

    ksrates plot-paralogs ${config_args} >> $logs_folder/${logs_names["doRateAdjustment"]} 2>&1
    
    RET_CODE=\$?
    echo "done [\${RET_CODE}] `date "+%T"`"

    echo "\n-----------------------------------------------------------\n" >> $logs_folder/${logs_names["doRateAdjustment"]}

    echo "Log can be found in:"
    echo " $logs_folder/${logs_names["doRateAdjustment"]}"

    cd \$processDir
    """
}



/*
 * Process that performs the WGM peak calling and generates the 
 * related figures (anchor pair Ks clustering and mixture model).
 * It is performed only once as last process after all the doRateAdjustment's,
 * when the adjustment table is certainly complete.
 */
process paralogsAnalyses {

    input:
        val species /// from species_channel
        val logs_folder /// from logs_folder_channel
        file config /// from configfile
        val trigger /// from trigger_paralogsAnalyses_from_doRateAdjustment_channel.collect()
        
    output:
        stdout emit: outParalogsAnalyses /// DONE

    script:
    logProcessInfo(task, species)
    """
    echo ""
    echo -n "[$species] `date "+%T"` Starting detection of potential WGD signatures... "

    processDir=\$PWD
    cd $PWD
    echo "NF internal work directory for [paralogsAnalyses (${task.index})] process:\n\$processDir\n" >> $logs_folder/${logs_names["paralogsAnalyses"]}

    ksrates paralogs-analyses ${config_args} >> $logs_folder/${logs_names["paralogsAnalyses"]} 2>&1
 
    RET_CODE=\$?
    echo "done [\${RET_CODE}] `date "+%T"`"

    echo "\n" >> $logs_folder/${logs_names["paralogsAnalyses"]}

    echo "[$species] Log can be found in:"
    echo "[$species]  $logs_folder/${logs_names["paralogsAnalyses"]}"

    cd \$processDir
    """
}



/*
 * Process that plots the input tree with branch length equal to 
 * the Ks distances.
 * It is performed only once as last process after all the doRateAdjustment's,
 * when the adjustment table is certainly complete.
 */
process drawTree {

    executor 'local'
    maxForks 1

    input:
        val species /// from species_channel
        val logs_folder /// from logs_folder_channel
        file config /// from configfile
        val trigger /// from trigger_drawTree_from_doRateAdjustment_channel.collect()
        
    output:
        stdout emit: outDrawTree /// DONE

    script:
    logProcessInfo(task)
    """
    echo ""
    echo -n "`date "+%T"` Plotting phylogram with Ks branch lengths... "

    processDir=\$PWD
    cd $PWD
    echo "NF internal work directory for [drawTree] process:\n\$processDir\n" >> $logs_folder/${logs_names["drawTree"]}

    ksrates plot-tree ${config_args} --nextflow >> $logs_folder/${logs_names["drawTree"]} 2>&1

    RET_CODE=\$?
    echo "done [\${RET_CODE}] `date "+%T"`"

    echo "Log can be found in:"
    echo " $logs_folder/${logs_names["drawTree"]}"

    cd \$processDir
    """
}



/*
 * On completion of workflow clean up any temporary files left behind.
 */
workflow.onComplete {
    if (LOG_OUTPUT) {
        log.info ""
        log.info ""
        // if a ksrates config file needed to be generated (newConfigFile == true,
        // see also process checkConfig) then there won't be anything to clean up
        if ( params.preserve == false && !newConfigFile) {
            log.info "Cleaning up any temporary files left behind..."

            focal_species_line = file("${configfile}").readLines().find{ it =~ /focal_species/ }
            if ( focal_species_line.strip().split("=").size() == 2 ) {
                focal_species = focal_species_line.strip().split("=")[1].strip()
            }
            else {
                log.info "Parameter focal_species in configuration file is empty: possible leftovers in paralogs_distributions won't be deleted"
            }

            paralog_dir_path = "${workflow.launchDir}/paralog_distributions/wgd_${focal_species}*"
            ortholog_dir_path = "${workflow.launchDir}/ortholog_distributions/wgd_*"

            // Clean both paralog and ortholog directories
            for ( dir_path : [ paralog_dir_path, ortholog_dir_path ] ) {
                file(dir_path, type: "dir")
                    .each { wgd_dir ->
                        // Remove BLAST temporary folder, if any
                        file("${wgd_dir}/*.blast_tmp", type: "dir")
                            .each { tmp ->
                                if ( tmp.exists() == true ) {
                                    result = tmp.deleteDir();
                                    log.info result ? "Deleted: ${tmp}" : "Can't delete: ${tmp}"
                                    // Remove associated incomplete BLAST TSV file
                                    file("${wgd_dir}/*.blast.tsv", type: "file")
                                        .each { tsv ->
                                            if ( tsv.exists() == true ) {
                                                result = tsv.delete()
                                                log.info result ? "Deleted: ${tsv}" : "Can't delete: ${tsv}"
                                            }
                                        }
                                }
                            }
                        // Remove Ks temporary directory, if any
                        file("${wgd_dir}/*.ks_tmp", type: "dir")
                            .each { tmp ->
                                if ( tmp.exists() == true ) {
                                    result = tmp.deleteDir();
                                    log.info result ? "Deleted: ${tmp}" : "Can't delete: ${tmp}"
                                }
                            }
                        // Remove i-ADHoRe temporary directory, if any (applicable only to paralog_distributions)
                        file("${wgd_dir}/*.ks_anchors_tmp", type: "dir")
                            .each { tmp ->
                                if ( tmp.exists() == true ) {
                                    result = tmp.deleteDir();
                                    log.info result ? "Deleted: ${tmp}" : "Can't delete: ${tmp}"
                                }
                            }
                        // Remove paralog and ortholog sub-directories if they ended up being empty
                        if ( (wgd_dir.list() as List).empty == true ) {
                            result = wgd_dir.deleteDir();
                            log.info result ? "Deleted: ${wgd_dir}" : "Can't delete: ${wgd_dir}"
                        }
                    }
            }

            // Remove "core" files generated when the submitted job doesn't have enough memory
            file("${workflow.launchDir}/core.[0-9]*")
                .each { core ->
                    if ( core.exists() == true ) {
                        result = core.delete();
                        log.info result ? "Deleted: ${core}" : "Can't delete: ${core}"
                    }
                }

            log.info "Done"
            log.info ""
        }

        log.info "Logs folder: logs_${workflow.sessionId.toString().substring(0,8)}"
        log.info ""
        log.info "Pipeline completed at: $workflow.complete taking $workflow.duration"
        log.info "Execution status: ${ workflow.success ? 'OK' : 'failed' }"
//        log.info ""
//        log.info "Project:    $workflow.projectDir"
//        log.info "Cmd line:   $workflow.commandLine"
//        log.info "Session ID: $workflow.sessionId"
//        log.info "Launch dir: $workflow.launchDir"
//        log.info "Work dir:   $workflow.workDir"
//        log.info "Stats:      $workflow.stats"
        log.info ""
    }
}


/*
 * Print a error box summing up the error message
 */
workflow.onError {
    // If user interrupts with ctrl+C
    if ( workflow.errorReport == "SIGINT" || workflow.errorMessage == "SIGINT" ) {
        error_box_sigint = "Pipeline interrupted by user (${workflow.errorMessage} signal)"
        log.error " "
        log.error "${'=' * error_box_sigint.length()}\n" + \
                  "${error_box_sigint}\n" + \
                  "${'=' * error_box_sigint.length()}"
    }

    // Else if pipeline interrupted by an error
    else {
        // Define the stopped process name
        for ( name : logs_names ) {
            if ( workflow.errorReport.split("\n").find{ it =~ name.key } ) {
                process = name.key
                break;
            }
        }

        // Define the focal species name
        focal_species_line = file("${configfile}").readLines().find{ it =~ /focal_species/ }
        if ( focal_species_line.strip().split("=").size() == 2 ) {
            focal_species = focal_species_line.strip().split("=")[1].strip()
        }
        else {
            focal_species = ""
        }

        // For process wgdOrthologs, the log filename depends on the two species names,
        // which are parsed from the errorReport under "Command executed" (this latter
        // prints the process "script" section where the two names appear)
        unknown_species_pair = false
        if ( process == ("wgdOrthologs") ) {
            // If there is a match in errorReport, use that to retrieve the names
            species_pair = workflow.errorReport =~ /\[(\w+?) – (\w+?)\]/
            if ( species_pair != null && species_pair.size() != 0 ) {
                (full, species1, species2) = species_pair[0]
            }
            // If errorReport is not useful, apply the same for errorMessage
            if ( species_pair == null || species_pair.size() == 0 || species1 == null || species1 == "" || species2 == null || species2 == "" ) {
                species_pair = workflow.errorMessage =~ /\[(\w+?) – (\w+?)\]/
                if ( species_pair != null && species_pair.size() != 0 ) {
                    (full, species1, species2) = species_pair[0]
                }
            }
            // If species1 and species2 are still unknown, will point to a generic "wgd_orthologs_" later
            if ( species1 != null && species1 != "" && species2 != null && species2 != "" ) {
                logs_names["wgdOrthologs"] = "${logs_names["wgdOrthologs"]}${species1}_${species2}.log"
            }
            else {
                unknown_species_pair = true
            }
        }

        // Define the stopped process log filename
        focal_species_dir = focal_species != "" ? "${focal_species}/" : ""
        log_filename = "rate_adjustment/${focal_species_dir}" + \
                       "logs_${workflow.sessionId.toString().substring(0,8)}/${logs_names[process]}"

        // Parse the "Caused by:" line from errorReport
        // Find in errorReport the line matching "Caused by:" and get its index
        index_cause = "${workflow.errorReport}".split("\n").findIndexOf{ it =~ /Caused by:/ }
        if ( index_cause != -1 ) {
            // Get the successive line, which contains the related message
            // If not defined, caused_by_line will be null
            caused_by_line = "${workflow.errorReport.split("\n")[index_cause + 1].trim()}\n"
        }

        // Print separator to highlight the beginning of the error box; print headline
        headline = "The pipeline terminated during process '${process}' with the following error message:"
        error_box = "${'=' * headline.length()}\n" + \
                    "${headline}\n\n"                     

        log_file = file("${log_filename}")
        if ( log_file.exists() == true ) {
            // If the process log file exists and there are "ERROR" lines,
            // print such lines (the error was caught by the script)
            if ( log_file.readLines().findAll{ it =~ /ERROR/ }.size() != 0 ) {
                for ( line : log_file.readLines().findAll { it =~ /ERROR/ } ) {
                        error_box += "${line}\n"
                }
                
                // Point to the complete output of the stopped process log file
                error_box += "\nMore details may be found in the following log file:\n" + \
                             "${log_filename}\n"
            }
            // If the process log file exists but there are no "ERROR" lines,
            // print a message from errorReport and errorMessage (it's an unexpected error)
            else {
                // Print the "Caused by:" line from errorReport
                if ( caused_by_line != null ) {
                    error_box += caused_by_line
                }
                // Print errorMessage, if any
                if ( workflow.errorMessage != null ) {
                    error_box += "${workflow.errorMessage}\n"
                }

                // Point to errorReport and to the complete output of the stopped process log file
                error_box += "\nMore details may be found in the error report above, in ./.nextflow.log, " + \
                             "and/or in the \n" + \
                             "following log file: ${log_filename}\n"
            }
        }

        // Else if the process log file doesn't exist, investigate the cause in errorReport and errorMessage
        else {
            // Print the "Caused by:" line from errorReport
            if ( caused_by_line != null ) {
                error_box += caused_by_line
            }
            // Print errorMessage, if any
            if ( workflow.errorMessage != null ) {
                error_box += "${workflow.errorMessage}\n"
            }

            // Point to the complete Nextflow errorReport, to nextflow.log (and to the wgd_ortholog_ log file)
            if ( unknown_species_pair == false ) {
                error_box += "\nMore details may be found in the error report above and/or in ./.nextflow.log\n"
            }
            else {
                error_box += "\nMore details may be found in the error report above, in ./.nextflow.log, and/or in one of\n" + \
                             "the log files named rate_adjustment/${focal_species}/" + \
                             "logs_${workflow.sessionId.toString().substring(0,8)}/${logs_names["wgdOrthologs"]} (if any exist)\n"
            }
        }

        // Separator to highlight the end of the error box
        error_box += "${'=' * headline.length()}"

        log.error " "
        log.error "${error_box}"
    }
}
