#!/usr/bin/env nextflow

LOG = true  // should probably use our own Logger...

/*
 * Pipeline input parameters
 */

// Parameter to automatically delete or not leftover folders at the end of the pipeline
params.preserve = false

// giving the configuration file through the "input" process section
configfile = file(params.config)

log.info ""
log.info """\
         K S R A T E S   -   N E X T F L O W   P I P E L I N E
         -----------------------------------------------------
         
         Configuration file:                    $params.config
         Logs folder:                           logs_${workflow.sessionId.toString().substring(0,8)}
         Preserve leftover files:               $params.preserve
         """
         .stripIndent()
log.info ""

log.info "Cmd line:       $workflow.commandLine"
log.info "Launch dir:     $workflow.launchDir"
log.info "Work dir:       $workflow.workDir"
log.info "Project:        $workflow.projectDir"
//log.info "PWD:            $PWD"
//log.info "Base dir:       $baseDir"
//log.info "Logs folder:    logs_${workflow.sessionId.toString().substring(0,8)}"
// log.info "Session ID:     $workflow.sessionId"
// log.info "Stats:          $workflow.stats"
log.info ""

/*
 * Check if Nextflow runtime environment variable $NXF_ANSI_LOG is set,
 * if not set it to true (the default) so it is accessible within script.
 */
try {
    NXF_ANSI_LOG
} catch (Exception e) {
    NXF_ANSI_LOG = true
}

// Dictionary mapping each process with the associated log file
logs_names = [
"checkConfig" : null,
"setupAdjustment" : "setup_adjustment.log",
"setParalogAnalysis" : "wgd_paralogs.log",
"setOrthologAnalysis": "set_orthologs.log",
"estimatePeak": "estimate_peak.log", 
"wgdParalogs": "wgd_paralogs.log", 
"wgdOrthologs": "wgd_orthologs_",
"plotOrthologDistrib": "plot_ortholog_distributions.log", 
"doRateAdjustment": "rate_adjustment.log",
"paralogsAnalyses": "paralogs_analyses.log",
"drawTree": "rate_adjustment.log" ]


process checkConfig {

    executor 'local'
    maxForks 1

    input:
        file config from configfile
    
    output:
        stdout outcheckConfig
        env trigger_pipeline into trigger_setupCorrection_channel

    script:
    """
    processDir=\$PWD
    cd $PWD

    if [ ! -f ${config} ]; then
        echo "Configuration file [${config}] not found: it will be now generated"
        echo "Please fill in with the required parameters:"
        echo "species, newick_tree, latin_names, fasta_filenames and if applicable gff_filename, gff_feature and gff_attribute"
        echo "Then rerun the Nextflow pipeline."
        ksrates generate-config ${config} >> \$processDir/generate_config.txt
        trigger_pipeline=false
    else
        trigger_pipeline=true
    fi
    cd \$processDir
    """
}

/*
 * Log output to stdout/console.
 * (output will be logged only after process checkConfig terminates)
 */
if (LOG) {
    outcheckConfig
        .splitText() {
            if ( it.endsWith("\n") ) {
                "[checkConfig] " + it.replace('\n', '')
            } else {
                "[checkConfig] $it"
            }
        }
        .subscribe { log.info it }
}



/* 
 * Process that extracts from the input tree the ortholog trios and
 * the ortholog species pairs used for the rate-adjustment.
 */
process setupAdjustment {

    executor 'local'
    maxForks 1

    input:
        file config from configfile
        val trigger_pipeline from trigger_setupCorrection_channel

    output:
        stdout outsetupAdjustment
        env species into species_channel
        env logs_folder into logs_folder_channel
        file "ortholog_pairs_*.tsv" into check_ortholog_pairs_channel
        env trigger_plot_orthologs into trigger_plotOrtholog_from_setupCorrection_channel

    when:
        trigger_pipeline == "true"

    script:
    if (LOG) {
        log.info "Process setupAdjustment (${task.index}) has " +\
                 "${task.cpus} ${ task.cpus == 1 ? 'CPU' : 'CPUs' } available."
    }
    """
    trigger_plot_orthologs=false

    processDir=\$PWD
    cd $PWD

    species=`grep "^[[:space:]]*focal_species[[:space:]]*=" ${config} | cut -d "=" -f 2 | xargs`

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


    echo "[\$species] Extracting ortholog pairs from Newick tree"
    echo "NF internal work directory for [setupAdjustment] process:\n\$processDir\n" > \${logs_folder}/${logs_names["setupAdjustment"]}

    ksrates init ${config} --nextflow >> \${logs_folder}/${logs_names["setupAdjustment"]} 2>&1
    cat rate_adjustment/\$species/ortholog_pairs_\${species}.tsv > \${processDir}/ortholog_pairs_\${species}.tsv

    # If all the ortholog data are present in the databases, already trigger plotOrthologs to plot the ortholog distributions
    if [ -z "`tail -n +2 rate_adjustment/\$species/ortholog_pairs_\${species}.tsv`" ]; then    # if file size NOT greater than 0, so if the file is empty and there aren't unknown pairs
        trigger_plot_orthologs=true
    else
        echo "Ortholog data in DB and / or TSV files not present for one or more species pair required for rate-adjustment"
        echo "Ortholog Ks analysis is scheduled"
        echo "The species pair list can be found in [rate_adjustment/\$species/ortholog_pairs_\${species}.tsv]"
    fi

    echo "[\$species] log can be found in: \${logs_folder}/${logs_names["setupAdjustment"]}"
    
    cd \$processDir
    """
}

/*
 * Log output to stdout/console.
 * (output will be logged only after process setupAdjustment terminates)
 */
if (LOG) {
    outsetupAdjustment
        .splitText() {
            if ( it.endsWith("\n") ) {
                "[setupAdjustment] " + it.replace('\n', '')
            } else {
                "[setupAdjustment] $it"
            }
        }
        .subscribe { log.info it }
}



/*
 * Process that checks if the .ks.tsv file containing the paralog Ks values
 * for the focal species is already present. If not, triggers 
 * wgdParalogs for their estimate.
 */ 
process setParalogAnalysis {

    executor 'local'
    maxForks 1

    input:
        val species from species_channel
        val logs_folder from logs_folder_channel
        file config from configfile

    output:
        stdout outsetParalogAnalysis
        env trigger_wgdPara into trigger_wgdPara_channel
        env trigger_doRateCorrection_from_para into trigger_doRateCorrection_from_setParalog_channel

    script:
    if (LOG) {
        log.info "Process setParalogAnalysis (${task.index}) has " +\
                 "${task.cpus} ${ task.cpus == 1 ? 'CPU' : 'CPUs' } available."
    }
    """
    echo "[$species] Organizing paralog wgd runs"

    trigger_doRateCorrection_from_para=false
    trigger_wgdPara=false

    paranome_status="not_required"
    colinearity_status="not_required"

    paranome=`grep "^[[:space:]]*paranome[[:space:]]*=" ${config} | cut -d "=" -f 2 | xargs | tr '[:upper:]' '[:lower:]'`
    colinearity=`grep "^[[:space:]]*collinearity[[:space:]]*=" ${config} | cut -d "=" -f 2 | xargs | tr '[:upper:]' '[:lower:]'`

    processDir=\$PWD
    cd $PWD

    echo "NF internal work directory for [setParalogAnalysis (${task.index})] process:\n\$processDir\n" > $logs_folder/${logs_names["setParalogAnalysis"]}

    # Triggering wgdParalog process only if ".ks.tsv" (and ".ks_anchors.tsv") files are missing

    if [ \${paranome} = "yes" ]; then
        if [ ! -f paralog_distributions/wgd_${species}/${species}.ks.tsv ]; then
            echo "[$species] Paralog TSV file not found [${species}.ks.tsv]; wgd pipeline will be started"
            echo "[$species] Paralog TSV file not found [${species}.ks.tsv]" >> $logs_folder/${logs_names["setParalogAnalysis"]}
            echo "[$species] Whole-paranome wgd pipeline will be started\n" >> $logs_folder/${logs_names["setParalogAnalysis"]}
            paranome_status="todo"
        else
            paranome_status="already_done"
        fi
    fi
    if [ \${colinearity} = "yes" ]; then
        if [ ! -f paralog_distributions/wgd_${species}/${species}.ks_anchors.tsv ]; then
            echo "[$species] Anchor pairs TSV file not found [${species}.ks_anchors.tsv]; wgd pipeline will be started"
            echo "[$species] Anchor pairs TSV file not found [${species}.ks_anchors.tsv]" >> $logs_folder/${logs_names["setParalogAnalysis"]}
            echo "[$species] Co-linearity wgd pipeline will be started\n" >> $logs_folder/${logs_names["setParalogAnalysis"]}
            colinearity_status="todo"
        else
            colinearity_status="already_done"
        fi
    fi

    if [ \$paranome_status = "todo" ] || [ \$colinearity_status = "todo" ]; then
        # Trigger wgdParalog process to get the missing tsv file(s)
        trigger_wgdPara=true
    else
        echo "[$species] Paralog TSV file(s) already present; skipping paralog wgd pipeline"
        echo "[$species] Paralog TSV file(s) already present; skipping paralog wgd pipeline\n" >> $logs_folder/${logs_names["setParalogAnalysis"]}
        
        # Trigger doRateAdjustment process to plot (at least) the paralog distribution in the mixed plot
        trigger_doRateCorrection_from_para=true
    fi

    echo "----------------------------------------------------------------\n" >> $logs_folder/${logs_names["setParalogAnalysis"]}
    echo "[$species] log can be found in: $logs_folder/${logs_names["setParalogAnalysis"]}"
    cd \$processDir
    """
}

/*
 * Log output to stdout/console.
 * (output will be logged only after process setParalogAnalysis terminates)
 */
if (LOG) {
    outsetParalogAnalysis
        .splitText() {
            if ( it.endsWith("\n") ) {
                "[setParalogAnalysis] " + it.replace('\n', '')
            } else {
                "[setParalogAnalysis] $it"
            }
        }
        .subscribe { log.info it }
}


/*
 * Process that receives the list of all species pairs missing in database and that 
 * 1) triggers wgdOrthologs for pairs without the ks.tsv file and/or 2) triggers EstimatePeak
 * for pairs without ortholog peak data in database but that do have the .ks.tsv file.
 */
process setOrthologAnalysis {

    executor 'local'
    maxForks 1

    input:
        file ortholog_pairs from check_ortholog_pairs_channel
        val logs_folder from logs_folder_channel
        file config from configfile

    output:
        stdout outsetOrthologAnalysis
        file "tmp_species_pairs_for_wgdOrtholog.txt" optional true into species_pairs_for_wgd_Orthologs_channel
        file "tmp_species_pairs_for_estimatePeak.txt" optional true into file_for_estimatePeak_channel
        env estimatePeak_not_needed into trigger_plotOrthologs_together_with_wgdOrtholog_channel
        env wgdOrtholog_not_needed into trigger_plotOrthologs_together_with_estimatePeak_channel

    script:
    if (LOG) {
        log.info "Process setOrthologAnalysis (${task.index}) has " +\
                 "${task.cpus} ${ task.cpus == 1 ? 'CPU' : 'CPUs' } available."
    }
    """
    processDir=\$PWD
    cd $PWD

    echo "NF internal work directory for [setOrthologAnalysis (${task.index})] process:\n\$processDir\n" >> $logs_folder/${logs_names["setOrthologAnalysis"]}

    wgdOrtholog_not_needed=false
    estimatePeak_not_needed=false

    # If the file for the species pair is empty, then setOrthologAnalysis skips the ortholog Ks estimate and peak estimate steps. 
    if [ -z "`tail -n +2 \${processDir}/$ortholog_pairs`" ]; then
        echo "No species pairs are listed for wgd ortholog pipeline or ortholog peak estimate."
        echo "No species pairs are listed for wgd ortholog pipeline or ortholog peak estimate." >> $logs_folder/${logs_names["setOrthologAnalysis"]}
    fi 

    while read -r species1 species2 || [ -n "\$species1" ]; do
        if [ ! -f ortholog_distributions/wgd_\${species1}_\${species2}/\${species1}_\${species2}.ks.tsv ]; then
            echo "\$species1\t\$species2" >> \${processDir}/tmp_species_pairs_for_wgdOrtholog.txt
            echo "[\$species1 – \$species2] Ortholog TSV file not present [\${species1}_\${species2}.ks.tsv]" >> $logs_folder/${logs_names["setOrthologAnalysis"]}
            echo "[\$species1 – \$species2] Ortholog wgd pipeline will be started" >> $logs_folder/${logs_names["setOrthologAnalysis"]}
        else
            echo "[\$species1 – \$species2] Ortholog TSV file already present; skipping ortholog wgd pipeline" >> $logs_folder/${logs_names["setOrthologAnalysis"]}
            echo "[\$species1 – \$species2] Update of ortholog database(s) will be started." >> $logs_folder/${logs_names["setOrthologAnalysis"]}

            if [ ! -f \${processDir}/tmp_species_pairs_for_estimatePeak.txt ]; then   # if the file does not exist yet, add headers
                echo "Species1\tSpecies2" >> \${processDir}/tmp_species_pairs_for_estimatePeak.txt
            fi
            echo "\$species1\t\$species2" >> \${processDir}/tmp_species_pairs_for_estimatePeak.txt   # add a species pair
        fi
        echo " \n----------------------------------------------------------------\n" >> $logs_folder/${logs_names["setOrthologAnalysis"]}
        echo "[\$species1 – \$species2] log can be found in: $logs_folder/${logs_names["setOrthologAnalysis"]}"
    done < <(tail -n +2 \${processDir}/$ortholog_pairs) # skipping the headers

    if [ ! -s \${processDir}/tmp_species_pairs_for_wgdOrtholog.txt ]; then    # if the file for wgdOrthologs is empty (size not > 0)
        wgdOrtholog_not_needed=true
    else
        if [ ! -d ortholog_distributions ]; then
            mkdir ortholog_distributions
        fi
    fi

    if [ -z "`tail -n +2 \${processDir}/tmp_species_pairs_for_estimatePeak.txt`" ]; then    # if the file content for estimatePeaks (skipping headers) is empty (empty string check)
        estimatePeak_not_needed=true
    fi

    cd \$processDir
    """
}

/*
 * Log output to stdout/console.
 * (output will be logged only after process setOrthologAnalysis terminates)
 */
if (LOG) {
    outsetOrthologAnalysis
        .splitText() { 
            if ( it.endsWith("\n") ) {
                "[setOrthologAnalysis] " + it.replace('\n', '') 
            } else {
                "[setOrthologAnalysis] $it" 
            }
        }
        .subscribe { log.info it }
}



/*
 * Process that estimates the ortholog Ks distribution peak for a species pair 
 * that has already its .ks.tsv file, but that for some reasons 
 * is not present anymore in the ortholog database(s) (i.e. has been deleted).
 */
process estimatePeak {

    input:
        file species_pairs_for_peak from file_for_estimatePeak_channel
        val logs_folder from logs_folder_channel
        file config from configfile
        
    output:
        stdout outestimatePeak
        val true into trigger_doRateCorrection_from_estimatePeak_channel
        val true into trigger_plotOrtholog_from_estimatePeak_channel
        val true into trigger_plotOrtholog_from_estimatePeak_together_with_wgdOrthologs_channel

    script:
    if (LOG) {
        log.info "Process estimatePeak (${task.index}) has " +\
                 "${task.cpus} ${ task.cpus == 1 ? 'CPU' : 'CPUs' } available."
    }
    """
    processDir=\$PWD
    cd $PWD

    echo "NF internal work directory for [estimatePeak (${task.index})] process:\n\$processDir\n" >> $logs_folder/${logs_names["estimatePeak"]}
    echo "Updating ortholog peak database" >> $logs_folder/${logs_names["estimatePeak"]}

    ksrates orthologs-analysis ${config} --ortholog-pairs=\$processDir/$species_pairs_for_peak >> $logs_folder/${logs_names["estimatePeak"]} 2>&1
    """
}

/*
 * Log output to stdout/console.
 * (output will be logged only after process estimatePeak terminates)
 */
if (LOG) {
    outestimatePeak
        .splitText() { 
            if ( it.endsWith("\n") ) {
                "[estimatePeak] " + it.replace('\n', '') 
            } else {
                "[estimatePeak] $it" 
            }
        }
        .subscribe { log.info it }
}



/*
 * Process that estimates the paranome and/or anchor pair Ks values
 * for the focal species.
 */
process wgdParalogs {

    input:
        val species from species_channel
        val trigger_wgdPara from trigger_wgdPara_channel
        val logs_folder from logs_folder_channel
        file config from configfile

    output:
        stdout outParalogs
        val true into trigger_doRateCorrection_from_wgdParalog_channel

    when:
        trigger_wgdPara == "true"

    script:
    if (LOG) {
        log.info "Process wgdParalogs (${task.index}) will use species $species and has " +\
                 "${task.cpus} ${ task.cpus == 1 ? 'CPU' : 'CPUs' } available."
    }
    """
    echo "[$species] `date`"
    echo "[$species] Starting paralog wgd analysis"

    processDir=\$PWD
    cd $PWD
    echo "NF internal work directory for [wgdParalogs] process:\n\$processDir\n" >> $logs_folder/${logs_names["wgdParalogs"]}

    echo "Using ${task.cpus} thread(s)\n">> $logs_folder/${logs_names["wgdParalogs"]}
    echo "[$species] Using ${task.cpus} thread(s)"

    ksrates paralogs-ks ${config} --n-threads=${task.cpus} >> $logs_folder/${logs_names["wgdParalogs"]} 2>&1

    RET_CODE=\$?
    echo "[$species] Done [\${RET_CODE}]"
    echo "[$species] wgd log can be found in: $logs_folder/${logs_names["wgdParalogs"]}"
    echo "[$species] `date`"

    cd \$processDir
    """
}

/*
 * Log output to stdout/console.
 * (output will be logged only after process wgdParalogs terminates)
 */
if (LOG) {
    outParalogs
        .splitText() { 
            if ( it.endsWith("\n") ) {
                "[wgdParalogs] " + it.replace('\n', '') 
            } else {
                "[wgdParalogs] $it" 
            }
        }
        .subscribe { log.info it }
}



/*
 * Process that estimates the ortholog Ks values for a species pair and
 * estimates also the ortholog Ks distribution peak.
 */
process wgdOrthologs {

    input:
        val species from species_channel
        tuple species1, species2 from species_pairs_for_wgd_Orthologs_channel.splitCsv(sep:'\t')
        val logs_folder from logs_folder_channel
        file config from configfile

    output:
        stdout outOrthologs
        val true into trigger_doRateCorrection_from_wgdOrtholog_channel
        val true into trigger_plotOrtholog_from_wgdOrtholog_channel
        val true into trigger_plotOrtholog_from_wgdOrtholog_together_with_estimatePeak_channel

    script:
    if (LOG) {
        log.info "Process wgdOrthologs (${task.index}) will use species pair $species1 – $species2 and has " +\
                 "${task.cpus} ${ task.cpus == 1 ? 'CPU' : 'CPUs' } available."
    }
    """
    echo "[$species1 – $species2] `date`"
    echo "[$species1 – $species2] Starting ortholog wgd analysis"

    processDir=\${PWD}
    cd $PWD
    echo "NF internal work directory for [wgdOrthologs (${task.index})] process:\n\$processDir\n" > $logs_folder/${logs_names["wgdOrthologs"]}${species1}_${species2}.log

    echo "Using ${task.cpus} thread(s)\n">> $logs_folder/${logs_names["wgdOrthologs"]}${species1}_${species2}.log
    echo "[$species1 – $species2] Using ${task.cpus} thread(s)"

    ksrates orthologs-ks ${config} $species1 $species2 --n-threads=${task.cpus} >> $logs_folder/${logs_names["wgdOrthologs"]}${species1}_${species2}.log 2>&1
    RET_CODE=\$?
    echo "[$species1 – $species2] wgd done [\${RET_CODE}]"
    
    echo "[$species1 – $species2] Computing ortholog peak and error"
    echo "Species1\tSpecies2\n$species1\t$species2" > \${processDir}/tmp_${species1}_${species2}.txt
    ksrates orthologs-analysis ${config} --ortholog-pairs=\${processDir}/tmp_${species1}_${species2}.txt >> $logs_folder/${logs_names["wgdOrthologs"]}${species1}_${species2}.log 2>&1

    RET_CODE=\$?
    echo "[$species1 – $species2] Compute peak done [\${RET_CODE}]"
    echo "[$species1 – $species2] wgd log and compute peak log can be found in: $logs_folder/${logs_names["wgdOrthologs"]}${species1}_${species2}.log"
    echo "[$species1 – $species2] `date`"

    cd \$processDir
    """
}

/*
 * Log output to stdout/console.
 * (output will be logged only after process wgdOrthologs terminates)
 */
if (LOG) {
    outOrthologs
        .splitText() { 
            if ( it.endsWith("\n") ) {
                "[wgdOrthologs] " + it.replace('\n', '') 
            } else {
                "[wgdOrthologs] $it" 
            }
        }
        .subscribe { log.info it }
}



/*
 * Process that generates figures collecting all the ortholog distributions 
 * used for the adjustment and that highlights their peaks.
 */
process plotOrthologDistrib {

    maxForks 1

    input:
        val species from species_channel
        val trigger from trigger_plotOrtholog_from_setupCorrection_channel.mix(trigger_plotOrtholog_from_estimatePeak_together_with_wgdOrthologs_channel.merge(trigger_plotOrtholog_from_wgdOrtholog_together_with_estimatePeak_channel.collect()), trigger_plotOrthologs_together_with_wgdOrtholog_channel.merge(trigger_plotOrtholog_from_wgdOrtholog_channel.collect()), trigger_plotOrthologs_together_with_estimatePeak_channel.merge(trigger_plotOrtholog_from_estimatePeak_channel))
        val logs_folder from logs_folder_channel
        file config from configfile

    output:
        stdout outplotOrthologDistrib

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
         * An accepted trigger comes from setupAdjustment if all ortholog data are already present and ready to be plotted; from estimatePeak together with wgdOrthologs 
         * when they have both finished with all the peaks; from setOrthologs and wgdOrthologs after this latter has done with all ortholog distributions;
         * from setOrtholog and estimatePeak when all missing peaks are computed.
         */

    script:
    if (LOG) {
        log.info "Process plotOrthologDistrib (${task.index}) has " +\
                 "${task.cpus} ${ task.cpus == 1 ? 'CPU' : 'CPUs' } available."
    }
    """
    echo "[$species] `date`"
    echo "[$species] Plotting ortholog distributions"

    processDir=\$PWD
    cd $PWD
    echo "NF internal work directory for [plotOrthologDistrib] process:\n\$processDir\n" > $logs_folder/${logs_names["plotOrthologDistrib"]}

    ksrates plot-orthologs ${config} >> $logs_folder/${logs_names["plotOrthologDistrib"]} 2>&1

    RET_CODE=\$?
    echo "[$species] Done [\${RET_CODE}]"
    echo "[$species] log can be found in: $logs_folder/${logs_names["plotOrthologDistrib"]}"
    echo "[$species] `date`"
    """
}

/*
 * Log output to stdout/console.
 * (output will be logged only after process plotOrthologDistrib terminates)
 */
if (LOG) {
    outplotOrthologDistrib
        .splitText() { 
            if ( it.endsWith("\n") ) {
                "[plotOrthologDistrib] " + it.replace('\n', '') 
            } else {
                "[plotOrthologDistrib] $it" 
            }
        }
        .subscribe { log.info it }
}



/*
 * Process that performs the rate-adjustment and generates the
 * mixed distribution plots.
 */
process doRateAdjustment {

    executor 'local'
    maxForks 1

    input:
        val species from species_channel
        val trigger from trigger_doRateCorrection_from_setParalog_channel.mix(trigger_doRateCorrection_from_estimatePeak_channel, trigger_doRateCorrection_from_wgdOrtholog_channel, trigger_doRateCorrection_from_wgdParalog_channel)
        val logs_folder from logs_folder_channel
        file config from configfile

    output:
        stdout outRateAdjustmentAnalysis
        val true into trigger_peakCalling_from_doRateCorrection_channel
        val true into trigger_drawTree_from_doRateCorrection_channel

    when:
        trigger == "true" || trigger == true
        /*
         * String "true" comes from setParalog process, which returns either a string "true" or a string "false" (it uses "env" variable, not a boolean);
         * Boolean true come from the other triggers, which return a boolean "val" true (trigger_doRateCorrection_from_estimatePeak_channel, trigger_doRateCorrection_from_wgdOrtholog_channel, trigger_doRateCorrection_from_wgdParalog_channel)
         *
         * An accepted trigger can come from setParalog (if paralog data is already present it returns string "true", so that at least the paralog distribution will be plotted),
         *     or from estimatePeak (updates the databases and sends a boolean true) or from wgdOrthologs (when finishes a single run returns boolean true) or finally
         *     from wgdParalogs (when finishes returns boolean true).
         */

    script:
    if (LOG) {
        log.info "Process doRateAdjustment (${task.index}) will use species $species and has " +\
                 "${task.cpus} ${ task.cpus == 1 ? 'CPU' : 'CPUs' } available."
    }
    """
    echo "[$species] `date`"
    echo "[$species] Starting rate-adjustment analysis"

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

    if [ \${missing_paranome} = "true" ] && [ \${missing_anchorpairs} = "true" ]; then
        echo "[$species] Whole-paranome and anchor pair data do not yet exist -> skipping rate-adjustment"
    fi
    if [ \${missing_paranome} = "true" ] && [ \${missing_anchorpairs} = "false" ]; then
        echo "[$species] Whole-paranome data [${species}.ks.tsv] does not yet exist -> skipping rate-adjustment"
    fi
    if [ \${missing_paranome} = "false" ] && [ \${missing_anchorpairs} = "true" ]; then
        echo "[$species] Anchor pair data [${species}.ks_anchors.tsv] does not yet exist -> skipping rate-adjustment"
    fi
    if [ \${missing_paranome} = "true" ] || [ \${missing_anchorpairs} = "true" ]; then
        echo "[$species] Done"
        echo "[$species] `date`"
        exit 0
    fi

    processDir=\$PWD
    cd $PWD
    
    echo "NF internal work directory for [doRateAdjustment (${task.index})] process:\n\$processDir\n" >> $logs_folder/${logs_names["doRateAdjustment"]}

    echo "[$species] Performing rate-adjustment"
    ksrates orthologs-adjustment ${config} >> $logs_folder/${logs_names["doRateAdjustment"]} 2>&1
    echo "\n" >> $logs_folder/${logs_names["doRateAdjustment"]}

    echo "[$species] Plotting mixed distributions"
    ksrates plot-paralogs ${config} >> $logs_folder/${logs_names["doRateAdjustment"]} 2>&1
    echo "\n" >> $logs_folder/${logs_names["doRateAdjustment"]}

    echo " \n-----------------------------------------------------------\n" >> $logs_folder/${logs_names["doRateAdjustment"]}

    RET_CODE=\$?
    echo "[$species] Done [\${RET_CODE}]"
    echo "[$species] log can be found in: $logs_folder/${logs_names["doRateAdjustment"]}"
    echo "[$species] `date`"
    """
}

/*
 * Log output to stdout/console.
 * (output will be logged only after process doRateAdjustment terminates)
 */
if (LOG) {
    outRateAdjustmentAnalysis
        .splitText() { 
            if ( it.endsWith("\n") ) {
                "[doRateAdjustment] " + it.replace('\n', '') 
            } else {
                "[doRateAdjustment] $it" 
            }
        }
        .subscribe { log.info it }
}



/*
 * Process that performs the WGM peak calling and generates the 
 * related figures (anchor pair Ks clustering and mixture model).
 */
process paralogsAnalyses {

    input:
        val species from species_channel
        val logs_folder from logs_folder_channel
        file config from configfile
        val trigger from trigger_peakCalling_from_doRateCorrection_channel.collect()

    output:
        stdout outPeakCalling

    script:
    /*
     * It is performed only once as last process after all the doRateAdjustment's,
     * when the adjustment table is certainly complete.
     */
    """
    echo "[$species] `date`"
    echo "[$species] Performing WGM peak calling"

    processDir=\$PWD
    cd $PWD
    echo "NF internal work directory for [paralogsAnalyses (${task.index})] process:\n\$processDir\n" >> $logs_folder/${logs_names["paralogsAnalyses"]}

    # If colinearity is on, perform peak calling only on anchor distribution with cluster-anchors
    # Else, perform peak calling only on paranome with mixture model(s)
    # Additional methods can be required through the expert configuration file

    ksrates paralogs-analyses ${config} >> $logs_folder/${logs_names["paralogsAnalyses"]} 2>&1
    echo "\n" >> $logs_folder/${logs_names["paralogsAnalyses"]}

    RET_CODE=\$?
    echo "[$species] Done [\${RET_CODE}]"
    echo "[$species] log can be found in: $logs_folder/${logs_names["paralogsAnalyses"]}"
    echo "[$species] `date`"
    """
}

/*
 * Log output to stdout/console.
 * (output will be logged only after process paralogsAnalyses terminates)
 */
if (LOG) {
    outPeakCalling
        .splitText() { 
            if ( it.endsWith("\n") ) {
                "[paralogsAnalyses] " + it.replace('\n', '') 
            } else {
                "[paralogsAnalyses] $it" 
            }
        }
        .subscribe { log.info it }
}



/*
 * Process that plots the input tree with branch length equal to 
 * the Ks distances.
 */
process drawTree {

    executor 'local'
    maxForks 1

    input:
        val species from species_channel
        val logs_folder from logs_folder_channel
        file config from configfile
        val trigger from trigger_drawTree_from_doRateCorrection_channel.collect()

    output:
        stdout outDrawTree

    script:
    /*
     * It is performed only once as last process after all the doRateAdjustment's,
     * when the adjustment table is certainly complete.
     */
    """
    echo "[$species] `date`"
    echo "[$species] Plotting tree with branch length equal to Ks distances"

    processDir=\$PWD
    cd $PWD
    echo "NF internal work directory for [drawTree (${task.index})] process:\n\$processDir\n" >> $logs_folder/${logs_names["drawTree"]}

    ksrates plot-tree ${config} --nextflow >> $logs_folder/${logs_names["drawTree"]} 2>&1

    RET_CODE=\$?
    echo "[$species] Done [\${RET_CODE}]"
    echo "[$species] log can be found in: $logs_folder/${logs_names["drawTree"]}"
    echo "[$species] `date`"
    """
}

/*
 * Log output to stdout/console.
 * (output will be logged only after process drawTree terminates)
 */
if (LOG) {
    outDrawTree
        .splitText() { 
            if ( it.endsWith("\n") ) {
                "[drawTree] " + it.replace('\n', '') 
            } else {
                "[drawTree] $it" 
            }
        }
        .subscribe { log.info it }
}



/*
 * On completion of workflow clean up any temporary files left behind.
 */
workflow.onComplete {
    if (LOG) {
        log.info ""
        log.info("Logs folder: logs_${workflow.sessionId.toString().substring(0,8)}")

        if ( params.preserve == false ) {
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

            log.info "Done."
            log.info ""
        }

        log.info "Pipeline completed at: $workflow.complete taking $workflow.duration"
        log.info "Execution status: ${ workflow.success ? 'OK' : 'failed' }"
        log.info ""
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
        log_filename = "rate_adjustment/${focal_species}/logs_${workflow.sessionId.toString().substring(0,8)}/${logs_names[process]}"

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
            }

            // In any case, point to the complete output of the stopped process log file
            error_box += "\nMore details can be found in the error report above or in the following log file:\n" + \
                         "${log_filename}\n"
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
                error_box += "\nMore details may be found in the error report above or in ./.nextflow.log.\n"
            }
            else {
                error_box += "\nMore details may be found in the error report above, in ./.nextflow.log, or in one of\n" + \
                             "the log files named rate_adjustment/${focal_species}/" + \
                             "logs_${workflow.sessionId.toString().substring(0,8)}/${logs_names["wgdOrthologs"]} (if any exists).\n"
            }
        }

        // Separator to highlight the end of the error box
        error_box += "${'=' * headline.length()}"

        log.error "${error_box}"
    }
}
