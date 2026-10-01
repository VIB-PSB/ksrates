#!/usr/bin/env nextflow

/*
 * Launches the shared paralog Ks database server (sqld) via "ksrates launch-paralog-ks-server".
 * This never finishes on success - the underlying process execs into sqld and runs indefinitely -
 * so submit this once, as its own long-running job, separately from any analysis pipeline run
 * (e.g. via main.nf), and leave it running. See docs/configuration.rst, section
 * "Paralog Ks database server", for the full setup instructions.
 *
 * Uses the same container/profile settings (-profile docker/apptainer/singularity, -c nextflow.config)
 * as main.nf - only the entry script and params differ. If nextflow.config sets executor.name to a
 * cluster scheduler (e.g. 'slurm'), override it back to 'local' with -process.executor, so this
 * single always-running process isn't submitted as its own nested job subject to the cluster's
 * default walltime:
 *
 *   nextflow run VIB-PSB/ksrates -main-script setup_database_server.nf -profile singularity \
 *       -c nextflow.config -process.executor=local \
 *       --location /path/to/central/dir
 *
 * Submit the command above as (or from within) your own long-running/walltime-unlimited job.
 */

// Central directory to hold the server's data/keys directories and address file (required)
params.location = false

// Port sqld listens on
params.port = 8080

// Filename (not path) of the address file written under params.location
params.address_filename = "paralog_ks_server_address.txt"

if (!params.location) {
    log.error "Missing required parameter --location (central directory for the server's data/keys/address file)."
    exit 1
}

process launchParalogKsServer {

    input:
        val location
        val port
        val address_filename

    script:
    """
    ksrates launch-paralog-ks-server --location ${location} --port ${port} --address-filename ${address_filename}
    """
}

workflow {
    launchParalogKsServer(params.location, params.port, params.address_filename)
}
