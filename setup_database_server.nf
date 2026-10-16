#!/usr/bin/env nextflow

/*
 * Launches the shared Ks database server (sqld) via "ksrates launch-ks-database-server".
 * Never finishes on success - execs into sqld and runs indefinitely - so submit this once, as its
 * own long-running job separate from any analysis run. See docs/ks_database.rst,
 * "Setting up the database server", for full setup instructions.
 *
 * Override executor.name back to 'local' via -process.executor if nextflow.config sets it to a
 * cluster scheduler, so this always-running process isn't itself submitted as a nested job
 * subject to the cluster's default walltime:
 *
 *   nextflow run VIB-PSB/ksrates -main-script setup_database_server.nf -profile singularity \
 *       -c nextflow.config -process.executor=local --location /path/to/central/dir
 */

// Central directory to hold the server's data/keys directories and address file (required)
params.location = false

// Port sqld listens on
params.port = 8080

// Filename (not path) of the address file written under params.location
params.address_filename = "ks_database_server_address.txt"

if (!params.location) {
    log.error "Missing required parameter --location (central directory for the server's data/keys/address file)."
    exit 1
}

process launchKsDatabaseServer {

    input:
        val location
        val port
        val address_filename

    script:
    """
    ksrates launch-ks-database-server --location ${location} --port ${port} --address-filename ${address_filename}
    """
}

workflow {
    launchKsDatabaseServer(params.location, params.port, params.address_filename)
}
