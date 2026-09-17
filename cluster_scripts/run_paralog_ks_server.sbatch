#!/bin/bash
# Launches the sqld server that hosts the shared paralog Ks database (see docs/configuration.rst,
# field "ks_list_paralog_database_path"). Submit this ONCE; with a high (or, if possible, unlimited) job walltime it can then
# run for a long time, and every ksrates analysis job on the cluster becomes a network client of it
# instead of opening a local database file directly.
#
# Copy this file and edit the variables below for your own cluster before submitting
# (sbatch run_paralog_ks_server.sbatch). If your cluster enforces a walltime limit, add a
# "#SBATCH --time=..." line here matching your longest available partition, and be prepared to
# resubmit this job (and update ADDRESS_FILE) if/when it gets killed on walltime.

#SBATCH --job-name=para_server
#SBATCH --output=paralog_ks_server.log
#SBATCH --nodes=1
#SBATCH --ntasks=1

# --- Edit these for your own setup ---------------------------------------------------
SQLD_VERSION=libsql-server-v0.24.32
SQLD_INSTALL_DIR=./sqld_bin
SQLD_DATA_DIR=./paralog_ks_sqld_data
ADDRESS_FILE=./paralog_ks_server_address.txt
PORT=8080
# Must match this ksrates project's "ks_list_paralog_database_path" config field.
# ---------------------------------------------------------------------------------------

SQLD_BIN="$SQLD_INSTALL_DIR/bin/sqld"

# Install sqld once, if not already present. Only works if this job's compute node has outbound
# internet access; many clusters restrict that to login nodes only. If this fails, run this same
# install command by hand on a node that does have internet access, pointed at the same
# SQLD_INSTALL_DIR, before submitting this job.
if [ ! -x "$SQLD_BIN" ]; then
    echo "sqld binary not found at $SQLD_BIN, attempting to install $SQLD_VERSION..."
    mkdir -p "$SQLD_INSTALL_DIR"
    CARGO_DIST_FORCE_INSTALL_DIR="$(realpath "$SQLD_INSTALL_DIR")" \
        curl --proto '=https' --tlsv1.2 -LsSf \
        "https://github.com/tursodatabase/libsql/releases/download/${SQLD_VERSION}/libsql-server-installer.sh" | sh
fi

if [ ! -x "$SQLD_BIN" ]; then
    echo "ERROR: sqld binary still not found/executable at $SQLD_BIN after attempting install."
    echo "This compute node likely has no internet access. Install sqld manually (see comment"
    echo "above) from a node that does, then resubmit this job."
    exit 1
fi

mkdir -p "$SQLD_DATA_DIR"

# Written once at startup: every ksrates client job reads this file to find the server.
echo "$(hostname):$PORT" > "$ADDRESS_FILE"

echo "Server node: $(hostname)"
echo "Port: $PORT"
echo "Data directory: $(realpath "$SQLD_DATA_DIR")"
echo "Address file: $(realpath "$ADDRESS_FILE")"
echo "Starting sqld..."

exec "$SQLD_BIN" --http-listen-addr=0.0.0.0:$PORT --db-path="$SQLD_DATA_DIR"
