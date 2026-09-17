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
KEY_DIR=./paralog_ks_sqld_keys
ADDRESS_FILE=./paralog_ks_server_address.txt
PORT=8080
# Must match this ksrates project's "ks_list_paralog_database_path" config field.
# ---------------------------------------------------------------------------------------

SQLD_BIN="$SQLD_INSTALL_DIR/bin/sqld"

# Install sqld once, if not already present. Only works if this job's compute node has outbound
# internet access; many clusters restrict that to login nodes only. If this fails, run the same
# download/extract steps by hand on a node that does have internet access, pointed at the same
# SQLD_INSTALL_DIR, before submitting this job.
#
# Fetches the release tarball directly rather than going through the official
# libsql-server-installer.sh: that installer refuses to use the glibc build on any system with
# glibc older than 2.35 and falls back to requesting a "musl-static" asset instead — but no such
# asset is actually published for Linux, so on a node with (for example) glibc 2.34 the installer
# always fails outright. The glibc build is normally compatible one or two point releases below
# its build target, so it's still worth trying directly here instead of trusting that gate.
if [ ! -x "$SQLD_BIN" ]; then
    echo "sqld binary not found at $SQLD_BIN, downloading $SQLD_VERSION..."
    mkdir -p "$SQLD_INSTALL_DIR/bin"
    case "$(uname -m)" in
        x86_64)  SQLD_ASSET="libsql-server-x86_64-unknown-linux-gnu.tar.xz" ;;
        aarch64) SQLD_ASSET="libsql-server-aarch64-unknown-linux-gnu.tar.xz" ;;
        *)
            echo "ERROR: unsupported architecture $(uname -m); no sqld release for it."
            exit 1
            ;;
    esac
    SQLD_TARBALL=$(mktemp)
    if ! curl --proto '=https' --tlsv1.2 -LsSf -o "$SQLD_TARBALL" \
            "https://github.com/tursodatabase/libsql/releases/download/${SQLD_VERSION}/${SQLD_ASSET}"; then
        echo "ERROR: failed to download $SQLD_ASSET. This compute node likely has no internet"
        echo "access. Download it manually from a node that does, extract it, and place the"
        echo "'sqld' binary at $SQLD_BIN before resubmitting this job."
        rm -f "$SQLD_TARBALL"
        exit 1
    fi
    SQLD_EXTRACT_DIR=$(mktemp -d)
    if ! tar -xJf "$SQLD_TARBALL" -C "$SQLD_EXTRACT_DIR"; then
        echo "ERROR: downloaded $SQLD_ASSET but could not extract it (this node's 'tar' may be"
        echo "missing xz support). Extract $SQLD_TARBALL manually and place the 'sqld' binary"
        echo "at $SQLD_BIN before resubmitting this job."
        rm -f "$SQLD_TARBALL"
        rm -rf "$SQLD_EXTRACT_DIR"
        exit 1
    fi
    rm -f "$SQLD_TARBALL"
    find "$SQLD_EXTRACT_DIR" -type f -name sqld -exec cp {} "$SQLD_BIN" \;
    rm -rf "$SQLD_EXTRACT_DIR"
    chmod +x "$SQLD_BIN" 2>/dev/null
fi

if [ ! -x "$SQLD_BIN" ]; then
    echo "ERROR: sqld binary still not found/executable at $SQLD_BIN after attempting download."
    echo "Download $SQLD_ASSET manually from a node with internet access, extract it, and place"
    echo "its 'sqld' binary at $SQLD_BIN before resubmitting this job."
    exit 1
fi

mkdir -p "$SQLD_DATA_DIR"

# Authentication: every client on the network can reach this port, so the server requires a
# JWT on every request. Kept deliberately simple: one shared, non-expiring,
# full-access token generated once and reused by every client job. This is safe because the
# key pair (and thus the token) lives in KEY_DIR, next to the data directory, not in this
# script or in version control, and is never printed to the job log.
#
# The private key is generated once and kept persistent across server restarts so old tokens
# (already handed out in old address files) keep working; the public key is what sqld itself
# uses to verify tokens (--auth-jwt-key-file).
mkdir -p "$KEY_DIR"
PRIVATE_KEY="$KEY_DIR/jwt_private.pem"
PUBLIC_KEY="$KEY_DIR/jwt_public.pem"
if [ ! -f "$PRIVATE_KEY" ]; then
    echo "No JWT signing key found at $PRIVATE_KEY, generating one (one-time setup)..."
    ( umask 077 && openssl genpkey -algorithm ed25519 -out "$PRIVATE_KEY" )
    openssl pkey -in "$PRIVATE_KEY" -pubout -out "$PUBLIC_KEY"
fi

# base64url = base64 with '+/' -> '-_' and no padding, per the JWT spec (RFC 7515 Appendix C).
b64url() { base64 | tr -d '=\n' | tr '+/' '-_'; }

# EdDSA header, empty payload (sqld treats an empty payload as valid and grants full
# read/write access, with no expiration required, so this token never needs to be rotated).
HEADER_B64=$(printf '%s' '{"alg":"EdDSA","typ":"JWT"}' | b64url)
PAYLOAD_B64=$(printf '%s' '{}' | b64url)
SIGNING_INPUT="${HEADER_B64}.${PAYLOAD_B64}"
# openssl's Ed25519 signer needs a seekable input (it computes the message length up front),
# so the signing input has to go through a temp file rather than a pipe.
SIGNING_INPUT_FILE=$(mktemp)
printf '%s' "$SIGNING_INPUT" > "$SIGNING_INPUT_FILE"
SIGNATURE_B64=$(openssl pkeyutl -sign -inkey "$PRIVATE_KEY" -rawin -in "$SIGNING_INPUT_FILE" | b64url)
rm -f "$SIGNING_INPUT_FILE"
AUTH_TOKEN="${SIGNING_INPUT}.${SIGNATURE_B64}"

# Written once at startup: every ksrates client job reads this file to find the server and,
# on the second line, the token needed to authenticate to it.
{
    echo "$(hostname):$PORT"
    echo "$AUTH_TOKEN"
} > "$ADDRESS_FILE"
chmod 600 "$ADDRESS_FILE"

echo "Server node: $(hostname)"
echo "Port: $PORT"
echo "Data directory: $(realpath "$SQLD_DATA_DIR")"
echo "Address file: $(realpath "$ADDRESS_FILE")"
echo "Starting sqld..."

exec "$SQLD_BIN" --http-listen-addr=0.0.0.0:$PORT --db-path="$SQLD_DATA_DIR" --auth-jwt-key-file="$PUBLIC_KEY"
