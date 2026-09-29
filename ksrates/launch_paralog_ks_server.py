import base64
import os
import socket
import sys

from cryptography.hazmat.primitives import serialization
from cryptography.hazmat.primitives.asymmetric.ed25519 import Ed25519PrivateKey

_SQLD_BIN = "/usr/local/bin/sqld"


def _b64url(data):
	"""base64url per the JWT spec (RFC 7515 Appendix C): standard base64 with '+/' -> '-_' and no padding."""
	return base64.urlsafe_b64encode(data).rstrip(b'=').decode('ascii')


def _generate_jwt_keys(key_dir):
	"""
	Load (or generate, if not already present) the Ed25519 signing keypair used to authenticate
	clients. The private key persists across server restarts (kept in key_dir, next to the data
	directory) so tokens already handed out in old address files keep working.

	Uses the `cryptography` package directly rather than shelling out to the openssl CLI: the
	CLI's raw-Ed25519-signing support (pkeyutl -rawin) was only added in OpenSSL 3.0, so shelling
	out is fragile across container images with older OpenSSL builds (seen firsthand: OpenSSL
	1.1.1f, common on Debian-buster-era base images, has no -rawin at all).

	:param key_dir: directory to hold jwt_private.pem/jwt_public.pem
	:return: (private_key, public_key_path) - private_key is an Ed25519PrivateKey object
	"""
	private_key_path = os.path.join(key_dir, "jwt_private.pem")
	public_key_path = os.path.join(key_dir, "jwt_public.pem")

	if os.path.isfile(private_key_path):
		with open(private_key_path, "rb") as f:
			private_key = serialization.load_pem_private_key(f.read(), password=None)
	else:
		print(f"No JWT signing key found at {private_key_path}, generating one (one-time setup)...")
		private_key = Ed25519PrivateKey.generate()
		private_bytes = private_key.private_bytes(
			encoding=serialization.Encoding.PEM,
			format=serialization.PrivateFormat.PKCS8,
			encryption_algorithm=serialization.NoEncryption(),
		)
		with open(private_key_path, "wb") as f:
			f.write(private_bytes)
		os.chmod(private_key_path, 0o600)

	# Always (re)written from the private key, so it self-heals if ever missing/out of sync.
	public_bytes = private_key.public_key().public_bytes(
		encoding=serialization.Encoding.PEM,
		format=serialization.PublicFormat.SubjectPublicKeyInfo,
	)
	with open(public_key_path, "wb") as f:
		f.write(public_bytes)

	return private_key, public_key_path


def _generate_jwt(private_key):
	"""
	Sign a non-expiring, full-access JWT: EdDSA header, empty payload (sqld treats an empty
	payload as valid and grants full read/write access, with no expiration required).

	:param private_key: Ed25519PrivateKey object from _generate_jwt_keys
	:return: the complete JWT string (header.payload.signature)
	"""
	header_b64 = _b64url(b'{"alg":"EdDSA","typ":"JWT"}')
	payload_b64 = _b64url(b'{}')
	signing_input = f"{header_b64}.{payload_b64}".encode('ascii')
	signature_b64 = _b64url(private_key.sign(signing_input))
	return f"{header_b64}.{payload_b64}.{signature_b64}"


def launch_server(location, port=8080, address_filename="paralog_ks_server_address.txt",
				   max_response_size="200MB", max_total_response_size="500MB"):
	"""
	Launch the sqld server that hosts the shared paralog Ks database, in the foreground - this
	call never returns (it execs into sqld), so it's meant to be run once as a long-running job
	(e.g. via setup_database_server.nf) and left running indefinitely. Every ksrates analysis then becomes
	a network client of it instead of opening a local database file directly.

	:param location: central directory to hold the data/keys directories and address file; every
	                  dataset's ks_list_paralog_database_path should point at
	                  <location>/<address_filename>
	:param port: port sqld listens on (default: 8080)
	:param address_filename: filename (not path) of the address file written under location
	                          (default: "paralog_ks_server_address.txt"); customizable so multiple
	                          independent servers can share the same location without colliding
	:param max_response_size: passed through to sqld's --max-response-size (default: "200MB")
	:param max_total_response_size: passed through to sqld's --max-total-response-size (default: "500MB")
	"""
	location = os.path.abspath(location)
	data_dir = os.path.join(location, "paralog_ks_sqld_data")
	key_dir = os.path.join(location, "paralog_ks_sqld_keys")
	address_file = os.path.join(location, address_filename)

	os.makedirs(data_dir, exist_ok=True)
	os.makedirs(key_dir, exist_ok=True)

	if not os.path.isfile(_SQLD_BIN) or not os.access(_SQLD_BIN, os.X_OK):
		print(f"ERROR: sqld binary not found or not executable at {_SQLD_BIN}.", file=sys.stderr)
		print("This command expects sqld to be baked into the container image at build time; "
			  "it does not download it. Use a ksrates image built from the current Dockerfile.", file=sys.stderr)
		sys.exit(1)

	private_key, public_key = _generate_jwt_keys(key_dir)
	token = _generate_jwt(private_key)

	hostname = socket.gethostname()
	with open(address_file, "w") as f:
		f.write(f"{hostname}:{port}\n")
		f.write(f"{token}\n")
	os.chmod(address_file, 0o600)

	print(f"Server node: {hostname}")
	print(f"Port: {port}")
	print(f"Data directory: {data_dir}")
	print(f"Address file: {address_file}")
	print("Starting sqld...")
	sys.stdout.flush()

	# Replaces this process with sqld (rather than spawning a subprocess), so the container's
	# main process is sqld itself - matching the sbatch script's "exec sqld" and keeping signal
	# handling/PID-1 semantics simple, whether launched by SLURM, a container runtime, or Nextflow.
	os.execv(_SQLD_BIN, [
		_SQLD_BIN,
		f"--http-listen-addr=0.0.0.0:{port}",
		f"--db-path={data_dir}",
		f"--auth-jwt-key-file={public_key}",
		f"--max-response-size={max_response_size}",
		f"--max-total-response-size={max_total_response_size}",
	])
