import base64
import os
import platform
import shutil
import socket
import ssl
import sys
import tarfile
import tempfile
import urllib.request

import certifi
from cryptography.hazmat.primitives import serialization
from cryptography.hazmat.primitives.asymmetric.ed25519 import Ed25519PrivateKey

# Prefer a sqld binary already baked into a container image (see Dockerfile); fall back to a
# user-level cache so a manually pip-installed ksrates (no container, likely no write access to
# /usr/local/bin) can still download and reuse its own copy across projects.
_SQLD_VERSION = "libsql-server-v0.24.32"
_SQLD_BAKED_IN_BIN = "/usr/local/bin/sqld"
_SQLD_ASSET_BY_ARCH = {
	"x86_64": "libsql-server-x86_64-unknown-linux-gnu.tar.xz",
	"aarch64": "libsql-server-aarch64-unknown-linux-gnu.tar.xz",
}


def _sqld_cache_dir():
	"""
	Resolve the user-level cache directory for a downloaded sqld binary ($XDG_CACHE_HOME, or
	~/.cache otherwise). os.makedirs(..., exist_ok=True) at the call site creates every missing
	directory in this path, including ~/.cache itself if it doesn't exist yet.

	:return: absolute path to <cache base>/ksrates/sqld_bin
	:raises RuntimeError: if neither $XDG_CACHE_HOME nor $HOME/a passwd entry is available to
	                       resolve a real home directory (expanduser("~") then returns "~"
	                       unexpanded, e.g. in some minimal/restricted container environments)
	"""
	base = os.environ.get("XDG_CACHE_HOME") or os.path.expanduser("~/.cache")
	if not os.path.isabs(base):
		raise RuntimeError(
			"Could not resolve a cache directory for sqld: neither $XDG_CACHE_HOME nor $HOME is "
			"set to a usable path. Set one of them, or place a sqld binary at "
			f"{_SQLD_BAKED_IN_BIN} directly."
		)
	return os.path.join(base, "ksrates", "sqld_bin")


def _download_sqld(target_path):
	"""
	Download and extract the sqld release binary to target_path. Used when sqld isn't already
	baked into a container image.

	:param target_path: final path to place the executable sqld binary at
	:raises RuntimeError: on unsupported architecture or any download/extraction failure
	"""
	arch = platform.machine()
	asset = _SQLD_ASSET_BY_ARCH.get(arch)
	if asset is None:
		raise RuntimeError(f"Unsupported architecture for sqld: {arch}")

	url = f"https://github.com/tursodatabase/libsql/releases/download/{_SQLD_VERSION}/{asset}"
	print(f"sqld binary not found, downloading {_SQLD_VERSION} for {arch} to {target_path}...")

	cache_dir = os.path.dirname(target_path)
	os.makedirs(cache_dir, exist_ok=True)
	try:
		with tempfile.TemporaryDirectory(dir=cache_dir) as tmp_dir:
			tarball_path = os.path.join(tmp_dir, "sqld.tar.xz")
			try:
				# Explicitly use certifi's CA bundle rather than relying on urllib's default
				# system trust store lookup, which fails with CERTIFICATE_VERIFY_FAILED on some
				# module-system/custom-built Python installs (seen on real HPC clusters) even
				# when the network connection itself is fine.
				ssl_context = ssl.create_default_context(cafile=certifi.where())
				with urllib.request.urlopen(url, context=ssl_context) as response, open(tarball_path, "wb") as out_file:
					shutil.copyfileobj(response, out_file)
			except Exception as e:
				raise RuntimeError(
					f"Failed to download {asset}: {e}. This node may have no internet access, or a "
					f"certificate/proxy issue; download it manually and place the 'sqld' binary at "
					f"{target_path}."
				) from e

			extract_dir = os.path.join(tmp_dir, "extracted")
			os.makedirs(extract_dir)
			with tarfile.open(tarball_path, "r:xz") as tar:
				tar.extractall(extract_dir)

			for root, _, files in os.walk(extract_dir):
				if "sqld" in files:
					extracted_bin = os.path.join(root, "sqld")
					break
			else:
				raise RuntimeError(f"No 'sqld' binary found inside downloaded {asset}")

			# Write to a temp path in the cache dir, then atomically move into place, so a
			# concurrent launch never sees a partially-written binary.
			tmp_target = target_path + ".tmp"
			with open(extracted_bin, "rb") as src, open(tmp_target, "wb") as dst:
				dst.write(src.read())
			os.chmod(tmp_target, 0o755)
			os.replace(tmp_target, target_path)
	except RuntimeError:
		raise
	except Exception as e:
		raise RuntimeError(f"Failed to install sqld: {e}") from e


def _resolve_sqld_bin():
	"""
	Find the sqld binary to run, downloading it to a user-level cache if it isn't already
	available (e.g. baked into a container image).

	:return: path to an executable sqld binary
	:raises RuntimeError: if sqld can't be found or installed
	"""
	if os.path.isfile(_SQLD_BAKED_IN_BIN) and os.access(_SQLD_BAKED_IN_BIN, os.X_OK):
		return _SQLD_BAKED_IN_BIN

	cached_bin = os.path.join(_sqld_cache_dir(), "sqld")
	if not (os.path.isfile(cached_bin) and os.access(cached_bin, os.X_OK)):
		_download_sqld(cached_bin)

	return cached_bin


def _b64url(data):
	"""base64url per the JWT spec (RFC 7515 Appendix C): standard base64 with '+/' -> '-_' and no padding."""
	return base64.urlsafe_b64encode(data).rstrip(b'=').decode('ascii')


def _generate_jwt_keys(key_dir):
	"""
	Load (or generate, if not already present) the Ed25519 signing keypair used to authenticate
	clients. The private key persists across server restarts (kept in key_dir, next to the data
	directory) so tokens already handed out in old address files keep working.

	Uses the `cryptography` package rather than shelling out to openssl: openssl's raw-Ed25519
	signing (pkeyutl -rawin) is only available from OpenSSL 3.0 onward and varies across
	container base images.

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

	try:
		sqld_bin = _resolve_sqld_bin()
	except RuntimeError as e:
		print(f"ERROR: {e}", file=sys.stderr)
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

	# Replaces this process with sqld (rather than spawning a subprocess) so sqld becomes the
	# container's main process, keeping signal handling/PID-1 semantics simple.
	os.execv(sqld_bin, [
		sqld_bin,
		f"--http-listen-addr=0.0.0.0:{port}",
		f"--db-path={data_dir}",
		f"--auth-jwt-key-file={public_key}",
		f"--max-response-size={max_response_size}",
		f"--max-total-response-size={max_total_response_size}",
	])
