"""
Test helper: stands in for the sqld server in the unit tests, so that they don't need a running server.

install_fake_sqld() replaces libsql_client.create_client_sync() (as used by ksrates.fc_ks_database_utilities) with a
client backed by a local SQLite file, which understands the same SQL we send (CREATE TABLE IF NOT EXISTS,
pragma_table_info(), INSERT ... ON CONFLICT DO UPDATE). Every "connection" opens the same file, like the
connections of independent jobs to the same server. Address files whose host is "unreachable" make every
query fail, to simulate a server that is down.
"""
import os
import sqlite3
from types import SimpleNamespace
from unittest import mock

import ksrates.fc_ks_database_utilities as fcKsDbUtil




class _FakeClient:
    def __init__(self, address, sqlite_path):
        self.address, self.sqlite_path = address, sqlite_path

    def __enter__(self):
        return self

    def __exit__(self, *args):
        return False

    def execute(self, sql, params=()):
        if self.address.startswith("unreachable"):
            raise ConnectionError(f"Cannot connect to host {self.address}")
        connection = sqlite3.connect(self.sqlite_path, timeout=30)
        try:
            cursor = connection.execute(sql, params)
            rows = cursor.fetchall()
            connection.commit()
            return SimpleNamespace(rows=rows)
        finally:
            connection.close()


def install_fake_sqld(testcase, tmpdir):
    """
    Patches the database client for the duration of a unittest.TestCase test and writes two address files in
    tmpdir. Returns (working_address_file, unreachable_address_file).
    """
    sqlite_path = os.path.join(tmpdir, "fake_sqld.sqlite")

    def create_client_sync(url, auth_token=None):
        return _FakeClient(url.replace("ws://", ""), sqlite_path)

    patcher = mock.patch.object(fcKsDbUtil.libsql_client, "create_client_sync", create_client_sync)
    patcher.start()
    testcase.addCleanup(patcher.stop)

    working = os.path.join(tmpdir, "ks_database_server_address.txt")
    unreachable = os.path.join(tmpdir, "unreachable_address.txt")
    with open(working, "w") as handle:
        handle.write("localhost:8080\n")
    with open(unreachable, "w") as handle:
        handle.write("unreachable:8080\n")
    return working, unreachable
