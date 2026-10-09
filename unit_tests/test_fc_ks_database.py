"""
Unit tests for the Ks database modules fc_ks_database_utilities, _paralog and _ortholog.
The sqld server is replaced by a local SQLite file (see fake_sqld.py).

Run with: python -m pytest unit_tests/test_fc_ks_database.py
"""
import os
import sys
import tempfile
import threading
import unittest

sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))
sys.path.insert(0, os.path.dirname(__file__))

import ksrates.fc_ks_database_utilities as fcKsDbUtil
import ksrates.fc_ks_database_paralog as fcKsDbPara
import ksrates.fc_ks_database_ortholog as fcKsDbOrtho
from fake_sqld import install_fake_sqld

SP_A, SP_B = "Zea mays", "Oryza sativa"
KEY = "Oryza sativa_Zea mays"  # alphabetically sorted


class KsDatabaseTest(unittest.TestCase):

    def setUp(self):
        self.tmpdir = tempfile.TemporaryDirectory()
        self.addCleanup(self.tmpdir.cleanup)
        self.db, self.unreachable = install_fake_sqld(self, self.tmpdir.name)
        fcKsDbUtil.initialize_database(self.db)

    # ---- core / ortholog table
    def test_pair_key_is_sorted_case_insensitive_and_symmetric(self):
        self.assertEqual(fcKsDbOrtho.pair_key(SP_A, SP_B), KEY)
        self.assertEqual(fcKsDbOrtho.pair_key(SP_B, SP_A), KEY)
        self.assertEqual(fcKsDbOrtho.pair_key("b sp", "A sp"), "A sp_b sp")

    def test_initialize_is_idempotent_and_adds_missing_columns_to_an_older_table(self):
        fcKsDbUtil.initialize_database(self.db)
        fcKsDbUtil.initialize_database(self.db)
        # Simulate a table created by an older schema (only the key column)
        with fcKsDbUtil.connection(self.db) as client:
            client.execute("DROP TABLE ortholog_ks")
            client.execute("CREATE TABLE ortholog_ks (pair_key TEXT PRIMARY KEY)")
        fcKsDbUtil.initialize_database(self.db)
        self.assertTrue(fcKsDbOrtho.write_to_ortholog_db(self.db, SP_A, SP_B, peak=(1.0, 0.1), ks_list=[0.5]))

    def test_partial_writes_never_wipe_other_columns(self):
        self.assertEqual(fcKsDbOrtho.existing_ortholog_data(self.db, SP_A, SP_B), {"peak": False, "ks_list": False})
        fcKsDbOrtho.write_to_ortholog_db(self.db, SP_A, SP_B, peak=(1.5, 0.1))
        self.assertEqual(fcKsDbOrtho.existing_ortholog_data(self.db, SP_B, SP_A), {"peak": True, "ks_list": False})
        fcKsDbOrtho.write_to_ortholog_db(self.db, SP_B, SP_A, ks_list=[0.1, 0.2, 0.3])
        self.assertEqual(fcKsDbOrtho.existing_ortholog_data(self.db, SP_A, SP_B), {"peak": True, "ks_list": True})
        self.assertEqual(fcKsDbOrtho.read_peak(self.db, SP_A, SP_B), (1.5, 0.1))
        self.assertEqual(fcKsDbOrtho.read_ks_list(self.db, SP_A, SP_B), [0.1, 0.2, 0.3])

    def test_nothing_to_write_returns_false(self):
        self.assertFalse(fcKsDbOrtho.write_to_ortholog_db(self.db, SP_A, SP_B))
        self.assertFalse(fcKsDbUtil.row_exists(self.db, fcKsDbOrtho.ORTHOLOG_TABLE, KEY))

    def test_reads_of_unknown_pair(self):
        self.assertIsNone(fcKsDbOrtho.read_peak(self.db, SP_A, SP_B))
        self.assertIsNone(fcKsDbOrtho.read_ks_list(self.db, SP_A, SP_B))

    def test_bulk_reads_have_legacy_shape(self):
        fcKsDbOrtho.write_to_ortholog_db(self.db, SP_A, SP_B, peak=(1.5, 0.1), ks_list=[0.1, 0.2])
        fcKsDbOrtho.write_to_ortholog_db(self.db, SP_A, "Aaa bbb", ks_list=[0.3])  # Ks list but no peak
        peaks = fcKsDbOrtho.read_all_peaks(self.db)
        self.assertEqual(list(peaks.index), [KEY])
        self.assertEqual(list(peaks.columns), ["Species1", "Species2", "Mode", "Mode_SD"])
        self.assertEqual(peaks.at[KEY, "Species1"], SP_B)  # alphabetically first
        ks_lists = fcKsDbOrtho.read_ks_lists(self.db, [KEY, "Aaa bbb_Zea mays", "unknown_pair"])
        self.assertEqual(set(ks_lists.index), {KEY, "Aaa bbb_Zea mays"})
        self.assertEqual(ks_lists.at[KEY, "Ks_Values"], [0.1, 0.2])  # real list, not a string
        self.assertEqual(fcKsDbOrtho.pair_keys_in_database(self.db, "peak"), {KEY})
        self.assertEqual(fcKsDbOrtho.pair_keys_in_database(self.db, "ks_list"), {KEY, "Aaa bbb_Zea mays"})

    def test_empty_database_gives_empty_frames_with_legacy_columns(self):
        peaks = fcKsDbOrtho.read_all_peaks(self.db)
        self.assertEqual(peaks.shape[0], 0)
        self.assertEqual(list(peaks.columns), ["Species1", "Species2", "Mode", "Mode_SD"])
        self.assertEqual(fcKsDbOrtho.read_ks_lists(self.db, [KEY]).shape[0], 0)

    def test_concurrent_writers_of_the_same_new_pair_leave_one_consistent_row(self):
        # The bug fixed by the database: concurrent jobs appending to TSV files.
        errors = []

        def write(**kwargs):
            try:
                fcKsDbOrtho.write_to_ortholog_db(self.db, SP_A, SP_B, **kwargs)
            except Exception as e:  # pragma: no cover
                errors.append(e)

        threads = [threading.Thread(target=write, kwargs={"peak": (1.5, 0.1)}) for _ in range(4)] + \
                  [threading.Thread(target=write, kwargs={"ks_list": [0.1, 0.2]}) for _ in range(4)] + \
                  [threading.Thread(target=write, kwargs={"peak": (1.5, 0.1), "ks_list": [0.1, 0.2]}) for _ in range(4)]
        for thread in threads:
            thread.start()
        for thread in threads:
            thread.join()
        self.assertEqual(errors, [])
        self.assertEqual(fcKsDbUtil.list_keys(self.db, fcKsDbOrtho.ORTHOLOG_TABLE), [KEY])
        self.assertEqual(fcKsDbOrtho.read_peak(self.db, SP_A, SP_B), (1.5, 0.1))
        self.assertEqual(fcKsDbOrtho.read_ks_list(self.db, SP_A, SP_B), [0.1, 0.2])

    def test_delete_keys(self):
        fcKsDbOrtho.write_to_ortholog_db(self.db, SP_A, SP_B, peak=(1.5, 0.1))
        fcKsDbUtil.delete_keys(self.db, fcKsDbOrtho.ORTHOLOG_TABLE, [KEY])
        self.assertEqual(fcKsDbUtil.list_keys(self.db, fcKsDbOrtho.ORTHOLOG_TABLE), [])

    def test_unknown_column_is_rejected(self):
        with self.assertRaises(ValueError):
            fcKsDbUtil.upsert(self.db, fcKsDbOrtho.ORTHOLOG_TABLE, KEY, {"nope; DROP TABLE x": 1})

    # ---- paralog table (ported code: same behavior as before)
    def test_paralog_roundtrip_and_partial_writes(self):
        self.assertFalse(fcKsDbPara.species_exists(self.db, SP_A))
        self.assertEqual(fcKsDbPara.existing_data_types(self.db, SP_A),
                         {"paranome": False, "anchors": False, "recret": False, "iadhore": False})
        ks_data = {"Ks": [1.0, 2.0], "Family": ["a", "b"]}
        self.assertTrue(fcKsDbPara.write_to_paralog_db(SP_A, {"paranome": ks_data, "anchors": None, "recret": None}, self.db))
        self.assertTrue(fcKsDbPara.species_exists(self.db, SP_A))
        self.assertEqual(fcKsDbPara.read_analysis_data(self.db, SP_A, "paranome"), ks_data)
        self.assertIsNone(fcKsDbPara.read_analysis_data(self.db, SP_A, "anchors"))
        # a later partial write must not wipe the paranome
        fcKsDbPara.write_to_paralog_db(SP_A, {"paranome": None, "anchors": {"Ks": [3.0]}, "recret": None}, self.db)
        self.assertEqual(fcKsDbPara.read_analysis_data(self.db, SP_A, "paranome"), ks_data)
        self.assertEqual(fcKsDbPara.read_analysis_data(self.db, SP_A, "anchors"), {"Ks": [3.0]})
        self.assertFalse(fcKsDbPara.write_to_paralog_db(SP_A, {"paranome": None, "anchors": None, "recret": None}, self.db))

    def test_paralog_iadhore_files(self):
        self.assertTrue(fcKsDbPara.write_anchor_iadhore_files(SP_A, {"anchorpoints": "text a", "segments": "text s"}, self.db))
        read = fcKsDbPara.read_anchor_iadhore_files(self.db, SP_A)
        self.assertEqual(read["anchorpoints"], "text a")
        self.assertIsNone(read["multiplicons"])
        self.assertFalse(fcKsDbPara.existing_data_types(self.db, SP_A)["iadhore"])  # needs all 5 files
        full = {key: key for key in fcKsDbPara.IADHORE_FILES}
        fcKsDbPara.write_anchor_iadhore_files(SP_A, full, self.db)
        self.assertTrue(fcKsDbPara.existing_data_types(self.db, SP_A)["iadhore"])

    # ---- fail-fast
    def test_unreachable_server_raises_everywhere(self):
        for call in (lambda: fcKsDbUtil.initialize_database(self.unreachable),
                     lambda: fcKsDbOrtho.read_peak(self.unreachable, SP_A, SP_B),
                     lambda: fcKsDbOrtho.read_all_peaks(self.unreachable),
                     lambda: fcKsDbOrtho.pair_keys_in_database(self.unreachable, "peak"),
                     lambda: fcKsDbPara.existing_data_types(self.unreachable, SP_A),
                     lambda: fcKsDbPara.read_analysis_data(self.unreachable, SP_A, "paranome"),
                     lambda: fcKsDbOrtho.write_to_ortholog_db(self.unreachable, SP_A, SP_B, peak=(1.0, 0.1)),
                     lambda: fcKsDbPara.write_to_paralog_db(SP_A, {"paranome": {"Ks": [1]}, "anchors": None, "recret": None}, self.unreachable)):
            with self.assertRaises(fcKsDbUtil.KsDatabaseError):
                call()

    def test_missing_or_malformed_address_file_raises(self):
        with self.assertRaises(fcKsDbUtil.KsDatabaseError):
            fcKsDbUtil.initialize_database(os.path.join(self.tmpdir.name, "does_not_exist.txt"))
        malformed = os.path.join(self.tmpdir.name, "malformed.txt")
        with open(malformed, "w") as handle:
            handle.write("no-port-here\n")
        with self.assertRaises(fcKsDbUtil.KsDatabaseError):
            fcKsDbUtil.initialize_database(malformed)

    def test_initialize_or_exit_exits_with_error_code_1(self):
        with self.assertRaises(SystemExit) as context:
            fcKsDbUtil.initialize_or_exit(self.unreachable)
        self.assertEqual(context.exception.code, 1)
        fcKsDbUtil.initialize_or_exit(self.db)  # a working one doesn't exit

    # ---- legacy TSV import
    def test_extract_from_legacy_tsvs_merges_and_never_recomputes(self):
        peak_tsv = os.path.join(self.tmpdir.name, "peak.tsv")
        ks_tsv = os.path.join(self.tmpdir.name, "ks.tsv")
        with open(peak_tsv, "w") as handle:
            handle.write("\tSpecies1\tSpecies2\tMode\tMode_SD\n")
            handle.write(f"{KEY}\t{SP_B}\t{SP_A}\t1.23\t0.1\n")
            handle.write("A_B\tA\tB\t2.0\t0.2\n")  # peak only
        with open(ks_tsv, "w") as handle:
            handle.write("\tSpecies1\tSpecies2\tKs_Values\n")
            handle.write(f"{KEY}\t{SP_B}\t{SP_A}\t[1.1, 1.2]\n")
            handle.write("C_D\tC\tD\t[1.1, 1.2\n")  # truncated (the old concurrent-append bug): skipped, not fatal
        pairs = fcKsDbOrtho.extract_from_legacy_tsvs(peak_tsv, ks_tsv)
        self.assertEqual(pairs[KEY]["peak"], (1.23, 0.1))
        self.assertEqual(pairs[KEY]["ks_list"], [1.1, 1.2])
        self.assertIsNone(pairs["A_B"]["ks_list"])
        self.assertNotIn("C_D", pairs)


if __name__ == "__main__":
    unittest.main()
