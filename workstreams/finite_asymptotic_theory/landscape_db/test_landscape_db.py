#!/usr/bin/env python3

from __future__ import annotations

import pathlib
import sqlite3
import sys
import tempfile
import unittest

sys.path.insert(0, str(pathlib.Path(__file__).resolve().parent))
import build_landscape_db


class LandscapeDatabaseTest(unittest.TestCase):
    def setUp(self) -> None:
        self.temporary = tempfile.TemporaryDirectory()
        self.path = pathlib.Path(self.temporary.name) / "test.sqlite3"
        build_landscape_db.build(self.path)
        self.db = sqlite3.connect(self.path)

    def tearDown(self) -> None:
        self.db.close()
        self.temporary.cleanup()

    def scalar(self, query: str) -> object:
        return self.db.execute(query).fetchone()[0]

    def test_rm2sub_only(self) -> None:
        self.assertEqual(self.scalar("SELECT COUNT(*) FROM inner_configs WHERE family!='rm2sub'"), 0)

    def test_random_means_outer_reference(self) -> None:
        self.assertGreater(self.scalar("SELECT COUNT(*) FROM outer_models WHERE family='random'"), 0)
        self.assertEqual(
            self.scalar("SELECT COUNT(*) FROM results WHERE result_class='reference'"),
            self.scalar("SELECT COUNT(*) FROM landscape WHERE outer_family='random'"),
        )

    def test_full_certificate_present(self) -> None:
        row = self.db.execute(
            """
            SELECT message_bits,output_bits,occupation_min,occupation_max,margin_bits
            FROM certified_results WHERE coverage_kind='full_distance'
            """
        ).fetchone()
        self.assertIsNotNone(row)
        self.assertEqual(row[:4], (65536, 131072, 1, 256))
        self.assertGreater(row[4], 42.57)

    def test_sources_are_hashed(self) -> None:
        self.assertEqual(self.scalar("SELECT COUNT(*) FROM sources WHERE length(sha256)!=64"), 0)


if __name__ == "__main__":
    unittest.main()
