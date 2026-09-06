#!/usr/bin/env python3

from __future__ import annotations

import pathlib
import sqlite3
import sys
import tempfile
import unittest
import json

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

    def test_historical_certificate_preserved_but_excluded(self) -> None:
        row = self.db.execute(
            """
            SELECT message_bits,output_bits,occupation_min,occupation_max,margin_bits
            FROM results_under_review WHERE coverage_kind='full_distance'
            """
        ).fetchone()
        self.assertIsNotNone(row)
        self.assertEqual(row[:4], (65536, 131072, 1, 256))
        self.assertGreater(row[4], 42.57)
        self.assertEqual(self.scalar("SELECT COUNT(*) FROM certified_results"), 0)
        self.assertEqual(self.scalar("SELECT COUNT(*) FROM landscape WHERE transfer_review_status='activation_under_reaudit'"), 7)

    def test_sources_are_hashed(self) -> None:
        self.assertEqual(self.scalar("SELECT COUNT(*) FROM sources WHERE length(sha256)!=64"), 0)

    def test_balanced_preferred_pilot(self) -> None:
        self.assertEqual(self.scalar("SELECT COUNT(*) FROM landscape WHERE comparison_eligible=1"), 715)
        self.assertEqual(self.scalar("SELECT COUNT(*) FROM landscape WHERE comparison_eligible=1 AND dominant_witness_at_grid_edge!=0"), 0)
        self.assertEqual(self.scalar("SELECT COUNT(*) FROM landscape WHERE comparison_eligible=1 AND outer_model_kind='fixed_certified_constraints'"), 0)
        rows = self.db.execute("SELECT message_exponent,COUNT(*) FROM landscape WHERE comparison_eligible=1 GROUP BY message_exponent").fetchall()
        self.assertEqual(rows, [(16, 143), (18, 143), (20, 143), (22, 143), (24, 143)])

    def test_historical_ids_stay_stable(self) -> None:
        self.assertEqual(self.scalar("SELECT COUNT(*) FROM results WHERE result_id='b08a619f03878eb7f7857ee6'"), 1)

    def test_changed_screen_csv_is_rejected(self) -> None:
        catalog = json.loads(build_landscape_db.CATALOG.read_text())
        spec = dict(catalog['activation_sources'][0])
        manifest = json.loads((build_landscape_db.HERE / spec['manifest']).read_text())
        manifest['csv_sha256'] = '0' * 64
        path = pathlib.Path(self.temporary.name) / 'bad_manifest.json'
        path.write_text(json.dumps(manifest))
        spec['manifest'] = str(path)
        fresh = sqlite3.connect(':memory:')
        try:
            fresh.executescript(build_landscape_db.SCHEMA.read_text())
            with self.assertRaisesRegex(ValueError, 'changed screen CSV'):
                build_landscape_db.import_csv(fresh, spec)
        finally:
            fresh.close()


if __name__ == "__main__":
    unittest.main()
