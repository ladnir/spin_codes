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
    @classmethod
    def setUpClass(cls) -> None:
        # Every test reads the same authenticated import; the corruption test
        # uses a separate in-memory database. Avoid rebuilding the large grid
        # for each independent read-only assertion.
        cls.temporary = tempfile.TemporaryDirectory()
        cls.path = pathlib.Path(cls.temporary.name) / "test.sqlite3"
        build_landscape_db.build(cls.path)

    def setUp(self) -> None:
        self.db = sqlite3.connect(f'{self.path.as_uri()}?mode=ro',uri=True)

    def tearDown(self) -> None:
        self.db.close()

    @classmethod
    def tearDownClass(cls) -> None:
        cls.temporary.cleanup()

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
        self.assertEqual(self.scalar("SELECT COUNT(*) FROM landscape WHERE comparison_eligible=1 AND occupation_min=1 AND source_path LIKE 'landscape_db/activation_pilot%'"), 1012)
        self.assertEqual(self.scalar("SELECT COUNT(*) FROM landscape WHERE comparison_eligible=1 AND dominant_witness_at_grid_edge!=0 AND source_path LIKE 'landscape_db/activation_pilot%'"), 0)
        self.assertEqual(self.scalar("SELECT COUNT(*) FROM landscape WHERE comparison_eligible=1 AND outer_model_kind='fixed_certified_constraints'"), 0)
        rows = self.db.execute("SELECT message_exponent,COUNT(*) FROM landscape WHERE comparison_eligible=1 AND occupation_min=1 AND source_path LIKE 'landscape_db/activation_pilot%' GROUP BY message_exponent").fetchall()
        self.assertEqual(rows, [(16, 242), (18, 242), (20, 242), (22, 143), (24, 143)])

    def test_historical_ids_stay_stable(self) -> None:
        self.assertEqual(self.scalar("SELECT COUNT(*) FROM results WHERE result_id='b08a619f03878eb7f7857ee6'"), 1)

    def test_random_conditional_event_is_explicit(self) -> None:
        row=dict(outer_model='random-simultaneous-spectrum-caps')
        with self.assertRaisesRegex(ValueError,'shared setup event'):
            build_landscape_db.outer_model_kind('random full-rank [32,16] reused',row)
        row.update(setup_event_id='test-event',setup_failure_bits=60)
        self.assertEqual(build_landscape_db.outer_model_kind('random full-rank [32,16] reused',row)[0],
                         'random_setup_spectrum_caps')
        self.assertEqual(self.scalar("SELECT COUNT(*) FROM landscape WHERE outer_model_kind='random_setup_spectrum_caps' AND (length(setup_event_id)!=64 OR setup_failure_bits!=60 OR setup_event_id IS NULL OR setup_failure_bits IS NULL)"),0)

    def test_small_state_maps_are_exact_prefixes(self) -> None:
        root=build_landscape_db.HERE
        for path in (root/'activation_pilot_small_state_v1/maps').glob('*.json'):
            payload=json.loads(path.read_text())
            parent=json.loads((root.parent.parent.parent/payload['parent_source']).read_text())
            self.assertEqual(payload['generator_words_hex'],parent['generator_words_hex'][:payload['state_bits']])

    def test_q2_has_realized_spectrum_and_pair_witness(self) -> None:
        rows=self.db.execute("SELECT outer_model_kind,dominant_weight,dominant_second_weight FROM landscape WHERE comparison_eligible=1 AND occupation_min=2 AND source_path LIKE 'landscape_db/activation_pilot%' ").fetchall()
        self.assertEqual(len(rows),36)
        self.assertTrue(all(kind=='fixed_exact_spectrum' and a>0 and b>0 for kind,a,b in rows))

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
