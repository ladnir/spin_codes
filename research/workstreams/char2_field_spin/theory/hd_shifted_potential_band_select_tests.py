from fractions import Fraction as F
import json
from pathlib import Path
import unittest

from hd_shifted_potential_band_select import select_interval
from hd_shifted_swap_screen import CONSTRUCTION


class SelectionTests(unittest.TestCase):
    def fixture(self):
        return dict(construction=CONSTRUCTION, directed=False, every_weight_evaluated=True,
            requested_weight_range=[1, 4], max_weights=4,
            results=[dict(active_groups=r, log2_upper_screen=-100.) for r in range(1, 5)],
            fixed_witness=dict(source_exact_witness=dict(z='1/2')))

    def test_selection_is_complete_and_preserves_source(self):
        source = self.fixture()
        result = select_interval(source, 2, 3)
        self.assertEqual(result['requested_weight_range'], [2, 3])
        self.assertEqual(result['weights_evaluated'], 2)
        self.assertEqual(result['union_log2_upper_screen'], -99.)
        self.assertEqual(result['new_fft_tables_built'], 0)
        self.assertFalse(result['directed'])
        self.assertEqual(result['fixed_witness'], source['fixed_witness'])
        self.assertEqual(len(source['results']), 4)

    def test_missing_duplicate_or_outside_weight_rejected(self):
        for rows in ([1, 3, 4], [1, 2, 2, 3, 4]):
            source = self.fixture()
            source['results'] = [dict(active_groups=r, log2_upper_screen=-100.) for r in rows]
            with self.assertRaises(ValueError):
                select_interval(source, 1, 4)
        with self.assertRaises(ValueError):
            select_interval(self.fixture(), 1, 5)

    def test_optional_broad_screen_has_one_maximal_passing_interval(self):
        path = Path(__file__).with_name('hd_shifted_potential_band_512_400_600.json')
        if not path.exists():
            self.skipTest('Optional broad FLOAT screen absent; tests never rebuild an atlas')
        source = json.loads(path.read_text(encoding='utf-8'))
        self.assertEqual([row['active_groups'] for row in source['results']], list(range(400, 601)))
        # A failing individual weight cannot belong to a passing union.
        self.assertEqual([row['active_groups'] for row in source['results']
                          if row['log2_upper_screen'] <= -80.], list(range(458, 567)))
        selected = select_interval(source, 458, 566)
        self.assertLess(selected['union_log2_upper_screen'], -80.)
        self.assertEqual(source['fft_tables_built'], 19)
        self.assertEqual(source['witness_searches'], 0)

    def test_optional_wide_receipts_preserve_precisions_and_narrow_replay(self):
        paths = [Path(__file__).with_name('hd_shifted_potential_band_512_458_566_directed%d.json' % p)
                 for p in (256, 384)]
        narrow = Path(__file__).with_name('hd_shifted_potential_band_512_496_528_directed256.json')
        if not all(path.exists() for path in (*paths, narrow)):
            self.skipTest('Optional directed receipts absent; tests never run a proof replay')
        receipts = [json.loads(path.read_text(encoding='utf-8')) for path in paths]
        for receipt in receipts:
            self.assertEqual(receipt['max_weights'], 109)
            self.assertEqual(receipt['weights_replayed'], 109)
            self.assertEqual(receipt['atlas_builds'], 19)
            self.assertEqual(receipt['integer_profile_triples'], 7212399)
            self.assertEqual([row['active_groups'] for row in receipt['results']], list(range(458, 567)))
            self.assertEqual(F(receipt['log2_upper_rational']), F(-49106949211489, 500000000000))
            self.assertTrue(receipt['directed'])
            self.assertFalse(receipt['whole_code_certificate'])
        maps = [{row['active_groups']: row['log2_upper_rational'] for row in r['results']} for r in receipts]
        self.assertEqual(*maps)
        prior = json.loads(narrow.read_text(encoding='utf-8'))
        for row in prior['results']:
            self.assertEqual(maps[0][row['active_groups']], row['log2_upper_rational'])

    def test_optional_two_central_banks_have_complete_maximal_intervals(self):
        cases = ((320, 257, 457, 280, 349), (768, 567, 974, 700, 809))
        paths = [Path(__file__).with_name('hd_shifted_potential_band_%d_%d_%d.json' % case[:3])
                 for case in cases]
        if not all(path.exists() for path in paths):
            self.skipTest('Optional bounded two-bank screens absent; no fitting occurs in tests')
        for (_, first, last, low, high), path in zip(cases, paths):
            source = json.loads(path.read_text(encoding='utf-8'))
            self.assertEqual([row['active_groups'] for row in source['results']], list(range(first, last+1)))
            self.assertEqual(source['weights_evaluated'], last-first+1)
            self.assertEqual(source['neighbor_witness_searches'], 0)
            self.assertEqual(source['fft_tables_built'], 22)
            self.assertEqual(source['central_marker_fits'], 22)
            self.assertEqual([row['active_groups'] for row in source['results']
                              if row['log2_upper_screen'] <= -80.], list(range(low, high+1)))
            self.assertLess(select_interval(source, low, high)['union_log2_upper_screen'], -80.)
            self.assertFalse(source['directed'])


if __name__ == '__main__':
    unittest.main()
