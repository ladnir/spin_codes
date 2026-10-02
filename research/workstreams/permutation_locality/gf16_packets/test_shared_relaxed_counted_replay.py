from fractions import Fraction as Q
import hashlib
import json
from pathlib import Path
import unittest
from unittest.mock import patch

import shared_relaxed_strategy as s
import shared_relaxed_alternative_counts as counts


class CountedReplayTests(unittest.TestCase):
    def record(self):
        caps = [0, 2, 1]
        return dict(pruned=True, zero_bits=64, cost_tilt='1/4',
            mixture=[dict(mass='4', activity='1/2')],
            count_witnesses=[dict(path='fake-exact-dual.json', sha256=hashlib.sha256(b'dual').hexdigest())],
            count_refinement_iterations=0,
            cap_sha256=hashlib.sha256(json.dumps(list(map(str, caps))).encode()).hexdigest())

    def test_fresh_replay_receives_saved_rationals_not_saved_caps(self):
        record = self.record()
        result = (['fresh-components'], [0, 2, 1], [(Q(4), Q(1, 2))],
                  dict(count_witnesses=record['count_witnesses']))
        with patch.object(Path, 'read_bytes', return_value=b'dual'), \
                patch.object(counts, 'build_components', return_value=result) as fresh:
            self.assertEqual(s.checked_comparison(record), ['fresh-components'])
            fresh.assert_called_once_with([Path('fake-exact-dual.json')], iterations=0,
                step=8, zero_bits=64, cost_tilt=Q(1, 4), saved_mixture=record['mixture'])

    def test_changed_file_rejected_before_count_reconstruction(self):
        with patch.object(Path, 'read_bytes', return_value=b'changed'), \
                patch.object(counts, 'build_components') as fresh:
            with self.assertRaises(ValueError): s.checked_comparison(self.record())
            fresh.assert_not_called()

    def test_changed_cap_digest_or_replayed_source_rejected(self):
        record = self.record()
        for caps, sources in (([0, 3, 1], record['count_witnesses']),
                             ([0, 2, 1], [dict(path='fake-exact-dual.json', sha256='f'*64)])):
            with patch.object(Path, 'read_bytes', return_value=b'dual'), \
                    patch.object(counts, 'build_components', return_value=([], caps, [], dict(count_witnesses=sources))):
                with self.assertRaises(ValueError): s.checked_comparison(record)

    def test_malformed_sources_or_iteration_scope_rejected(self):
        for field, value in (('count_witnesses', None), ('count_witnesses', [{}]),
                ('count_refinement_iterations', True), ('count_refinement_iterations', -1),
                ('count_refinement_iterations', 31), ('pruned', False)):
            record = self.record(); record[field] = value
            with self.assertRaises(ValueError): s.checked_comparison(record)
        record = self.record(); record['count_witnesses'] = []; record['count_refinement_iterations'] = 1
        with self.assertRaises(ValueError): s.checked_comparison(record)
        record = self.record(); record['count_witnesses'] *= 2
        with patch.object(Path, 'read_bytes', return_value=b'dual'):
            with self.assertRaises(ValueError): s.checked_comparison(record)


if __name__ == '__main__': unittest.main()
