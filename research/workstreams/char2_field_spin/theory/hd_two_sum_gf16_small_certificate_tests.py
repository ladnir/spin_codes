"""Exact/light all-source adapter tests; no production regional bank builds."""
from collections import Counter
from copy import deepcopy
from fractions import Fraction as F
from hashlib import sha256
from math import comb
from pathlib import Path
import unittest
from unittest.mock import patch

from flint import arb, ctx

import hd_two_sum_gf16_small_certificate as target
from hd_two_sum_affine_sparse_certificate_tests import encloses_fraction
from hd_xor_dense_holder_certificate import dyadic_power_dominates
from hd_xor_joint_columns import categories


class Gf16SmallCertificateTests(unittest.TestCase):
    def test_frozen322_partition_excludes_only_r1_dense(self):
        selected = target.checked_selection()
        self.assertEqual(len(target.ALL_KEYS), 323)
        self.assertEqual(len(selected), 322)
        self.assertNotIn((1, 0, 0, 1), selected)
        self.assertEqual({key for key in selected if key[0] == 1}, {(1, 1, 0, 0), (1, 0, 1, 0)})
        self.assertEqual(sum(key[0] >= 2 for key in selected), 320)
        self.assertEqual(Counter(selected.values()), {'dense-inherited': 5, 'h2-inherited': 317})
        for r, n1, n2, dense in selected:
            fields = target.cell_scope((r, n1, n2, dense))
            self.assertEqual(fields['source_supports'], comb(8192, r))
            self.assertEqual(fields['source_colorings'], comb(r, dense)*comb(r-dense, n1))
        with patch.object(target, 'PROPOSAL_SHA256', '0'*64), self.assertRaises(ValueError):
            target.checked_selection()

    def test_combined_exact_moments_and_root_directions(self):
        nums, den, _ = target.combined_moments()
        self.assertEqual(den, 255**128)
        with ctx.workprec(256):
            weights, _ = target.combined_weights_arb()
            for w in range(17):
                exact = sum((F(target.arity_holder_moments(a, 16, 16, 8)[0][w], 255**(8*a))
                             for a in range(3, 17)), F(0))
                self.assertEqual(F(nums[w], den), exact)
                if nums[w]:
                    self.assertTrue(dyadic_power_dominates(weights[w], nums[w], den, 16))

    def test_maximum_uses_upper_endpoints_for_every_feasible_f(self):
        raw = {(m, f): tuple(arb(f'{1+m+f+i} +/- 0.3') for i in range(4))
               for m in range(5) for f in range(m//2+1)}
        envelope = target.all_source_envelope(raw, 4)
        self.assertEqual(len(envelope), 16)
        for (m, f), matrix in raw.items():
            for i, value in enumerate(matrix):
                self.assertTrue(envelope[0][m][i] >= value.upper())
        self.assertTrue(all(row == envelope[0] for row in envelope))

    def test_dense_and_sparse_floor_normalizations_once(self):
        caps = {kind: arb(1) for kind in categories(17)}
        z, potential = F(1, 2), arb(1)/4
        with ctx.workprec(256):
            with patch.object(target, 'bivariate_tail_arb', return_value=(arb(7), [])) as tail:
                actual = target.cell_value((2, 0, 1, 1), z, potential, caps, (caps, caps))
            self.assertEqual(tail.call_args.kwargs, dict(columns=16, x_floor=18, y_floor=9))
            expected = F(comb(8192, 2)*2*120*2**734003*4*7, 255**17)
            self.assertTrue(encloses_fraction(actual, expected))
            with patch.object(target, 'sparse_tail_arb', return_value=arb(11)) as sparse_tail:
                actual = target.cell_value((2, 1, 1, 0), z, potential, caps, None)
            self.assertEqual(sparse_tail.call_args.kwargs, dict(columns=16, dimension=8))
            expected = F(comb(8192, 2)*2*16*120*2**734003*4*11, 255)
            self.assertTrue(encloses_fraction(actual, expected))
        with self.assertRaises(ValueError):
            target.cell_value((1, 0, 0, 1), z, potential, caps, (caps, caps))

    def test_only_selected_banks_built_and_no_float_bound_consumed(self):
        with patch.object(target, 'combined_weights_arb', return_value=({0: arb(1)}, 'toy')), \
                patch.object(target, 'arb_pair_regions', return_value=None) as raw, \
                patch.object(target, 'all_source_envelope', return_value=None), \
                patch.object(target, 'category_caps_arb', return_value=None) as sparse, \
                patch.object(target, 'weight_vectors_arb', return_value=None) as vectors, \
                patch.object(target, 'holder_branch_caps_arb', return_value=None), \
                patch.object(target, 'cell_value', return_value=arb(2)**-100) as cell:
            result = target.receipt()
        self.assertEqual(raw.call_count, 2)
        self.assertEqual(sparse.call_count, 2)
        self.assertEqual(vectors.call_count, 2)
        self.assertEqual(cell.call_count, 322)
        self.assertEqual({call.args[0] for call in cell.call_args_list}, set(target.KEYS))
        self.assertEqual([row['cells'] for row in result['bank_costs']], [5, 317])
        self.assertFalse(result['saved_float_estimates_used'])
        with ctx.workprec(512):
            expected = target.ceiling((arb(322)*arb(2)**-100).log()/arb(2).log())
        self.assertEqual(F(result['failure_log2_upper']), expected)

    def fixture(self):
        selection = target.checked_selection()
        entries = []
        for key, name in selection.items():
            z, p = target.LITERALS[name]
            entries.append(dict(target.cell_scope(key), witness=name, mark=str(z), log2_potential=str(p),
                                failure_log2_upper='-100'))
        return dict(target.EXPECTED_SCOPE, precision_bits=256, frozen_proposal_sha256=target.PROPOSAL_SHA256,
            reviewed_helper_sha256=target.HELPER_PINS, combined_moments_sha256=target.combined_moments()[2],
            root_upper_checks_passed=True, source_sha256=sha256(Path(target.__file__).read_bytes()).hexdigest(),
            entries=entries, by_active_groups=[dict(active_groups=r, failure_log2_upper='-95') for r in range(1, 18)])

    def test_scope_and_cell_assignments_fail_closed(self):
        data = self.fixture()
        self.assertEqual(len(target.checked_rows(data, 256)[0]), 322)
        for key, value in (('whole_code_certificate', True), ('one_dense_group_range', [1, 17]),
                           ('merged_arity_moments', False), ('selected_cells', 323)):
            with self.assertRaises(ValueError):
                target.checked_rows(dict(data, **{key: value}), 256)
        for field, value in (('cohort_counts', [0, 0, 1]), ('witness', 'sparse-inherited'),
                             ('source_supports', 1), ('sparse_normalizer_exponent', 0)):
            changed = deepcopy(data)
            changed['entries'][0][field] = value
            with self.assertRaises(ValueError):
                target.checked_rows(changed, 256)
        changed = deepcopy(data)
        changed['entries'][-1] = changed['entries'][0]
        with self.assertRaises(ValueError):
            target.checked_rows(changed, 256)

    def test_budget_and_precision_fail_before_operator(self):
        for p in (True, 128, 512):
            with self.assertRaises(ValueError):
                target.receipt(p)
        for seconds in (True, 0, 601, float('nan'), float('inf')):
            with self.assertRaises(ValueError):
                target.receipt(max_seconds=seconds)
        with patch.object(target, 'arb_pair_regions') as bank, self.assertRaises(TimeoutError):
            target.receipt(max_seconds=1e-12)
        bank.assert_not_called()


if __name__ == '__main__':
    unittest.main()
