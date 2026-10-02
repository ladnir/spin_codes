from fractions import Fraction as Q
from itertools import product
from pathlib import Path
import subprocess
import sys
import unittest
from unittest.mock import patch

import numpy as np
from flint import arb, ctx
import count_caps_fast as fast


def mgfs(atom=False):
    return [dict(tilt=str(t), dual=['0', '0', '0'], **(
        dict(tilted_atom=True, tilted_variance_dual=['0', '0', '0']) if atom else {}))
        for t in (Q(-4), Q(-1), Q(0), Q(1), Q(4))]


def exact_atoms(probabilities):
    result = [Q(0)] * (len(probabilities)+1)
    for bits in product((0, 1), repeat=len(probabilities)):
        value = Q(1)
        for bit, p in zip(bits, probabilities):
            value *= p if bit else 1-p
        result[sum(bits)] += value
    return result


def endpoint(value):
    m, e = value.upper().man_exp()
    return Q(int(m)) * Q(2)**int(e)


class CountCapsFastTests(unittest.TestCase):
    def setUp(self):
        previous = ctx.prec
        self.addCleanup(setattr, ctx, 'prec', previous)
        ctx.prec = 256

    def arguments(self, atom=False):
        return ([Q(1, 4)], [1], (Q(1, 5), Q(3, 10)),
            (Q(1, 8), Q(1, 4)), 0, 8, arb(8), mgfs(atom))

    def test_cold_import_preserves_path_and_packed_dense_module(self):
        # A cold child is needed: reloading after discovery would reuse the
        # already-imported legacy dependencies and miss their path changes.
        script = """
from pathlib import Path
import sys
packed = Path(sys.argv[1]).resolve()
sys.path.insert(0, str(packed))
before = sys.path[:]
import count_caps_fast
assert sys.path == before, (before, sys.path)
import dense_cover
assert Path(dense_cover.__file__).resolve() == packed/'dense_cover.py'
assert dense_cover.SCHEMA == 'packed-canonical-gl32-dense-cover-1'
assert callable(dense_cover.fingerprint)
"""
        result = subprocess.run([sys.executable, '-B', '-c', script,
            str(Path(__file__).resolve().parent)], capture_output=True, text=True, timeout=30)
        self.assertEqual(result.returncode, 0, result.stdout+result.stderr)

    def assert_upper(self, bounds, probabilities):
        for j, value in enumerate(exact_atoms(probabilities)):
            self.assertGreaterEqual(endpoint(bounds[j]), value, j)

    def test_small_exact_enumerations_homogeneous_and_heterogeneous(self):
        for probabilities in ([Q(1, 4)]*6, [Q(1, 3), Q(1, 2), Q(2, 3)]*2,
                [Q(1, 8), Q(7, 8), Q(1, 4), Q(3, 4)],
                [Q(0), Q(1), Q(1, 4), Q(3, 4)]):
            groups = len(probabilities)
            mean = sum(probabilities)/groups
            variance = sum(p*(1-p) for p in probabilities)/groups
            args = (probabilities, [1]*groups, (mean, mean),
                (variance, variance), 0, groups, arb(2)**groups, mgfs(True))
            result = fast.count_mass_caps(*args)
            self.assert_upper(result, probabilities)

    def test_empty_family_preserves_exact_baseline(self):
        args = list(self.arguments())
        args[-1] = []
        with patch.object(fast, '_select') as select:
            actual = fast.count_mass_caps(*args)
        expected = fast.rc.count_mass_caps(*args)
        self.assertEqual(actual, expected)
        select.assert_not_called()

    def test_adversarial_valid_selector_can_only_weaken(self):
        args = self.arguments(True)
        for choices in ([0]*9, [4]*9, [4, 0, 3, 1, 2, 4, 0, 3, 1]):
            with patch.object(fast, '_select', return_value=choices):
                result = fast.count_mass_caps(*args)
            self.assert_upper(result, [Q(1, 4)]*8)
            self.assertTrue(all(value.is_finite() and value > 0 for value in result))

    def test_high_precision_original_candidates_are_below_selected_upper(self):
        args = self.arguments(True)
        chosen = fast.count_mass_caps(*args)
        ctx.prec = 768
        ordinary = fast.rc.count_mass_caps(*args)
        # Selected candidates form a subset of the independently checked
        # upper bounds. Higher precision makes endpoint-rounding effects
        # negligible; this fixture has strict rounding room even at ties.
        for actual, reference in zip(chosen, ordinary):
            self.assertGreaterEqual(endpoint(actual), endpoint(reference))
        self.assert_upper(chosen, [Q(1, 4)]*8)

    def test_tied_and_float_overflow_scores_remain_safe(self):
        args = list(self.arguments())
        args[-1] = [dict(tilt='0', dual=['0', '0', '0'])]*3
        result = fast.count_mass_caps(*args)
        self.assert_upper(result, [Q(1, 4)]*8)
        choices = fast._select([arb(10000).exp(), arb(1)], [Q(1000), Q(0)], 12)
        self.assertEqual(int(choices[0]), 1)
        self.assertEqual(int(choices[12]), 0)
        # All unusable float slopes fall back to a valid candidate index.
        huge = Q(10)**400
        self.assertEqual(list(fast._select([arb(1)], [huge], 4)), [0]*5)

    def test_extreme_output_mgf_no_float_moment_conversion(self):
        args = list(self.arguments())
        args[-1] = [dict(tilt=str(t), dual=['0', '0', '0']) for t in (-1000, 0, 1000)]
        result = fast.count_mass_caps(*args)
        self.assert_upper(result, [Q(1, 4)]*8)

    def test_invalid_selector_indices_fail_closed(self):
        for choices in ([0]*8, [0]*10, [-1]*9, [5]*9, [0.0]*9, [True]*9):
            with patch.object(fast, '_select', return_value=choices), self.assertRaises(ValueError):
                fast.count_mass_caps(*self.arguments())

    def test_every_candidate_checked_even_when_unselected(self):
        args = list(self.arguments())
        args[-1].append(dict(tilt='1', dual=['0', '0', '1']))
        with patch.object(fast, '_select', return_value=np.zeros(9, dtype=np.int64)) as select:
            with self.assertRaises(ValueError):
                fast.count_mass_caps(*args)
        select.assert_not_called()

    def test_malformed_domains_witnesses_and_caps_rejected(self):
        changes = [(0, []), (1, []), (2, (Q(0), Q(1, 2))),
            (3, (Q(0), Q(1, 3))), (4, 9), (5, 0), (5, True),
            (6, arb(0)), (6, arb(-1)), (6, arb('inf')), (6, arb('nan')),
            (7, [dict(tilt='0', dual=['0', '0'])]),
            (7, [dict(tilt='0', dual=['0', '0', '0'], tilted_atom=1)]),
            (7, [dict(tilt='0', dual=['0', '0', '0'], tilted_atom=True,
                      tilted_variance_dual=['0', '0', '-1'])])]
        for index, value in changes:
            args = list(self.arguments())
            args[index] = value
            with self.subTest(index=index, value=value), self.assertRaises((ValueError, ArithmeticError)):
                fast.count_mass_caps(*args)

    def test_invalid_candidate_arithmetic_rejected_even_if_not_selected(self):
        for value in (arb('inf'), arb('nan')):
            with patch.object(fast.rc, 'mgf_upper', return_value=value):
                with self.assertRaises((ValueError, ArithmeticError)):
                    fast.count_mass_caps(*self.arguments())
        with patch.object(fast.rc, 'tilted_atom_upper', return_value=arb(0)):
            with self.assertRaises(ArithmeticError):
                fast.count_mass_caps(*self.arguments(True))


if __name__ == '__main__':
    unittest.main()
