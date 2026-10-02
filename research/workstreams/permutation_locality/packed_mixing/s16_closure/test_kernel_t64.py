"""Exact small-state checks of the explicit two-physical-step macro kernel."""
from copy import deepcopy
from fractions import Fraction as Q
from itertools import combinations
import json
from math import comb
from pathlib import Path
from tempfile import TemporaryDirectory
import unittest

import numpy as np
from flint import arb, arb_mat, ctx
import kernel_t64 as kernel
import kernel_maps
import scalar_cover as sc
import sparse_kernel


def endpoint(value):
    mantissa, exponent = value.upper().man_exp()
    return Q(int(mantissa))*Q(2)**int(exponent)


def feedback(x, columns):
    value = 0
    for bit, column in enumerate(columns):
        if x >> bit & 1:
            value ^= column
    return value


def physical_exact(images, columns, occupancy, z):
    """Four actual states, one GF16 packet, exact rational arithmetic."""
    inputs = [0] if occupancy == 0 else range(1, 16)
    S = len(images)
    result = [[Q(0)]*S for _ in range(S)]
    for q in range(S):
        refreshed = range(1, S) if q else [0]
        for x in inputs:
            mass = z**(x ^ images[q]).bit_count()/len(inputs)/len(refreshed)
            for fresh in refreshed:
                result[q][fresh ^ feedback(x, columns)] += mass
    return result


def macro_direct(images, columns, occupancy, z):
    """Enumerate both input words and both refreshes without matrix composition."""
    S = len(images)
    words = [(x, y) for x in range(16) for y in range(16)
             if bool(x)+bool(y) == occupancy]
    result = [[Q(0)]*S for _ in range(S)]
    for q in range(S):
        first = range(1, S) if q else [0]
        for x, y in words:
            for fresh1 in first:
                middle = fresh1 ^ feedback(x, columns)
                second = range(1, S) if middle else [0]
                mass = z**((x ^ images[q]).bit_count()+(y ^ images[middle]).bit_count())
                mass /= len(words)*len(first)*len(second)
                for fresh2 in second:
                    result[q][fresh2 ^ feedback(y, columns)] += mass
    return result


def arb_matrix(rows):
    return arb_mat([[kernel.aq(value) for value in row] for row in rows])


class MacroKernelTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        ctx.prec = 192
        cls.images, cls.columns = (0, 1, 6, 7), [1, 2, 3, 1]
        cls.physical = {
            density: kernel_maps.prepare_maps(cls.images, cls.columns, bits=2,
                distribution='uniform_gl', birth_density=density)
            for density in ('classes', 'capped')}

    def setUp(self):
        ctx.prec = 192

    def assert_near(self, left, right):
        self.assertLess(abs(float(left-right)), 1e-48)

    def assert_matrix_near(self, left, right):
        self.assertEqual((left.nrows(), left.ncols()), (right.nrows(), right.ncols()))
        for i in range(left.nrows()):
            for j in range(left.ncols()):
                self.assert_near(left[i, j], right[i, j])

    def test_exhaustive_two_step_law(self):
        z = Q(3, 5)
        actual_local = [arb_matrix(physical_exact(self.images, self.columns, j, z)) for j in range(2)]
        composed = kernel.convolve(actual_local)
        for occupancy, matrix in enumerate(composed):
            exact = macro_direct(self.images, self.columns, occupancy, z)
            for q in range(4):
                for r in range(4):
                    self.assertGreaterEqual(endpoint(matrix[q, r]), exact[q][r])
                    self.assert_near(matrix[q, r], kernel.aq(exact[q][r]))

    def test_full_physical_envelopes_dominate_exact_macro_law(self):
        z = Q(3, 5)
        for density, physical in self.physical.items():
            local = sparse_kernel.outward_at_z(physical, kernel.aq(z))
            if density == 'capped':
                local = kernel.kernel_birth_density.refine_local(physical, local, kernel.aq(z), Q(1, 3))
            macro = kernel.convolve(local)
            levels = list(map(int, physical['birth_class_levels']))
            for occupancy, matrix in enumerate(macro):
                exact = macro_direct(self.images, self.columns, occupancy, z)
                for q in range(4):
                    source = 0 if q == 0 else 3+levels.index(self.images[q].bit_count())
                    self.assertGreaterEqual(endpoint(matrix[source, 0]), exact[q][0])
                    self.assertEqual(matrix[source, 1], 0)
                    uniform_cap = endpoint(matrix[source, 2])/3
                    for target, level in enumerate(levels, 3):
                        residual = sum(max(Q(0), exact[q][r]-uniform_cap)
                            for r in range(1, 4) if self.images[r].bit_count() == level)
                        self.assertGreaterEqual(endpoint(matrix[source, target]), residual)

    def test_chronology_and_exact_hypergeometric_split(self):
        # Different half families make an accidental reversal observable.
        left = [arb_mat([[1, j+1], [2*j, 1]]) for j in range(3)]
        right = [arb_mat([[2, j], [j+1, 3]]) for j in range(3)]
        macro = kernel.convolve(left, right)
        for j in range(5):
            direct = arb_mat(2, 2)
            for support in combinations(range(4), j):
                a = sum(p < 2 for p in support)
                direct += left[a]*right[j-a]/comb(4, j)
            self.assert_matrix_near(macro[j], direct)
        self.assertNotEqual(macro[0], right[0]*left[0])
        self.assert_matrix_near(macro[0], left[0]*right[0])
        self.assert_matrix_near(macro[4], left[2]*right[2])

    def test_iid_binomial_mixture_and_normalization(self):
        left = [arb_mat([[1, j+1], [2*j, 1]]) for j in range(3)]
        right = [arb_mat([[2, j], [j+1, 3]]) for j in range(3)]
        macro = kernel.convolve(left, right)
        def mixture(family, p):
            W = len(family)-1
            result = arb_mat(2, 2)
            for j, matrix in enumerate(family):
                result += kernel.aq(comb(W, j)*p**j*(1-p)**(W-j))*matrix
            return result
        for p in (Q(0), Q(1, 5), Q(3, 4), Q(1)):
            self.assert_matrix_near(mixture(macro, p), mixture(left, p)*mixture(right, p))
        for W in (1, 2, 16):
            normalized = kernel.convolve([arb_mat([[1]]) for _ in range(W+1)])
            self.assertTrue(all(endpoint(m[0, 0]) >= 1 for m in normalized))
            self.assertTrue(all(abs(float(m[0, 0]-1)) < 1e-48 for m in normalized))

    def test_zero_start_no_reset_no_flush(self):
        wrapper = kernel.wrap(self.physical['classes'])
        self.assertEqual(wrapper['initial_state'], 'zero')
        self.assertFalse(wrapper['flush'])
        self.assertEqual(wrapper['state_continuity'], 'retained_between_halves')
        z = Q(3, 5)
        exact = macro_direct(self.images, self.columns, 1, z)
        reset_moment = sum(z**x.bit_count() for x in range(1, 16))/15
        self.assertLess(sum(exact[0]), reset_moment)
        self.assertGreater(sum(exact[0][1:]), 0)
        quiet = macro_direct(self.images, self.columns, 0, z)
        self.assertEqual(quiet[0], [Q(1), Q(0), Q(0), Q(0)])
        self.assertTrue(all(quiet[q][0] == 0 for q in range(1, 4)))

    def test_public_iid_and_occupancy_routes(self):
        tilt, activity = Q(1, 4), Q(1, 5)
        probabilities = sc.probabilities(activity)
        for physical in self.physical.values():
            wrapper = kernel.wrap(physical)
            actual = kernel.local_operators(wrapper, tilt, activity)
            self.assertEqual(len(actual), 3)
            matrix = kernel_maps.outward(physical, probabilities, tilt)
            self.assert_matrix_near(kernel.outward(wrapper, probabilities, tilt), matrix*matrix)
            floating = kernel.floating(wrapper, np.array(list(map(float, probabilities))), float(tilt))
            physical_float = kernel_maps.floating(physical, np.array(list(map(float, probabilities))), float(tilt))
            np.testing.assert_array_equal(floating, physical_float@physical_float)

    def test_scope_geometry_and_map_tamper_rejected(self):
        wrapper = kernel.wrap(self.physical['classes'])
        changes = dict(schema='legacy', physical_steps=1, distribution='uniform_gl',
            physical_windows=2, macro_windows=3, windows=3, physical_step_bits=8,
            macro_step_bits=16, state_continuity='reset', initial_state='arbitrary',
            physical_updates_independent=False, flush=True, birth_density='capped',
            map_sha256='0'*64, columns=self.columns, updates=8)
        for key, value in changes.items():
            with self.subTest(key=key), self.assertRaises((ValueError, ArithmeticError)):
                kernel.authenticate(dict(wrapper, **{key: value}))
        changed = deepcopy(wrapper)
        changed['physical_data']['columns'][0] ^= 1
        with self.assertRaises(ArithmeticError):
            kernel.authenticate(changed)
        for family in ([], [arb_mat([[1]])], [arb_mat([[1]]), arb_mat([[1, 0], [0, 1]])],
                       [arb_mat([[1]]), arb_mat([[-1]])], [arb_mat([[1]]), arb_mat([[arb('nan')]])]):
            with self.assertRaises(ValueError):
                kernel.convolve(family)


class SelectedMapTests(unittest.TestCase):
    def test_fresh_selected_map_enumeration_and_untrusted_claims(self):
        images, columns, record = kernel._source_maps(kernel.SELECTED_MAP)
        self.assertEqual(len(images), 65536)
        self.assertEqual(len(columns), 64)
        self.assertEqual(record['expansion_rank'], 16)
        self.assertEqual(record['feedback_rank'], 16)
        self.assertTrue(record['feedback_times_expansion_zero'])
        self.assertEqual(record['packet_ranks'], [4]*16)
        self.assertEqual(record['minimum_expansion_weight'], 16)
        self.assertEqual(record['minimum_feedback_kernel_weight'], 6)
        self.assertEqual(record['expansion_spectrum'],
            {'0': 1, '16': 20, '24': 4640, '28': 13824, '32': 28566,
             '36': 13824, '40': 4640, '48': 20, '64': 1})
        self.assertEqual(record['feedback_kernel_spectrum'][6], '1984')
        self.assertFalse(record['whole_code_certificate'])
        source = json.loads(kernel.SELECTED_MAP.read_text())
        source.update(a_counts={'0': 9}, kernel_counts={'4': 100}, dA=1, dC=1)
        with TemporaryDirectory() as directory:
            path = Path(directory)/'map.json'
            path.write_text(json.dumps(source))
            _, _, fresh = kernel._source_maps(path)
            self.assertEqual(fresh['expansion_spectrum'], record['expansion_spectrum'])
            self.assertEqual(fresh['feedback_kernel_spectrum'], record['feedback_kernel_spectrum'])
            source['columns'][0] ^= 1
            path.write_text(json.dumps(source))
            with self.assertRaises(ArithmeticError):
                kernel._source_maps(path)

    def test_production_wrapper_has_no_fabricated_macro_map(self):
        wrapper, record = kernel.prepare()
        self.assertEqual(wrapper['windows'], 32)
        self.assertEqual(wrapper['physical_windows'], 16)
        self.assertEqual(len(wrapper['physical_data']['columns']), 64)
        self.assertNotIn('columns', wrapper)
        self.assertNotIn('map_images', wrapper)
        self.assertEqual(wrapper['physical_step_bits'], 64)
        self.assertEqual(wrapper['macro_step_bits'], 128)
        self.assertEqual(wrapper['map_sha256'], record['map_sha256'])
        self.assertIs(kernel.authenticate(wrapper), wrapper['physical_data'])


if __name__ == '__main__':
    unittest.main()
