"""Independent finite checks for the affine safe-count diagnostics."""
from itertools import product
import unittest

import numpy as np

from hd_affine_cut import (ROWS, STRIDES, candidate_sets, grouped_candidate_sets,
                           grouped_safe_counts, hillclimb, safe_counts,
                           spectral_safe_upper, stationary_pressure)
from hd_two_state_cut_screen import stationary_statistics


class AffineCutTests(unittest.TestCase):
    def test_safe_counts_against_literal_global_inverse_support(self):
        for rows, strides, offsets in ((4, (1, 3), (0, 0)),
                                       (8, (1, 3, 5), (0, 0, 0)),
                                       (8, (1, 3, 5), (2, 7, 1))):
            labels = [(stride*j+offset) % rows
                      for stride, offset in zip(strides, offsets) for j in range(rows)]
            for bits in product((0, 1), repeat=rows):
                expected = [0]*len(strides)
                for j, row in enumerate(labels):
                    safe = bits[row] and (j+1 == len(labels) or bits[labels[j+1]])
                    expected[j//rows] += int(safe)
                self.assertEqual(safe_counts(bits, strides, offsets), tuple(expected))

    def test_zero_offsets_only_have_one_terminal_correction(self):
        rng = np.random.default_rng(619)
        for _ in range(8):
            selected = rng.integers(0, 2, size=ROWS).astype(bool)
            counts = safe_counts(selected)
            cyclic = tuple(int(np.count_nonzero(selected & np.roll(selected, -s)))
                           for s in STRIDES)
            self.assertEqual(counts[:-1], cyclic[:-1])
            self.assertEqual(counts[-1]-cyclic[-1],
                             int(selected[-STRIDES[-1]] and not selected[0]))

    def test_uniform_spectral_cardinality_range(self):
        for a in range(1, 8568):
            self.assertLess(spectral_safe_upper(a), 32*a)
        self.assertGreaterEqual(spectral_safe_upper(8568), 32*8568)
        self.assertEqual(spectral_safe_upper(8192), 258081)

    def test_pressure_against_independent_moment_chain(self):
        for safe, num, den in ((ROWS//2, 1, 1), (ROWS//4, 1, 3),
                               (ROWS//8, 2, 3), (3*ROWS//4, 2, 1)):
            expected = stationary_statistics(safe/ROWS, num, den, 2)["score"]
            self.assertAlmostEqual(stationary_pressure(safe, num, den), expected, places=11)

    def test_bounded_hillclimb_preserves_size_and_improves_edges(self):
        rng = np.random.default_rng(421)
        selected = np.zeros(ROWS, dtype=bool)
        selected[rng.choice(ROWS, ROWS//3, replace=False)] = True
        before = sum(safe_counts(selected))
        after, moves = hillclimb(selected, 12)
        self.assertEqual(int(after.sum()), int(selected.sum()))
        self.assertEqual(moves, 12)
        self.assertGreaterEqual(sum(safe_counts(after)), before)

    def test_candidate_counts_obey_uniform_spectral_receipt(self):
        candidates = list(candidate_sets())
        self.assertEqual(len(candidates), 92)
        for name, selected in candidates:
            self.assertLessEqual(sum(safe_counts(selected)),
                                 spectral_safe_upper(int(selected.sum())), name)

    def test_grouped_counts_against_literal_global_inverse_support(self):
        for rows, group, strides in ((8, 2, (1, 3)), (8, 4, (1, 3)),
                                     (12, 4, (1, 2))):
            size = rows//group
            labels = [group*((stride*(j//group)) % size)+j % group
                      for stride in strides for j in range(rows)]
            for bits in product((0, 1), repeat=rows):
                expected = [0]*len(strides)
                for j, row in enumerate(labels):
                    safe = bits[row] and (j+1 == len(labels) or bits[labels[j+1]])
                    expected[j//rows] += int(safe)
                self.assertEqual(grouped_safe_counts(bits, group, strides), tuple(expected))

    def test_grouped_full_block_identity_and_menu(self):
        rng = np.random.default_rng(212)
        for group, count in ((2, 69), (4, 129)):
            bits = rng.integers(0, 2, size=ROWS//group).astype(bool)
            actual = grouped_safe_counts(np.repeat(bits, group), group)
            expected = [(group-1)*int(bits.sum())+
                        int(np.count_nonzero(bits & np.roll(bits, -stride)))
                        for stride in STRIDES]
            expected[-1] += int(bits[-STRIDES[-1] % len(bits)] and not bits[0])
            self.assertEqual(actual, tuple(expected))
            self.assertEqual(len(list(grouped_candidate_sets(group))), count)

    def test_retained_group4_borderline_refinements(self):
        targets = {"group4-fourier-1599-819", "group4-fourier-87-819",
                   "group4-fourier-619-819"}
        seen = set()
        for name, selected in grouped_candidate_sets(4):
            if name not in targets:
                continue
            seen.add(name)
            a = int(selected.sum())
            capacity = a/(2*ROWS)
            pressure = sum(stationary_pressure(z, 1, 12)
                           for z in grouped_safe_counts(selected, 4))/64
            gap = pressure+(0.4-capacity)/12-capacity
            self.assertLess(gap, -0.0026)
        self.assertEqual(seen, targets)


if __name__ == "__main__":
    unittest.main()
