"""Exact small-field checks for the outward local comparison arithmetic."""
from fractions import Fraction
from math import comb, nextafter, inf
import unittest

import numpy as np
import local_outward24 as o


def exact_small(data, z):
    b, W, S = data['packet_bits'], data['windows'], 1 << data['bits']
    labels = (1 << b)-1
    images = [o.maps24.apply(data['rows'], a) for a in range(S)]
    weighted = [[Fraction(0) for _ in range(W+1)] for _ in range(S)]
    emission = [[Fraction(0) for _ in range(W+1)] for _ in range(S)]
    for x in range(1 << (b*W)):
        j = sum(bool((x >> (b*h)) & labels) for h in range(W))
        den = comb(W, j)*labels**j
        weighted[o.maps24.apply(data['columns'], x)][j] += z**x.bit_count()/den
        for a in range(S):
            emission[a][j] += z**(x ^ images[a]).bit_count()/den
    local = [[[Fraction(0) for _ in range(W+2)] for _ in range(W+2)] for _ in range(W+1)]
    for j in range(W+1):
        local[j][0][0] = weighted[0][j]
        if j:
            local[j][0][j+1] = sum(row[j] for row in weighted[1:])
        local[j][1][1] = sum(row[j] for row in emission[1:])/(S-1)
        for i in range(1, W+1):
            mass = sum(row[i] for row in weighted[1:])
            local[j][i+1][1] = sum(weighted[a][i]*emission[a][j] for a in range(1, S))/mass
        if j:
            for i in range(1, W+2):
                local[j][i][0] = local[j][i][1]/(S-1)
    return weighted, emission, local


class OutwardLocal(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.data = o.maps24.make_maps(packet_bits=2, modulus=7, windows=4)
        cls.census = o.maps24.census(cls.data, chunk_bits=3)
        cls.z = Fraction(47, 64)
        cls.weighted, cls.emission, cls.exact = exact_small(cls.data, cls.z)

    def test_rational_conversion_checked_exactly(self):
        for value in (Fraction(0), Fraction(1, 3), Fraction(-7, 19), Fraction(1, 1 << 1000), Fraction(1 << 999)):
            lower, upper = o.enclose(value)
            self.assertLessEqual(o.fraction(lower), value)
            self.assertGreaterEqual(o.fraction(upper), value)
            self.assertTrue(lower == upper or nextafter(lower, inf) == upper)

    def test_exact_profile_polynomials(self):
        selected = np.flatnonzero(self.census['histogram'])
        nums, dens = o.exact_profiles(self.data, self.z, selected)
        for a in range(64):
            profile = int(self.census['expansion_indices'][a])
            for j in range(5):
                self.assertEqual(Fraction(nums[profile][j], dens[j]), self.emission[a][j])

    def test_walsh_and_compression_enclose_exact_masses(self):
        selected = np.unique(self.census['character_indices'])
        nums, dens = o.exact_profiles(self.data, self.z, selected, character=True)
        for j in range(5):
            centers, error, diagnostic = o.walsh_enclosure({p: row[j] for p, row in nums.items()},
                dens[j], self.census['character_indices'], 6, block_bits=2)
            self.assertTrue(diagnostic['no_subnormal_or_overflow_by_lattice_and_magnitude'])
            for a in range(64):
                self.assertLessEqual(abs(o.fraction(centers[a])-self.weighted[a][j]), error)
            bounds, zero, _ = o.compress_enclosure(centers, error, self.census['expansion_indices'],
                                                  self.census['histogram'], chunk=7)
            self.assertLessEqual(zero[0], self.weighted[0][j])
            self.assertGreaterEqual(zero[1], self.weighted[0][j])
            for p, (lower, upper) in bounds.items():
                exact = sum(self.weighted[a][j] for a in range(1, 64)
                            if int(self.census['expansion_indices'][a]) == p)
                self.assertLessEqual(lower, exact)
                self.assertGreaterEqual(upper, exact)

    def test_every_local_entry_contains_exact_comparison(self):
        lower, upper, diagnostic = o.local_operators(self.data, self.census, self.z)
        for j in range(5):
            for source in range(6):
                for destination in range(6):
                    exact = self.exact[j][source][destination]
                    self.assertLessEqual(o.fraction(lower[j, source, destination]), exact)
                    self.assertGreaterEqual(o.fraction(upper[j, source, destination]), exact)
        self.assertLess(diagnostic['max_relative_width'], 1e-9)
        self.assertTrue(diagnostic['exact_uniform_means'])

    def test_invalid_weight_and_unjustified_range_rejected(self):
        for z in (0, 1, Fraction(2, 3), Fraction(-1, 4)):
            with self.assertRaises(ValueError):
                o.validate(self.data, z)
        with self.assertRaises(FloatingPointError):
            o.walsh_enclosure({0: 1}, 1 << 1020, np.zeros(64, dtype=np.uint16), 6)


if __name__ == '__main__':
    unittest.main()
