"""Exhaustive tiny-field checks for the new refresh-order bounds."""
from math import comb
import unittest

import numpy as np

import refresh_gate as gate


def tiny():
    # GF4^2 state, four two-bit packets. CA=0 and all packet pairs full rank.
    mul = lambda a, b: _mul(a, b)
    rows = tuple(sum((1 << i if i < 2 else mul(h, 1 << (i-2))) << (2*h)
                     for h in range(4)) for i in range(4))
    columns = tuple((1 << bit) | (mul(h, 1 << bit) << 2) for h in range(4) for bit in range(2))
    return gate.maps.prepare_maps(rows, columns, packet_bits=2)[0]


def _mul(a, b):
    out = 0
    while b:
        if b & 1:
            out ^= a
        b >>= 1
        a <<= 1
        if a & 4:
            a ^= 7
    return out


class RefreshTests(unittest.TestCase):
    def test_positive_return_census(self):
        data = tiny()
        census = gate.exact_returns(data)
        expected = np.zeros((3, 9), dtype=np.int64)
        for x in range(256):
            j = sum(bool((x >> (2*h)) & 3) for h in range(4))
            syndrome = gate.maps.apply(data['columns'], x)
            if 0 < j <= 2 and syndrome:
                expected[j, (x ^ int(data['images'][syndrome])).bit_count()] += 1
        np.testing.assert_array_equal(census['histograms'], expected)

    def test_shared_refresh_operator_domination(self):
        data, z = tiny(), .67
        weighted, emission, _ = gate.prior.moments(data, z)
        local, _ = gate.operators(data, weighted, emission, z, mode='shared',
            census=gate.exact_returns(data), holder=gate.prepare_holder(data, block_bits=4))
        exact = np.zeros_like(local)
        for x in range(256):
            j = sum(bool((x >> (2*h)) & 3) for h in range(4))
            denominator = comb(4, j)*3**j
            syndrome = gate.maps.apply(data['columns'], x)
            exact[j, 0, int(syndrome != 0)] += z**x.bit_count()/denominator
            for state in range(1, 16):
                destination = int((state ^ syndrome) != 0)
                exact[j, 1, destination] += z**(x ^ int(data['images'][state])).bit_count()/(15*denominator)
        self.assertTrue(np.all(exact <= local+2e-14))
        np.testing.assert_allclose(local[:3], exact[:3], atol=2e-14, rtol=1e-12)

    def test_independent_refresh_operator_domination(self):
        data, z = tiny(), .83
        weighted, emission, _ = gate.prior.moments(data, z)
        local, _ = gate.operators(data, weighted, emission, z, mode='independent')
        exact = np.zeros_like(local)
        for x in range(256):
            j = sum(bool((x >> (2*h)) & 3) for h in range(4))
            denominator = comb(4, j)*3**j
            syndrome = gate.maps.apply(data['columns'], x)
            exact[j, 0, int(syndrome != 0)] += z**x.bit_count()/denominator
            for state in range(1, 16):
                mass = z**(x ^ int(data['images'][state])).bit_count()/(15*denominator)
                returned = (1/15) if syndrome else 0.
                exact[j, 1] += mass*np.array([returned, 1-returned])
        self.assertTrue(np.all(exact <= local+2e-14))

    def test_holder_return_majorant(self):
        data = tiny()
        prepared = gate.prepare_holder(data, block_bits=4)
        for z in (.2, .7, 1.):
            exact = np.zeros(5)
            for x in range(256):
                j = sum(bool((x >> (2*h)) & 3) for h in range(4))
                output = x ^ int(data['images'][gate.maps.apply(data['columns'], x)])
                exact[j] += z**output.bit_count()/(comb(4, j)*3**j)
            upper = gate.holder_returns(prepared, z)
            self.assertTrue(np.all(exact <= upper+2e-14))


if __name__ == '__main__':
    unittest.main()
