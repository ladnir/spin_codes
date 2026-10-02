"""Small exact checks of coset dominance and outward transfer arithmetic."""
from fractions import Fraction as Q
from itertools import combinations
from math import comb
import unittest

from flint import arb, ctx
import packet_outer_block_transfer as transfer


def span(rows):
    values = {0}
    for row in rows:
        values |= {value ^ row for value in values}
    return values


def moment(values, shift, z):
    return sum(z**((value ^ shift).bit_count()) for value in values)/len(values)


class OuterTransferTests(unittest.TestCase):
    def test_coset_dominance_exhaustive_small_subspaces(self):
        for n in range(1, 5):
            spaces = {tuple(sorted(span(rows))) for size in range(n+1)
                      for rows in combinations(range(1, 1 << n), size)}
            for values in spaces:
                for z in (Q(1, 4), Q(1, 2), Q(3, 4), Q(1)):
                    zero = moment(values, 0, z)
                    for shift in range(1 << n):
                        self.assertLessEqual(moment(values, shift, z), zero)

    def test_nonuniform_nz_input_does_not_have_coset_property(self):
        self.assertGreater(moment({1}, 1, Q(1, 2)), moment({1}, 0, Q(1, 2)))

    def test_continuous_state_block_induction_exact_toy(self):
        # Two output bits and two state bits per step; no state reset.
        def block(state):
            outcomes = []
            for x in (0, 1):
                for y in (0, 2):
                    emitted = x ^ state
                    middle = ((state << 1) | (state >> 1)) & 3
                    middle ^= x
                    emitted |= (y ^ middle) << 2
                    outgoing = (((middle << 1) | (middle >> 1)) & 3) ^ y
                    outcomes.append((outgoing, emitted.bit_count()))
            return outcomes
        z = Q(1, 2)
        zero = sum(z**weight for _, weight in block(0))/4
        for entering in range(4):
            self.assertLessEqual(sum(z**weight for _, weight in block(entering))/4, zero)
        mass = {0: Q(1)}
        for count in range(1, 5):
            next_mass = {}
            for entering, probability in mass.items():
                for outgoing, weight in block(entering):
                    next_mass[outgoing] = next_mass.get(outgoing, Q(0)) + probability*z**weight/4
            mass = next_mass
            self.assertLessEqual(sum(mass.values()), zero**count)

    def test_beta_ratio_exact(self):
        for m in (2, 4):
            actual = Q(1 << (256*m), ((1 << 32)-1)**(4*m)) / transfer.BASE_BETA**m
            self.assertEqual(transfer.beta_ratio(m), actual)
            self.assertLess(actual, 1)

    def test_cutoff_difference(self):
        for m in (2, 4):
            self.assertEqual((524288*m)//10-m*(524288//10), m-1)

    def test_outward_endpoint_against_higher_precision(self):
        previous = ctx.prec
        try:
            for m in (2, 4):
                for q in (3, 89, 385, 2048):
                    ctx.prec = 192
                    endpoint = transfer.transfer_endpoint([17, -50], q=q, tilt='17/200', multiplier=m)
                    ctx.prec = 768
                    truth = ((arb(17)/200*(m-1)).exp() * transfer._aq(transfer.beta_ratio(m))**q *
                             (arb(17)*arb(2)**-50)**m / arb(comb(2048, q))**(m-1))
                    bound = arb(endpoint[0])*arb(2)**endpoint[1]
                    self.assertGreater(bound, truth)
        finally:
            ctx.prec = previous

    def test_invalid_inputs(self):
        for m in (1, 3, True):
            with self.assertRaises(ValueError):
                transfer.beta_ratio(m)
        for q in (0, 1, 2, 2049, True):
            with self.assertRaises(ValueError):
                transfer.transfer_endpoint([1, -2], q=q, tilt='1/10', multiplier=2)
        for endpoint in ([0, 0], [-1, 0], [1.0, 0], [1, 10000001]):
            with self.assertRaises(ValueError):
                transfer.transfer_endpoint(endpoint, q=3, tilt='1/10', multiplier=2)


if __name__ == '__main__':
    unittest.main()
