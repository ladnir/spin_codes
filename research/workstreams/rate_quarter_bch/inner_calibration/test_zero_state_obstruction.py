"""Exact regression checks for deterministic and first-moment obstructions."""
import json
import math
import random
import unittest

import zero_state_rank_ceiling as rank_bound
import zero_state_obstruction as lower


class ObstructionTests(unittest.TestCase):
    def test_rank_ceiling_targets(self):
        expected = {18:589888,19:622656,20:655424,21:688192,22:720960,24:786496,25:819264}
        for s,bound in expected.items():
            row = rank_bound.ceiling(128,32,32768,256,s)
            self.assertEqual(row['minimum_distance_upper'],bound)
            self.assertEqual(row['zero_state_subspace_dimension_lower'],32)
        self.assertIsNone(rank_bound.ceiling(128,32,32768,256,64))

    def test_rank_bound_against_toy_setups(self):
        # Full-rank [4,2] outer: (a,b,a,b). Linear parity syndrome in each
        # 8-bit epoch. Enumerate the entire restricted message subspace.
        geometry = rank_bound.ceiling(4,2,16,8,1)
        q = geometry['restricted_outer_rows']
        rng = random.Random(719)
        for _ in range(12):
            permutation = list(range(64)); rng.shuffle(permutation)
            generators = [sum(1<<permutation[4*row+j] for j in (bit,bit+2))
                          for row in range(q) for bit in range(2)]
            word = 0; accepted = []
            for index in range(1,1<<len(generators)):
                word ^= generators[(index&-index).bit_length()-1]
                if all(((word>>(8*e))&255).bit_count()%2==0 for e in range(8)):
                    accepted.append(word)
            self.assertGreaterEqual(len(accepted),(1<<geometry['zero_state_subspace_dimension_lower'])-1)
            self.assertLessEqual(min(w.bit_count() for w in accepted),geometry['minimum_distance_upper'])

    def test_coefficient_term_is_positive_polynomial_term(self):
        kernel = [1,0,3,0,1]
        coefficients = [1]
        for _ in range(3):
            out = [0]*(len(coefficients)+len(kernel)-1)
            for i,x in enumerate(coefficients):
                for j,y in enumerate(kernel): out[i+j]+=x*y
            coefficients = out
        for j in range(0,13,2):
            self.assertLessEqual(lower.coefficient_term(kernel,3,j),coefficients[j])

    def test_integer_chernoff_against_exact_small_tail(self):
        lower.exact_tail_check(100,1,2,10,40,True)
        lower.exact_tail_check(100,1,2,90,40,False)
        numerator = sum(math.comb(100,j) for j in range(11))
        self.assertLessEqual(numerator<<40,1<<100)
        with self.assertRaises(AssertionError):
            lower.exact_tail_check(100,1,2,10,99,True)

    def test_retained_integer_certificates(self):
        folder = lower.check.calibration.HERE
        payload = json.loads((folder/'ZERO_STATE_OBSTRUCTION.json').read_text())
        counts = lower.check.calibration.smaller_outer.spectrum()
        for name,digest in payload['source_sha256'].items():
            self.assertEqual(lower.check.fixed.sha(lower.check.fixed.ROOT/name),digest,name)
        for row in payload['results']:
            record = json.loads((folder/'maps'/f"{row['tag']}.json").read_text())
            replay = dict(tag=row['tag'],**lower.certify(record,counts,lower.check.Fraction(row['distance_target'])))
            self.assertEqual(replay,row)
            self.assertGreater(row['expected_bad_messages_lower_power_of_two'],0)


if __name__=='__main__': unittest.main()
