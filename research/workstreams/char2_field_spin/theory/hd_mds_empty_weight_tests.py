"""Exact GF4 refresh/emission controls; no regional provider or gate."""
from collections import Counter
from fractions import Fraction as F
from itertools import product
from unittest import TestCase, main

from hd_blockwise_maps import mul
from hd_grouped_quiet import weight_transition
from hd_mds_empty_weight import emission_by_weight, empty_kernel, flat_cap_from_kernel, local_positive_counts

EMISSION = ((1, 0), (1, 1), (1, 2), (1, 3))
REFRESH = ((3, 2), (2, 3))


def apply(matrix, vector):
    return tuple(mul(row[0], vector[0], 4) ^ mul(row[1], vector[1], 4) for row in matrix)


def refresh_counts(vector):
    result = Counter()
    for before0, before1, after0, after1 in product(range(1, 4), repeat=4):
        middle = apply(REFRESH, (mul(before0, vector[0], 4), mul(before1, vector[1], 4)))
        result[mul(after0, middle[0], 4), mul(after1, middle[1], 4)] += 1
    return result


class EmptyWeightTests(TestCase):
    def test_mds_transition_every_fixed_input_and_target(self):
        transition = weight_transition(4, 2)
        for vector in product(range(4), repeat=2):
            counts = refresh_counts(vector)
            a = sum(bool(x) for x in vector)
            if not a:
                self.assertEqual(counts, {(0, 0): 81})
                continue
            self.assertNotIn((0, 0), counts)
            for target in product(range(4), repeat=2):
                b = sum(bool(x) for x in target)
                if b:
                    cardinality = (2 if b == 1 else 1)*3**b
                    self.assertEqual(F(counts[target], 81), transition[a-1][b-1]/cardinality)

    def test_actual_basis_emission_counts(self):
        counts = local_positive_counts(EMISSION, 4)
        expected = [[0, 0] for _ in range(3)]
        for vector in product(range(4), repeat=2):
            weight = sum(bool(x) for x in vector)
            if weight:
                output = apply(EMISSION, vector)
                expected[weight][0] += 1
                expected[weight][1] += sum(x == 0 for x in output)
        self.assertEqual(counts, tuple(map(tuple, expected)))
        means = emission_by_weight(F(3, 5), q=4, blocks=1, emission=EMISSION)
        for weight in range(3):
            vectors = [x for x in product(range(4), repeat=2) if sum(bool(v) for v in x) == weight]
            direct = sum((F(3, 5)**sum(bool(v) for v in apply(EMISSION, x)) for x in vectors), F(0))/len(vectors)
            self.assertEqual(means[weight], direct)

    def test_marked_kernel_factorization_and_weighted_history(self):
        z = F(1, 2)
        kernel, means, transition = empty_kernel(z, q=4, blocks=1, emission=EMISSION)
        for a in range(3):
            vectors = [x for x in product(range(4), repeat=2) if sum(bool(v) for v in x) == a]
            direct = [F(0)]*3
            for vector in vectors:
                mark = z**sum(bool(v) for v in apply(EMISSION, vector))
                for target, count in refresh_counts(vector).items():
                    direct[sum(bool(v) for v in target)] += mark*F(count, 81*len(vectors))
            self.assertEqual(tuple(direct), kernel[a])
        # Prior weighted histories need not give equal mass to input states.
        # Fresh refresh still makes the resulting density equal within each
        # output-weight class. This checks the actual conditional interface.
        output = Counter()
        for vector in product(range(4), repeat=2):
            if not any(vector):
                continue
            weight = F(1+2*vector[0]+vector[1])*z**sum(bool(v) for v in apply(EMISSION, vector))
            for target, count in refresh_counts(vector).items():
                output[target] += weight*F(count, 81)
        for b in (1, 2):
            values = {mass for target, mass in output.items() if sum(bool(v) for v in target) == b}
            self.assertEqual(len(values), 1)

    def test_mark_one_and_block_convolution(self):
        means = emission_by_weight(1, q=4, blocks=2, emission=EMISSION)
        self.assertEqual(means, (F(1),)*5)
        means = emission_by_weight(F(1, 2), q=4, blocks=2, emission=EMISSION)
        totals, counts = [F(0)]*5, [0]*5
        for vector in product(range(4), repeat=4):
            a = sum(bool(x) for x in vector)
            output = apply(EMISSION, vector[:2])+apply(EMISSION, vector[2:])
            totals[a] += F(1, 2)**sum(bool(x) for x in output)
            counts[a] += 1
        self.assertEqual(means, tuple(total/count for total, count in zip(totals, counts)))

    def test_inverse_scalar_is_pointwise_density_not_total_mass(self):
        z = F(3, 5)
        kernel, means, transition = empty_kernel(z, q=4, blocks=1, emission=EMISSION)
        scalar, targets, pi = flat_cap_from_kernel(means, transition, 4)
        direct = Counter()
        for vector in product(range(4), repeat=2):
            if any(vector):
                mark = z**sum(bool(x) for x in apply(EMISSION, vector))
                for target, count in refresh_counts(vector).items():
                    direct[target] += mark*F(count, 81*15)
        for target, mass in direct.items():
            b = sum(bool(x) for x in target)
            self.assertEqual(15*mass, targets[b-1])
            self.assertLessEqual(15*mass, scalar)
        self.assertEqual(scalar, max(15*mass for mass in direct.values()))


if __name__ == '__main__':
    main()
