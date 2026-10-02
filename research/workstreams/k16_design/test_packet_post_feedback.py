"""Small exact GL(2,2) and return-enumerator checks; no full map census."""
from fractions import Fraction as Q
from itertools import product
from math import comb
import unittest
from unittest.mock import patch

from flint import arb, ctx, fmpq, fmpq_mat
import packet_post_feedback as post


def aq(value):
    value = Q(value)
    return arb(value.numerator)/value.denominator


def endpoint(value):
    return value.upper().fmpq()


def maps(windows=2, *, constant=True):
    width = 4*windows
    rows = ((1 << width)-1, 15) if constant else (3, 12)
    images = tuple((rows[0] if state & 1 else 0) ^ (rows[1] if state & 2 else 0)
                   for state in range(4))
    columns = [sum(((row >> i) & 1) << j for j, row in enumerate(rows)) for i in range(width)]
    return post.q1.kernel_t64.kernel_maps.prepare_maps(images, columns, bits=2,
        distribution='uniform_gl', birth_density='classes')


def gl2():
    return [(a, b) for a in range(1, 4) for b in range(1, 4) if a != b]


def apply(matrix, value):
    return (matrix[0] if value & 1 else 0) ^ (matrix[1] if value & 2 else 0)


def feedback(data, word):
    result = 0
    for bit, column in enumerate(data['columns']):
        if word >> bit & 1:
            result ^= column
    return result


def literal(data, occupancy, z=Q(3, 4), *, post_transvection=False):
    W = data['windows']
    denominator = comb(W, occupancy)*15**occupancy
    result = [[Q(0)]*4 for _ in range(4)]
    sampled = gl2()
    transvections = [(u, v) for u in range(1, 4) for v in range(4)
                    if not (u & v).bit_count() & 1]
    for word in range(1 << (4*W)):
        if sum(bool((word >> (4*p)) & 15) for p in range(W)) != occupancy:
            continue
        c = feedback(data, word)
        for state in range(4):
            mass = z**(word ^ data['map_images'][state]).bit_count()/denominator
            for matrix in sampled:
                if post_transvection:
                    before = apply(matrix, state) ^ c
                    for u, v in transvections:
                        target = before ^ (u if (v & before).bit_count() & 1 else 0)
                        result[state][target] += mass/(len(sampled)*len(transvections))
                else:
                    target = apply(matrix, state ^ c)
                    result[state][target] += mass/len(sampled)
    return result


def collapsed(matrix):
    return [[matrix[0][0], sum(matrix[0][1:])],
            [sum(matrix[a][0] for a in range(1, 4))/3,
             sum(sum(matrix[a][1:]) for a in range(1, 4))/3]]


class PostFeedbackTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.data = maps()
        cls.model = post.PostFeedbackModel(cls.data, exact_through=2)

    def setUp(self):
        self.precision = ctx.prec
        ctx.prec = 192

    def tearDown(self):
        ctx.prec = self.precision

    def test_exact_return_counts_through_three(self):
        data = maps(3)
        counts = post.return_histograms(data, exact_through=3)
        for j in range(4):
            returns, zeros = [0]*13, [0]*13
            for x in range(4096):
                if sum(bool(x >> (4*p) & 15) for p in range(3)) != j:
                    continue
                c = feedback(data, x)
                transformed = x ^ data['map_images'][c]
                (returns if c else zeros)[transformed.bit_count()] += 1
            self.assertEqual(counts['return_counts'][j], tuple(returns))
            self.assertEqual(counts['zero_counts'][j], tuple(zeros))
            self.assertEqual(sum(returns)+sum(zeros), comb(3, j)*15**j)

    def test_exact_two_state_law_matches_all_gl2_maps(self):
        for precision in (192, 256):
            ctx.prec = precision
            family = self.model.physical_operators_at_z(aq(Q(3, 4)))
            for occupancy, bound in enumerate(family):
                full = literal(self.data, occupancy)
                exact = collapsed(full)
                for i in range(2):
                    for j in range(2):
                        truth = fmpq(exact[i][j].numerator, exact[i][j].denominator)
                        self.assertGreaterEqual(endpoint(bound[i, j]), truth)
                        self.assertLess(abs(bound[i, j]-aq(exact[i][j])), arb(2)**(-precision//2))
                # Each weighted outgoing nonzero coordinate is identical,
                # even from an individually fixed, nonuniform entering state.
                for row in full:
                    self.assertEqual(row[1], row[2])
                    self.assertEqual(row[2], row[3])

    def test_weighted_invariant_composes_for_many_steps(self):
        occupancy = 1
        exact = literal(self.data, occupancy)
        reduced = collapsed(exact)
        full = fmpq_mat([[fmpq(v.numerator, v.denominator) for v in row] for row in exact])
        small = fmpq_mat([[fmpq(v.numerator, v.denominator) for v in row] for row in reduced])
        for steps in (1, 2, 3, 7):
            lhs, rhs = full**steps, small**steps
            self.assertEqual(sum(lhs[0, j] for j in range(4)), sum(rhs[0, j] for j in range(2)))

    def test_conservative_return_caps_dominate_literal_entries(self):
        for constant in (False, True):
            data = maps(constant=constant)
            model = post.PostFeedbackModel(data, exact_through=0)
            for z in (Q(1), Q(3, 4), Q(1, 3)):
                matrices, details = model.physical_operators_at_z(aq(z), diagnostics=True)
                for j, matrix in enumerate(matrices):
                    exact = collapsed(literal(data, j, z))
                    for a in range(2):
                        for b in range(2):
                            truth = fmpq(exact[a][b].numerator, exact[a][b].denominator)
                            self.assertGreaterEqual(endpoint(matrix[a, b]), truth)
                    g = fmpq(exact[1][0].numerator, exact[1][0].denominator)
                    self.assertLessEqual(details[j]['G_lower'].fmpq(), g)
                    self.assertGreaterEqual(details[j]['G_upper'].fmpq(), g)

    def test_old_H_over_m_cap_does_not_transfer(self):
        x = 1 ^ self.data['map_images'][feedback(self.data, 1)]
        c = feedback(self.data, x)
        z = Q(3, 4)
        G = z**(x ^ self.data['map_images'][c]).bit_count()/3
        H = sum(z**(x ^ a).bit_count() for a in self.data['map_images'][1:])/3
        self.assertGreater(G, H/3)

    def test_zero_feedback_branch_and_macro_chronology(self):
        local = self.model.physical_operators_at_z(aq(Q(3, 4)))
        self.assertEqual(local[0][0, 0], 1)
        self.assertEqual(local[0][0, 1], 0)
        self.assertEqual(local[0][1, 0], 0)
        macro = post.q1.kernel_t64.convolve(local)
        expected = (local[0]*local[2] + 4*local[1]*local[1] + local[2]*local[0])/6
        for i in range(2):
            for j in range(2):
                self.assertGreaterEqual(macro[2][i, j], expected[i, j].lower())
                self.assertLess(abs(macro[2][i, j]-expected[i, j]), arb(2)**-100)

    def test_extra_transvection_positive_composition(self):
        data = post.q1.kernel_t64.wrap(self.data)
        with patch.object(post.sparse, 'validated', return_value=2):
            model = post.PostTransvectionModel(data, {}, updates=1)
        # A rational z permits exact comparison with the full sampled kernel.
        z = aq(Q(3, 4))
        old = post.q1.kernel_t64.sparse_kernel.outward_at_z(self.data, z)
        size = old[0].nrows()
        refresh = post.arb_mat(size, size)
        refresh[0, 0], refresh[2, 2] = 1, 1
        for i in (1, *range(3, size)):
            refresh[i, i], refresh[i, 2] = arb(1)/2, arb(1)/2
        operator = post.q1.rounded(old[1]*refresh)
        exact = literal(self.data, 1, post_transvection=True)
        full = fmpq_mat([[fmpq(v.numerator, v.denominator) for v in row] for row in exact])
        for steps in (1, 2, 4):
            truth = full**steps
            bound = operator**steps
            self.assertGreaterEqual(endpoint(sum(bound[0, j] for j in range(size))),
                                    sum(truth[0, j] for j in range(4)))
        self.assertEqual(model.metadata['lazy_weight'], '1/2')
        self.assertFalse(model.metadata['whole_code_certificate'])
        # Exercise the public method as well as the exact rational composition.
        tilt = Q(1, 100)
        physical = model.physical_operators(tilt)
        base = post.q1.kernel_t64.sparse_kernel.outward_at_z(self.data, (-post.aq(tilt)).exp())
        self.assertEqual(physical, [post.q1.rounded(operator*refresh) for operator in base])
        self.assertEqual(model.local_operators(tilt), post.q1.kernel_t64.convolve(physical))

    def test_metadata_does_not_relabel_original_certificate(self):
        record = self.model.metadata
        self.assertEqual(record['recurrence'], 'y=x+A*a; next_a=M*(a+C*x)')
        self.assertFalse(record['whole_code_certificate'])
        self.assertTrue(record['maps_shared_with_original_only'])
        self.assertEqual(record['matrix_coordinates'], ['zero', 'uniform_nonzero'])

    def test_invalid_maps_cutoffs_and_weights(self):
        for cutoff in (-1, 3, True):
            with self.assertRaises(ValueError):
                post.return_histograms(self.data, exact_through=cutoff)
        for z in (0, -1, 2):
            with self.assertRaises(ValueError):
                self.model.physical_operators_at_z(arb(z))
        changed = dict(self.data, columns=[2, *self.data['columns'][1:]])
        with self.assertRaises((ValueError, ArithmeticError)):
            post.PostFeedbackModel(changed, exact_through=1)

    def test_receipt_authentication_checks_direct_constructor_records(self):
        model = post.PostFeedbackModel(self.data, exact_through=1, source_map_record={'wrong': 'map'})
        with patch.object(post.sparse, 'validated', side_effect=ValueError('wrong source map')) as validate:
            with self.assertRaises(ValueError):
                post._authenticate_model(model)
            validate.assert_called_once()
        model.counts = dict(model.counts, exact_through=0)
        with self.assertRaises(ValueError):
            post._authenticate_model(model)


if __name__ == '__main__':
    unittest.main()
