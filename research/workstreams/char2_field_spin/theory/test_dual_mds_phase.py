"""Replayable exact checks for the new dual-MDS phase maps only."""
import itertools
from random import Random
import unittest
from unittest.mock import patch

import dual_mds_phase as dual
from rate_design_coset import mul, power
from region64_reference import (
    IDENTITY, base_maps as ordinary_base_maps, dot, flatten, inverse,
    matmul, matvec, phase_geometry, random_orientation, translation,
    translation_matrix, transpose,
)


class DualMdsPhaseTests(unittest.TestCase):
    def test_translation_preserves_reversal_pairing(self):
        points, _ = phase_geometry()
        self.assertEqual(matmul(dual.J, dual.J), IDENTITY)
        self.assertEqual(transpose(dual.J), dual.J)
        for beta in range(256):
            change = translation_matrix(beta)
            self.assertEqual(matmul(transpose(change), matmul(dual.J, change)), dual.J)
            self.assertEqual(matmul(dual.J, change), matmul(transpose(change), dual.J))
            for u in points:
                self.assertEqual(matvec(change, dual.curve(u)), dual.curve(u ^ beta))
            for coordinate in range(4):
                values = tuple(255 if i == coordinate else 0 for i in range(4))
                self.assertEqual(translation(values, beta), matvec(change, values))
                self.assertEqual(translation(values, beta, True), matvec(transpose(change), values))

    def test_actual_labels_and_pairing(self):
        points, phases = phase_geometry()
        self.assertEqual(phases, (0, 2, 6, 4))
        self.assertTrue(all(power(u, 16) == u for u in points))
        labels = dual.period_labels()
        self.assertEqual(len(labels), 64)
        self.assertEqual(len(set(labels)), 64)
        for a in labels:
            for b in labels:
                self.assertEqual(dot(dual.curve(a), matvec(dual.J, dual.curve(b))), power(a ^ b, 3))

    def test_local_CA_zero_and_direct_global_maps(self):
        ordinary_emission, injection = ordinary_base_maps()
        emission, actual_injection = dual.base_maps()
        self.assertEqual(actual_injection, injection)
        self.assertEqual(emission, matmul(ordinary_emission, dual.J))
        self.assertEqual(matmul(injection, emission), dual.ZERO_MATRIX)
        rng = Random(20261006)
        orientation = random_orientation(rng)
        for beta in range(256):
            change = translation_matrix(beta)
            direct_emission, direct_injection = dual.global_maps(orientation, beta)
            self.assertEqual(direct_emission, matmul(matmul(emission, change), inverse(orientation)))
            self.assertEqual(direct_injection, matmul(matmul(orientation, change), injection))
            self.assertEqual(matmul(direct_injection, direct_emission), dual.ZERO_MATRIX)

    def test_vandermonde_rank_certificate_and_small_subsets(self):
        # Distinctness proves every pairwise difference is nonzero.  The
        # determinant product therefore certifies ALL four-label minors,
        # not just the representative minors explicitly eliminated below.
        labels = dual.period_labels()
        for pair in itertools.combinations(labels, 2):
            self.assertNotEqual(dual.vandermonde_determinant(pair), 0)
        # Four labels from each phase: exhaust all subsets of sizes1..4.
        selected = tuple(labels[16*phase+offset] for phase in range(4) for offset in range(4))
        for size in range(1, 5):
            for subset in itertools.combinations(selected, size):
                matrix = tuple(dual.curve(label)[:size] for label in subset)
                expected = dual.vandermonde_determinant(subset)
                self.assertNotEqual(expected, 0)
                self.assertEqual(dual.determinant(matrix), expected)
        rng = Random(6404)
        orientation = random_orientation(rng)
        emission, injection = dual.period_maps(orientation)
        for _ in range(256):
            indices = tuple(rng.sample(range(64), 4))
            self.assertNotEqual(dual.determinant(tuple(emission[i] for i in indices)), 0)
            self.assertNotEqual(dual.determinant(tuple(tuple(row[i] for i in indices) for row in injection)), 0)
        # Repeated labels must not accidentally receive a full-rank claim.
        repeated = (labels[0], labels[1], labels[2], labels[0])
        self.assertEqual(dual.vandermonde_determinant(repeated), 0)
        self.assertEqual(dual.determinant(tuple(dual.curve(label) for label in repeated)), 0)

    def test_generalized_support_bounds_are_attained(self):
        labels = dual.period_labels()
        emission, _ = dual.period_maps()
        for dimension in range(1, 5):
            # P has 4-d distinct prescribed roots.  The coefficient vectors
            # of P,X*P,...,X^(d-1)*P give a d-dimensional state subspace.
            roots = labels[:4-dimension]
            polynomial = [1]
            for root in roots:
                following = [0]*(len(polynomial)+1)
                for index, coefficient in enumerate(polynomial):
                    following[index] ^= mul(root, coefficient)
                    following[index+1] ^= coefficient
                polynomial = following
            states = []
            for shift in range(dimension):
                coefficients = [0]*shift + polynomial
                coefficients += [0]*(4-len(coefficients))
                states.append(tuple(reversed(coefficients)))
            outputs = tuple(matvec(emission, state) for state in states)
            support = {i for i in range(64) if any(output[i] for output in outputs)}
            self.assertEqual(len(support), 60+dimension)
            self.assertEqual({labels[i] for i in range(64) if i not in support}, set(roots))

    def test_old_selfadjoint_shortcut_is_rejected(self):
        dual_output = (1,)+(0,)*15  # Label zero.
        terminal = (0, 0, 0, 0)
        correct = dual.adjoint_local_step(dual_output, terminal)
        incorrect = dual.local_step(dual_output, terminal)
        self.assertEqual(correct[1], (0, 0, 0, 1))
        self.assertEqual(incorrect[1], (1, 0, 0, 0))
        self.assertNotEqual(correct, incorrect)
        regions = (((0,)*16,),)
        dual_regions = ((dual_output,),)
        initial = (0, 0, 0, 1)
        encoded, _ = dual.encode_local(regions, (IDENTITY,), initial)
        transposed = dual.adjoint_local(dual_regions, (IDENTITY,))
        lhs = dot(flatten(encoded), flatten(dual_regions))
        self.assertEqual(lhs, dot(initial, transposed[1]))
        with patch.object(dual, "adjoint_local_step", side_effect=dual.local_step):
            bad = dual.adjoint_local(dual_regions, (IDENTITY,))
        self.assertNotEqual(lhs, dot(initial, bad[1]))

    def test_local_global_and_extended_adjoint(self):
        rng = Random(20261006)
        sequences = ((1,), (2,), (3,), (4,), (5,), (7,), (1, 2), (3, 4, 5), (7, 1, 6))
        for lengths in sequences:
            for with_diagonal in (False, True):
                for trial in range(2):
                    with self.subTest(lengths=lengths, diagonal=with_diagonal, trial=trial):
                        regions = tuple(tuple(tuple(rng.randrange(256) for _ in range(16))
                                              for _ in range(length)) for length in lengths)
                        values = tuple(tuple(tuple(rng.randrange(256) for _ in range(16))
                                             for _ in range(length)) for length in lengths)
                        orientations = tuple(IDENTITY if trial == 0 else random_orientation(rng) for _ in lengths)
                        diagonals = (tuple(tuple(tuple(rng.randrange(1, 256) for _ in range(16))
                                                 for _ in range(length)) for length in lengths)
                                     if with_diagonal else None)
                        initial = tuple(rng.randrange(1, 256) for _ in range(4))
                        terminal = tuple(rng.randrange(1, 256) for _ in range(4))
                        global_result = dual.encode_global(regions, orientations, initial, diagonals)
                        self.assertEqual(dual.encode_local(regions, orientations, initial, diagonals), global_result)
                        global_dual = dual.adjoint_global(values, orientations, terminal, diagonals)
                        self.assertEqual(dual.adjoint_local(values, orientations, terminal, diagonals), global_dual)
                        output, final = global_result
                        dual_input, dual_initial = global_dual
                        self.assertEqual(dot(flatten(output), flatten(values)) ^ dot(final, terminal),
                                         dot(flatten(regions), flatten(dual_input)) ^ dot(initial, dual_initial))
                        # The actual zero-initial/discarded-final interface.
                        output, _ = dual.encode_local(regions, orientations, diagonals=diagonals)
                        dual_input, _ = dual.adjoint_local(values, orientations, diagonals=diagonals)
                        self.assertEqual(dot(flatten(output), flatten(values)), dot(flatten(regions), flatten(dual_input)))

    def test_invalid_external_state_and_geometry(self):
        region = (((0,)*16,),)
        for bad in ((0, 0, 0), (0, 0, 0, 256)):
            with self.assertRaises(ValueError):
                dual.encode_local(region, (IDENTITY,), initial=bad)
            with self.assertRaises(ValueError):
                dual.adjoint_global(region, (IDENTITY,), terminal=bad)
        with self.assertRaises(ValueError):
            dual.encode_global((), ())
        with self.assertRaises(ValueError):
            dual.encode_local(region, (dual.ZERO_MATRIX,))


if __name__ == "__main__":
    unittest.main()
