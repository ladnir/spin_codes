"""Exact bounded tests: no native build or timing measurements."""

import random
import unittest

from hd_blockwise_maps import matmul, matvec, mul, power, xor_sum
from hd_blockwise_shared_maps import evaluation, identity, inverse, minor_count, rank
from hd_blockwise_shared_factor import (
    MOMENT_MASKS, MOMENT_TO_HIGH, accumulate_moments, cube_coefficients,
    cube_differences, cube_feedthrough, emission_coefficients, evaluation_at,
    matrix_xor, moment_coefficients, packet_circuit, pair_factor,
    polynomial_factor, product_receipt, reference_packet, tile_matrix,
)


class SharedFactorTests(unittest.TestCase):
    def test_both_rank_four_factorizations_and_feedback_defect(self):
        for beta in (8, 64, 128):
            n = matrix_xor(tile_matrix(beta), identity(8))
            for u, v in (pair_factor(beta)[:2], polynomial_factor(beta)):
                self.assertEqual(matmul(u,v), n)
                self.assertEqual(matmul(v,u), ((0,)*4,)*4)
                self.assertEqual(rank(u),4)
                self.assertEqual(rank(v),4)
                self.assertEqual(rank(v+inverse(evaluation(0))[4:]),6)
            self.assertEqual(matmul(tile_matrix(beta),tile_matrix(beta)),identity(8))

    def test_beta64_branch_nine_all_minors(self):
        self.assertEqual(minor_count(tile_matrix(64)),12869)

    def test_cube_feedthrough_all_basis_vectors(self):
        for beta in (8,64,128):
            for j in range(8):
                values = tuple(int(i==j) for i in range(8))
                self.assertEqual(cube_feedthrough(values,beta),
                                 matvec(tile_matrix(beta),values))
            self.assertEqual(len(cube_differences(tuple(range(8)))),19)
            self.assertEqual(cube_coefficients(beta)[0],1)

    def test_shared_global_moments_and_exact_feedback_identity(self):
        values = tuple((x*x+19*x+7)%256 for x in range(64))
        differences = tuple(cube_differences(values[i:i+8]) for i in range(0,64,8))
        moments = accumulate_moments(differences)
        for mask in MOMENT_MASKS:
            self.assertEqual(moments[mask],xor_sum(values[x] for x in range(64)
                                                 if x&mask==mask))
        for alpha in range(0,64,8):
            low = tuple(tuple(row[:4]) for row in evaluation(alpha))
            self.assertEqual(matmul(MOMENT_TO_HIGH,tuple(zip(*low))),
                             inverse(evaluation(alpha))[4:])

    def test_complete_packet_all_input_and_state_basis_vectors(self):
        for beta in (8,64,128):
            for j in range(68):
                values = tuple(int(i==j) for i in range(64))
                state = tuple(int(i+64==j) for i in range(4))
                self.assertEqual(packet_circuit(values,state,beta),
                                 reference_packet(values,state,beta))

    def test_nonzero_state_random_and_gf65536_extension(self):
        rng = random.Random(0xF864)
        for beta in (8,64,128):
            for _ in range(6):
                values = tuple(rng.randrange(65536) for _ in range(64))
                state = tuple(rng.randrange(65536) for _ in range(4))
                # Embedded GF256 multiplication acts independently on the two
                # coefficients in the existing GF65536 quadratic tower.
                for shift in (0,8):
                    v = tuple((x>>shift)&255 for x in values)
                    q = tuple((x>>shift)&255 for x in state)
                    self.assertEqual(packet_circuit(v,q,beta),reference_packet(v,q,beta))

    def test_product_accounting_and_degree_two_support(self):
        for beta,total,emission in ((8,58,27),(64,59,32),(128,59,32)):
            result = product_receipt(beta)
            self.assertEqual(result['feedthrough_record_products_per_tile'],19)
            self.assertEqual(result['feedback_record_products_per_packet'],28)
            self.assertEqual(result['emission_record_products_per_packet'],emission)
            self.assertEqual(result['packed_arithmetic_GFNI_per_packet'],total)
            for row in moment_coefficients()+emission_coefficients(beta):
                self.assertFalse(any(c for mask,c in enumerate(row) if mask.bit_count()>2))


if __name__ == '__main__':
    unittest.main()
