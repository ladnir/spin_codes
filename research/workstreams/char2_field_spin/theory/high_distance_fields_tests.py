from fractions import Fraction as F
import unittest
from gf8_punctured_screen import local_matrices as frozen
from high_distance_restart_screen import adjusted_bound
from high_distance_fields import alphabet,bound,local_matrices,roots,rectangular_cauchy_maps


class FieldScreenTests(unittest.TestCase):
    def test_exact_gf256_local_agreement(self):
        for y in (F(1,3),F(9,10),F(1)):
            for tau in (F(1),F(17),F(255)):
                for outputs in (1,2,3,4):
                    self.assertEqual(local_matrices(y,tau,4,outputs,4,8),frozen(y,tau,4,outputs,4))

    def test_gf256_complete_screen_agreement(self):
        for active in (1,8,16):
            self.assertAlmostEqual(bound((5.,2.),active),
                adjusted_bound((5.,2.),active,distance=F(9,25),outputs=2,dimension=16),places=7)
        point=(-2.,-.3,2.,0.)
        self.assertAlmostEqual(bound(point,100000,False),
            adjusted_bound(point,100000,False,F(9,25),outputs=2,dimension=16),places=7)

    def test_minimum_weight_marker_zero_agrees(self):
        self.assertEqual(bound((5.,2.,0.),8),bound((5.,2.),8))
        with self.assertRaises(ValueError):bound((5.,2.,-1.),8)

    def test_larger_alphabet_improves_roots_at_matched_marks(self):
        previous=None
        for bits in (8,16,32,64,128):
            current=roots(F(17),4,4,alphabet(bits))
            self.assertEqual(current[0],1)
            self.assertTrue(all(0<=v<=1 for v in current))
            if previous is not None:
                self.assertTrue(all(a<=b for a,b in zip(current,previous)))
            previous=current

    def test_large_field_moments_do_not_overflow(self):
        import math
        values=roots(17.0,4,16,1<<128)
        self.assertTrue(all(math.isfinite(v) and 0<=v<=1 for v in values))
        self.assertEqual(values,tuple(map(float,roots(F(17),4,16,1<<128))))

    def test_dense_endpoint_is_not_middle(self):
        from high_distance_fields import search
        with self.assertRaises(ValueError):search(1,False,bits=32)
        with self.assertRaises(ValueError):bound((-2.,-.3,2.,0.),4194304,False)

    def test_full_field_law_and_zero_state_normalization(self):
        with self.assertRaises(ValueError):local_matrices(F(1,2),F(256),bits=8)
        self.assertTrue(local_matrices(F(1,2),F(256),bits=16))
        with self.assertRaises(ValueError):alphabet(12)
        for bits in (8,16,32):
            matrices=local_matrices(F(1),F(1),bits=bits)
            self.assertEqual(matrices[0][0],[1,0])
            self.assertTrue(all(v>=0 for matrix in matrices for row in matrix for v in row))

    def test_rectangular_state8_cauchy_invariants(self):
        from gf8_packet_maps import matmul,rank
        from gf8_mds_packet import all_minors,identity,add
        from random import Random
        f,a,c,raw,e,ei=rectangular_cauchy_maps()
        transformed=tuple(f[i]+a[i] for i in range(4))+tuple(c[i]+identity(8)[i] for i in range(8))
        normalizer=identity(12)[:4]+tuple((0,)*4+row for row in ei)
        self.assertEqual(matmul(normalizer,raw),transformed)
        for matrix in (f,a,c,add(f,matmul(a,c))):
            self.assertGreater(all_minors(matrix),0)
        rng=Random(48012)
        for _ in range(256):
            rows=tuple(i for i in range(4) if rng.randrange(2))
            if rng.randrange(2):rows+=tuple(range(4,12))
            cols=tuple(i for i in range(12) if rng.randrange(2))
            self.assertEqual(rank(tuple(tuple(transformed[i][j] for j in cols) for i in rows)),min(len(rows),len(cols)))
        for bits in (8,32):
            terms=local_matrices(F(4,5),F(17),4,2,8,bits)
            self.assertTrue(all(v>=0 for matrix in terms for row in matrix for v in row))


if __name__=='__main__':
    unittest.main()
