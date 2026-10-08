"""Small exact controls for the variable-rate sparse replay."""
from fractions import Fraction as F
import unittest

from flint import arb, ctx

from field128_sparse_certificate import certify_zone as frozen_zone
from gf8_outer_rate_sparse import (certify_cover, certify_point, certify_zone, metadata,
    rank_polynomial_power, rank_terminal_power, rate_geometry)
from joint_outer_screen import WITNESSES_128_87, packet
from state8_sparse_certificate import ball, polynomial_power, terminal_power


class SparseRateTests(unittest.TestCase):
    def test_geometry(self):
        expected = {(16,88):(768,98304),(16,94):(704,90112),
                    (18,88):(3008,385024),(20,88):(11968,1531904)}
        for key, value in expected.items():
            g = rate_geometry(*key)
            self.assertEqual((g['rows'], g['n']), value)
            self.assertGreaterEqual(g['k'], 2**key[0])
            self.assertEqual(g['threshold'], (g['n']+9)//10-1)
        for dimension in (0,128,129):
            with self.assertRaises(ValueError):
                rate_geometry(20,dimension)

    def test_frozen_dimension87(self):
        mark = WITNESSES_128_87[1]
        old, _ = frozen_zone(1,1,1,mark,rate_geometry(20,87),precision=1024)
        new = certify_zone(1,1,1,mark,20,87,precision=1024)
        p,q = map(int,new['failure_upper_dyadic'])
        old_exact=old.fmpq()
        self.assertEqual(F(p,q), F(int(old_exact.p),int(old_exact.q)))

    def test_generic_polynomial_matches_two_state(self):
        previous=ctx.prec
        try:
            ctx.prec=512
            local=[packet(w,F(3,4),F(8),joint=True) for w in range(17)]
            for packets,degree in ((1,3),(2,5),(3,7),(5,4)):
                old=polynomial_power(local,packets,degree)
                new=rank_polynomial_power(local,packets,degree)
                for a,b in zip(old,new):
                    for h in range(degree+1):
                        self.assertTrue(a[h].overlaps(b[h]))
        finally:
            ctx.prec=previous

    def test_generic_terminal_matches_two_state(self):
        values=tuple(map(arb,(2,3,5,7)))
        for exponent in (0,1,2,3,8,17):
            self.assertEqual(rank_terminal_power(values,exponent),terminal_power(values,exponent))

    def test_rank_cone_exact(self):
        from gf8_rank_history_screen import local_matrix, cone_preserved
        for weight in (0,1,4,8,16):
            self.assertTrue(cone_preserved(local_matrix(weight,F(3,4),F(8),state=8)))

    def test_metadata_scope(self):
        for size in (16,18,20):
            receipt=metadata(size,88)
            self.assertEqual(receipt['outer'],[128,88])
            self.assertTrue(receipt['rigorous'])
            self.assertFalse(receipt['whole_code'])
            self.assertEqual(receipt['requested_k'],2**size)

    def test_complete_k18_dimension88(self):
        result=certify_cover(375,18,88,precision=4096,cell_bits=64)
        self.assertEqual(result['active_range'],[1,375])
        self.assertTrue(result['complete'])
        self.assertTrue(result['rigorous'])
        self.assertFalse(result['whole_code'])
        self.assertGreaterEqual(result['margin_bits_lower'],72)
        self.assertEqual(result['interval_count'],67)
        self.assertEqual(result['partition_sha256'],
                         '0f633ee9ea876ec634d2fcf0916117c6610df7bd75c9f9ddb504a10f5f37ea38')
        p,q=map(int,result['failure_upper_dyadic'])
        self.assertLess(F(p,q),F(1,2**72))


if __name__=='__main__':
    unittest.main()
