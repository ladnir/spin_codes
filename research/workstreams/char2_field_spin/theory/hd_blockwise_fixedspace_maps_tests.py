import itertools
import random
import unittest

from hd_blockwise_maps import matmul, matvec
from hd_blockwise_shared_maps import rank, identity
from hd_blockwise_fixedspace_maps import (
    fixedspace_maps, local_output_rank, singleton_cancellation,
    global_quiet_not_fixed, z_value, weight,
)


class FixedSpaceMapsTests(unittest.TestCase):
    def test_factorization_and_mds_identity(self):
        f, a, c, k = fixedspace_maps()
        self.assertEqual(k, ((17,114,0,240),(114,0,240,0),
                             (0,240,0,0),(240,0,0,0)))
        self.assertEqual(rank(k), 4)
        self.assertEqual(len({z_value(x) for x in range(64)}), 64)
        for offset in range(0,64,8):
            at = a[offset:offset+8]
            ct = tuple(row[offset:offset+8] for row in c)
            residual = tuple(tuple(f[offset+i][offset+j] ^ int(i == j)
                                   for j in range(8)) for i in range(8))
            self.assertEqual(matmul(at,ct),residual)
            self.assertEqual(matmul(ct,at),((0,)*4,)*4)
        self.assertEqual(matmul(f,a),a)
        self.assertEqual(matmul(c,f),c)
        self.assertEqual(matmul(c,a),((0,)*4,)*4)
        residual = tuple(tuple(v ^ int(i == j) for j,v in enumerate(row))
                         for i,row in enumerate(f))
        self.assertEqual(rank(c+residual),32)
        # Vandermonde nodes prove all global four-row minors; bounded exact
        # checks include every choice inside each tile and across two tiles.
        for offset in range(0,64,8):
            for rows in itertools.combinations(range(offset,offset+8),4):
                self.assertEqual(rank(tuple(a[i] for i in rows)),4)
        for rows in itertools.combinations(range(16),4):
            self.assertEqual(rank(tuple(a[i] for i in rows)),4)

    def test_exact_local_rank_formula(self):
        rng = random.Random(41)
        cases = [((j,),tuple(i for i in range(8) if i != j)) for j in range(8)]
        cases += [(tuple(i for i in range(8) if rng.randrange(2)),
                   tuple(i for i in range(8) if rng.randrange(2))) for _ in range(512)]
        for tile in (0,1,7):
            for support,zeros in cases:
                result=local_output_rank(support,zeros,tile)
                self.assertEqual(result['exact_rank'],result['formula'])

    def test_local_singleton_does_not_imply_packet_singleton(self):
        for j in range(64):
            result=singleton_cancellation(j)
            self.assertEqual(result['local_joint_weight'],2)
            self.assertGreaterEqual(result['packet_joint_weight'],37)

    def test_global_feedback_zero_is_not_local_feedback_zero(self):
        result=global_quiet_not_fixed()
        self.assertEqual(result['global_feedback'],(0,)*4)
        self.assertEqual(result['input_weight'],5)
        self.assertNotEqual(result['input'],result['output'])
        self.assertTrue(any(any(v) for v in result['local_feedback']))

    def test_packet_joint_weight_bound_sample(self):
        f,a,c,_=fixedspace_maps()
        rng=random.Random(73)
        for _ in range(64):
            q=tuple(rng.randrange(256) for _ in range(4))
            x=tuple(rng.randrange(256) if rng.randrange(3)==0 else 0 for _ in range(64))
            y=tuple(u^v for u,v in zip(matvec(f,x),matvec(a,q)))
            self.assertGreaterEqual(weight(x)+weight(y),37)
            self.assertEqual(matvec(c,y),matvec(c,x))


if __name__=='__main__':
    unittest.main()
