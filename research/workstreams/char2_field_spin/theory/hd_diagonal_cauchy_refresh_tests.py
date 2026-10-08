import unittest
from collections import Counter
from itertools import product

from hd_blockwise_shared_transfer import multiply,power,dot,rank
from hd_diagonal_cauchy_refresh import transfer


class DiagonalCauchyTests(unittest.TestCase):
    def test_every_small_field_vector_and_diagonal(self):
        # Over GF4, this 2x2 H has no zero entry and nonzero determinant.
        q=4
        h=((1,1),(1,2))
        self.assertEqual(rank(h,q),2)
        diagonals=tuple(product(range(1,q),repeat=2))
        for x in product(range(q),repeat=2):
            if x==(0,0):continue
            counts=Counter()
            for d1 in diagonals:
                a=tuple(multiply(t,v,q) for t,v in zip(d1,x))
                a=tuple(dot(row,a,q) for row in h)
                for d2 in diagonals:
                    counts[tuple(multiply(t,v,q) for t,v in zip(d2,a))]+=1
            self.assertEqual(counts[(0,0)],0)
            self.assertLessEqual(max(counts.values()),(q-1)**2)
            self.assertEqual(sum(counts.values()),(q-1)**4)

    def test_directed_transfer_closes(self):
        from flint import arb,ctx
        with ctx.workprec(256):
            r=transfer()
            self.assertTrue(r['whole_code_certificate'])
            self.assertTrue(arb(r['loss_bits'])<12)
            self.assertTrue(arb(r['union_log2_upper']) < -48)
            self.assertFalse(r['underlying_certificate_replayed'])


if __name__=='__main__':unittest.main()
