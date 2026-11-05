from itertools import product
import unittest
import numpy as np
from hd_pair_singleton_word import whole_word_roots


class SingletonWordTests(unittest.TestCase):
    def test_scalar_direct_order_counts_and_floor(self):
        blocks = [np.arange(1.,33.).reshape(2,4,2,2),
                  np.arange(33.,65.).reshape(2,4,2,2)]
        expected = []
        for source in range(2):
            total = 0.
            for bits in product((0,1),repeat=4):
                if sum(bits) < 3:
                    continue
                matrix = blocks[0][source,2*bits[0]+bits[1]] @ blocks[1][source,2*bits[2]+bits[3]]
                total += 3.**sum(bits)*np.sum(matrix[0,:])**.5
            expected.append(total)
        got = whole_word_roots([np.log2(x) for x in blocks],2,minimum_activity=3,q=4)
        np.testing.assert_allclose(2.**got,expected,rtol=2e-14)

    def test_zero_word_removed_and_tiny_term_preserved(self):
        block = np.full((1,4,2,2),-np.inf)
        block[0,0,0,0] = 0.
        block[0,3,0,0] = -18000.
        got = whole_word_roots([block],18,minimum_activity=2,q=2)
        self.assertEqual(got[0],-1000.)

    def test_r1_matches_direct_valid_polynomial_sum(self):
        rng = np.random.default_rng(314)
        blocks = [rng.uniform(.1,2.,size=(1,4,2,2)) for _ in range(3)]
        expected = 0.
        for mask in range(64):
            if mask.bit_count()<3:
                continue
            row = np.array([1.,0.])
            for b,block in enumerate(blocks):
                index = 2*((mask>>(2*b))&1)+((mask>>(2*b+1))&1)
                row = row@block[0,index]
            expected += 7.**mask.bit_count()*sum(row)
        got = whole_word_roots([np.log2(x) for x in blocks],1,minimum_activity=3,q=8)
        self.assertAlmostEqual(got[0],np.log2(expected),places=12)


if __name__ == '__main__':
    unittest.main()
