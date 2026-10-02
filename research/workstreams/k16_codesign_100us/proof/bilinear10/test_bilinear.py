"""Fresh independent algebra and physical/macro checks for the s10 candidate."""
from collections import Counter
from fractions import Fraction as Q
import unittest

import bilinear_maps as maps
from flint import arb,ctx


def polynomial_product(a,b):
    product=0
    for i in range(3):
        for j in range(3):
            if (a>>i&1) and (b>>j&1):product^=1<<(i+j)
    for degree in range(4,2,-1):
        if product>>degree&1:product^=0xb<<(degree-3)
    return product


class BilinearTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        ctx.prec=256
        cls.wrapper,cls.record=maps.prepare()
        cls.rows=[int(x,16) for x in cls.record['expansion_rows_hex']]
        cls.columns=cls.record['feedback_columns']

    def test_field_and_product_truth_table(self):
        for x in range(8):
            for y in range(8):
                value=polynomial_product(x,y)
                self.assertEqual(value,maps.multiply8(x,y))
                self.assertEqual(value==0,x==0 or y==0)
        for z in range(64):
            x=(z&1)|((z>>1)&2)|((z>>2)&4)
            y=((z>>1)&1)|((z>>2)&2)|((z>>3)&4)
            expected=[1,*[(z>>i)&1 for i in range(6)],
                      *[(polynomial_product(x,y)>>i)&1 for i in range(3)]]
            self.assertEqual(sum(v<<i for i,v in enumerate(expected)),self.columns[z])

    def test_every_state_weight_and_polar_form(self):
        weights=Counter(sum((state&column).bit_count()%2 for column in self.columns)
                        for state in range(1024))
        self.assertEqual(weights,Counter({0:1,28:448,32:126,36:448,64:1}))
        # Nondegeneracy directly: for every nonzero direction, the polar
        # derivative of each nonzero quadratic combination is nonconstant.
        for selector in range(1,8):
            q=[(maps.product_at(z)&selector).bit_count()%2 for z in range(64)]
            for direction in range(1,64):
                self.assertEqual({q[z]^q[z^direction]^q[direction]^q[0] for z in range(64)},{0,1})

    def test_CA_zero_and_every_packet_injective(self):
        self.assertTrue(all((a&b).bit_count()%2==0 for a in self.rows for b in self.rows))
        for start in range(0,64,4):
            images=set()
            for label in range(16):
                value=0
                for bit in range(4):
                    if label>>bit&1:value^=self.columns[start+bit]
                images.add(value)
            self.assertEqual(len(images),16)

    def test_two_physical_steps_not_one_t128(self):
        kernel=maps.screen.q1.kernel_t64
        physical=kernel.authenticate(self.wrapper)
        z=(-kernel.aq(Q(1,200))).exp()
        local=kernel.sparse_kernel.outward_at_z(physical,z)
        local=kernel.kernel_birth_density.refine_local(physical,local,z,Q(1,2))
        macro=kernel.local_operators(self.wrapper,Q(1,200))
        q1=maps.screen.q1
        direct=q1.placement(local,epochs=32,windows=16,maximum_groups=2,rounding=q1.rounded)
        composed=q1.placement(macro,epochs=16,windows=32,maximum_groups=2,rounding=q1.rounded)
        for q in range(3):
            for i in range(direct[q].nrows()):
                for j in range(direct[q].ncols()):
                    self.assertLess(abs(direct[q][i,j]-composed[q][i,j]),arb(2)**-220)


if __name__=='__main__':
    unittest.main()
