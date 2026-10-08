"""Independent rank/zero-mask checks of the global-MDS moment formulas."""
from fractions import Fraction as F
import random
import unittest
from hd_blockwise_shared_transfer import mm,rank
from hd_sumrank_maps import Q
from hd_systematic_mds import maps,moments,matrix


class SystematicMdsTests(unittest.TestCase):
    def test_packet_sort_preserves_count_distribution(self):
        from itertools import permutations,product
        from collections import Counter
        from math import comb,factorial,prod
        for r in range(1,6):
            observed=Counter()
            for perm in permutations(range(6)):
                packets=[tuple(sorted(perm[i:i+2])) for i in range(0,6,2)]
                counts=tuple(sum(x<r for x in packet) for packet in packets)
                observed[counts]+=1
            for counts in product(range(3),repeat=3):
                expected=factorial(6)*prod(comb(2,j) for j in counts)//comb(6,r) if sum(counts)==r else 0
                self.assertEqual(observed[counts],expected)

    def test_self_orthogonality(self):
        for width,state in ((8,4),(64,8),(64,16),(64,32)):
            a,c=maps(state,width)
            self.assertEqual(mm(c,a,Q),((0,)*state,)*state)
            self.assertEqual(rank(c,Q),state)

    def test_actual_global_mask_ranks(self):
        a,c=maps();rng=random.Random(146271)
        for _ in range(100):
            u=sorted(rng.sample(range(64),rng.randrange(65)))
            z=sorted(rng.sample(range(64),rng.randrange(65)))
            overlap=len(set(u)&set(z))
            f=tuple(tuple(int(i==j) for j in u) for i in z)
            quiet=f+tuple(tuple(row[j] for j in u) for row in c)
            flat=tuple(row+a[i] for row,i in zip(f,z))
            self.assertEqual(rank(quiet,Q),overlap+min(16,len(u)-overlap))
            self.assertEqual(rank(flat,Q),overlap+min(16,len(z)-overlap))

    def test_all_small_masks(self):
        width,state,z=8,4,F(3,4);a,c=maps(state,width)
        for j in range(width+1):
            u=tuple(range(j));expected=[F(0)]*3
            for mask in range(1<<width):
                out=tuple(i for i in range(width) if mask>>i&1)
                f=tuple(tuple(int(i==x) for x in u) for i in out)
                quiet=rank(f+tuple(tuple(row[x] for x in u) for row in c),Q)
                joint=rank(tuple(row+a[i] for row,i in zip(f,out)),Q)
                base=len(set(u)&set(out))
                weight=z**(width-len(out))*(1-z)**len(out)
                expected[0]+=weight*F(1,Q**quiet)
                expected[1]+=weight*F(1,Q**base)
                expected[2]+=weight*F(Q**(state-joint+base)-1,Q**base*(Q**state-1))
            self.assertEqual(tuple(expected),moments(j,z,state,width))

    def test_float_screen_and_grouping(self):
        import numpy as np
        from math import comb
        from hd_systematic_mds_screen import packet_logs
        from hd_sumrank_paired_certificate import grouped_packets
        logs=packet_logs(.75,16)
        groups=grouped_packets(F(3,4),feedthrough='systematic_mds')
        self.assertEqual(len(groups),65)
        for j,count,rounded in groups:
            self.assertEqual(count,comb(64,j))
            exact=matrix(j,F(3,4))
            self.assertTrue(all(x<=y for x,y in zip(exact,rounded)))
            for x,y in zip(exact,np.exp2(logs[j]).reshape(-1)):
                if x:self.assertAlmostEqual(float(x)/y,1,places=11)
                else:self.assertEqual(y,0)

    def test_smaller_packet_screen(self):
        import numpy as np
        from hd_systematic_mds_screen import packet_logs
        for width,state in ((8,4),(32,8)):
            logs=packet_logs(.75,state,width=width)
            self.assertEqual(logs.shape,(width+1,2,2))
            for j in range(width+1):
                for x,y in zip(matrix(j,F(3,4),state,width),np.exp2(logs[j]).reshape(-1)):
                    if x:self.assertAlmostEqual(float(x)/y,1,places=11)
                    else:self.assertEqual(y,0)

    def test_wider_packet_screen(self):
        import numpy as np
        from hd_systematic_mds_screen import packet_logs
        from hd_sumrank_paired_certificate import grouped_packets
        for width,state in ((128,20),(128,24),(256,32)):
            logs=packet_logs(.75,state,width=width)
            groups=grouped_packets(F(3,4),state=state,feedthrough='systematic_mds',width=width)
            self.assertEqual(len(groups),width+1)
            for j in (0,1,state,width//2,width):
                exact=matrix(j,F(3,4),state,width)
                self.assertTrue(all(x<=y for x,y in zip(exact,groups[j][2])))
                for x,y in zip(exact,np.exp2(logs[j]).reshape(-1)):
                    if x:self.assertAlmostEqual(float(x)/y,1,places=10)
                    else:self.assertEqual(y,0)

    def test_general_pair_coefficient_screen(self):
        import numpy as np
        from math import comb
        from hd_systematic_mds_screen import coefficients
        rng=np.random.default_rng(425)
        for size in (2,3):
            a=rng.uniform(.1,.9,(9,size,size));b=rng.uniform(.1,.9,(9,size,size))
            actual=np.exp2(coefficients(np.log2(a),.5,np.log2(b)))
            expected=np.zeros((17,size,size))
            for j in range(9):
                for k in range(9):expected[j+k]+=comb(8,j)*comb(8,k)*np.sqrt(a[j]@b[k])
            np.testing.assert_allclose(actual,expected,rtol=1e-12)

    def test_general_outer_counting_normalization(self):
        import numpy as np
        from math import comb,log2
        from hd_sumrank_profile_pairs import exact_region_screen
        width=128;transition=np.array(((.6,.4),(.2,.8)))
        coeff=np.array([log2(comb(width,j))+np.log2(transition) for j in range(width+1)])
        for outer_length in (64,128,256):
            for r in (3,13):
                actual=exact_region_screen(r,.75,.5,coeff,width=width,outer_length=outer_length)
                length=2**21//outer_length
                expected=log2(comb(length,r))+.5*(r*(16*outer_length-(outer_length//2)*log2(65535))-log2(65535)-838860*log2(.75))
                self.assertAlmostEqual(actual['log2_upper_screen'],expected,places=6)


if __name__=='__main__':unittest.main()
