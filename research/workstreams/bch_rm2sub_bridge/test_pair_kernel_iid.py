import itertools
import unittest
from fractions import Fraction as F
from inner_pair_spectrum import pair_counts
from pair_kernel_iid import probability,single


class KernelTests(unittest.TestCase):
    def test_fourier_against_direct_kernel(self):
        generators=[3,5];length=4;pairs,spectrum=pair_counts(generators)
        kernel=[x for x in range(1<<length) if all((x&g).bit_count()%2==0 for g in generators)]
        for masses in [(F(1,4),)*4,(F(1,2),F(1,6),F(1,6),F(1,6)),(F(0),F(1,2),F(1,2),F(0))]:
            direct=sum((__import__('math').prod(masses[2*((x>>j)&1)+((y>>j)&1)] for j in range(length))
                for x in kernel for y in kernel),F(0))
            self.assertEqual(probability(masses,length,2,pairs),direct)
        for p,q in itertools.product((F(1,8),F(1,3),F(3,4)),repeat=2):
            masses=((1-p)*(1-q),(1-p)*q,p*(1-q),p*q)
            self.assertEqual(probability(masses,length,2,pairs),single(p,spectrum,2)*single(q,spectrum,2))


if __name__=='__main__':unittest.main()
