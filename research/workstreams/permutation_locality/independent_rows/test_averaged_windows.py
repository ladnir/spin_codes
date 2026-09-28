"""Exact enumeration checks for categorical distinct-window moments."""
from fractions import Fraction as Q
from itertools import combinations, product
from math import comb, prod
import unittest

from flint import arb, ctx
from averaged_windows import histogram_moments


class WindowMoments(unittest.TestCase):
    def test_direct_packets(self):
        for pattern in ((0,4),(1,3),(0,2,4),(1,2,3)):
            hist = tuple(pattern.count(r) for r in range(5))
            image = sum(((1<<r)-1)<<(4*i) for i,r in enumerate(pattern))
            z = Q(7,8)
            for theta in ((Q(1),Q(0),Q(0),Q(0)),(Q(0),Q(0),Q(0),Q(1)),(Q(1,10),Q(2,10),Q(3,10),Q(4,10))):
                values = histogram_moments(hist,theta,z)
                for j in range(len(pattern)+1):
                    direct = Q(0)
                    for slots in combinations(range(len(pattern)),j):
                        for masks in product(range(1,16),repeat=j):
                            probability = prod(theta[mask.bit_count()-1]/comb(4,mask.bit_count()) for mask in masks)
                            word = sum(mask<<(4*i) for i,mask in zip(slots,masks))
                            direct += probability*z**((image^word).bit_count())
                    direct /= comb(len(pattern),j)
                    self.assertEqual(values[j],direct)

    def test_outward_coefficients(self):
        ctx.prec = 128
        hist,theta,z = (1,1,1,1,1),(Q(1,10),Q(2,10),Q(3,10),Q(4,10)),Q(11,13)
        exact = histogram_moments(hist,theta,z)
        outward = histogram_moments(hist,[arb(p.numerator)/p.denominator for p in theta],
                                    arb(z.numerator)/z.denominator,lambda x: arb(x.upper()))
        for a,b in zip(exact,outward):
            endpoint = b.fmpq()
            self.assertGreaterEqual(Q(int(endpoint.p),int(endpoint.q)),a)


if __name__ == '__main__':
    unittest.main()
