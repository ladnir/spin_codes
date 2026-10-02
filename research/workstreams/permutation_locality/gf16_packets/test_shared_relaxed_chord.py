from fractions import Fraction as Q
from math import comb
import unittest

from flint import arb,arb_mat,ctx
import shared_relaxed_chord as c


class ChordTests(unittest.TestCase):
    def setUp(self):
        self.precision=ctx.prec;ctx.prec=256
    def tearDown(self):ctx.prec=self.precision

    def test_concave_maximum_against_exact_stationary_point(self):
        value,witness=c.concave_maximum(3,10,0,(Q(1,10),Q(9,10)))
        exact=3*c.aq(Q(3,10)).log()+7*c.aq(Q(7,10)).log()
        self.assertTrue(c.aq(value)>=exact)
        self.assertLess(float(c.aq(value)-exact),1e-20)
        for penalty in (-100,100):
            value,_=c.concave_maximum(3,10,penalty,(Q(1,10),Q(9,10)))
            for k in range(1,10):
                m=Q(k,10)
                exact=3*c.aq(m).log()+7*c.aq(1-m).log()-c.aq(penalty*m)
                self.assertTrue(c.aq(value)>=exact)

    def test_scalar_polynomial_chord_bounds_whole_interval(self):
        local=[arb_mat([[c.aq(Q(1))]]),arb_mat([[c.aq(Q(1,3))]]),arb_mat([[c.aq(Q(2,7))]])]
        cell=(Q(1,20),Q(4,5));x=Q(3,16);epochs=4
        hlo,hhi=c.transition_endpoints(local,epochs,cell,x)
        bound,_=c.chord_from_endpoints(hlo,hhi,8,cell,x,3)
        for j in range(31):
            m=cell[0]+(cell[1]-cell[0])*Q(j,30)
            weighted=(1-m)**2+Q(2,3)*(m/x)*(1-m)+Q(2,7)*(m/x)**2
            exact=epochs*c.aq(weighted).log()-c.aq(3*m)
            self.assertTrue(c.aq(bound)>=exact)

    def test_noncommuting_positive_matrices_use_one_polynomial(self):
        local=[arb_mat([[1,1],[0,1]]),arb_mat([[0,2],[1,0]]),arb_mat([[1,0],[1,1]])]
        cell=(Q(1,10),Q(3,5));x=Q(1,4);epochs=3
        hlo,hhi=c.transition_endpoints(local,epochs,cell,x)
        bound,_=c.chord_from_endpoints(hlo,hhi,6,cell,x)
        for k in range(11):
            m=cell[0]+(cell[1]-cell[0])*Q(k,10)
            matrix=sum((mat*c.aq(comb(2,j)*(m/x)**j*(1-m)**(2-j))
                        for j,mat in enumerate(local)),arb_mat(2,2))
            power=matrix**epochs
            exact=(power[0,0]+power[0,1]).log()
            self.assertTrue(c.aq(bound)>=exact)

    def test_clipped_slope_preserves_degree_growth_bound(self):
        # Deliberately inflate one endpoint bound. Clipping above P is
        # justified by H(v)=v^P, not by the noisy chord slope.
        cell=(Q(1,4),Q(3,4));x=Q(1,2);degree=4
        vlo=cell[0]/(x*(1-cell[0]));vhi=cell[1]/(x*(1-cell[1]))
        bound,info=c.chord_from_endpoints(vlo**degree,100*vhi**degree,degree,cell,x)
        self.assertEqual(Q(info['slope']),degree)
        self.assertTrue(c.aq(bound)>=degree*c.aq(cell[1]/x).log())
        # An inflated low endpoint can make the chord slope negative.
        bound,info=c.chord_from_endpoints(100,1,degree,cell,x)
        self.assertEqual(Q(info['slope']),0)
        self.assertTrue(c.aq(bound)>=degree*c.aq(1-cell[0]).log())

    def test_singleton_and_invalid_geometry(self):
        bound,_=c.chord_from_endpoints(4,4,8,(Q(1,3),Q(1,3)),Q(1,4),2)
        exact=arb(4).log()+8*c.aq(Q(2,3)).log()-c.aq(Q(2,3))
        self.assertTrue(c.aq(bound)>=exact)
        for cell in ((Q(0),Q(1,2)),(Q(1,2),Q(1)),(Q(3,4),Q(1,4))):
            with self.assertRaises(ValueError):c.chord_from_endpoints(1,2,8,cell,Q(1,4))
        with self.assertRaises(ValueError):c.transition_endpoints([arb_mat([[1]]),arb_mat([[-1]])],4,(Q(1,4),Q(3,4)),Q(1,4))

    def test_large_endpoint_metadata_remains_compact(self):
        self.assertEqual(c.compact_endpoint(Q(3)*Q(2)**100000),dict(dyadic=[3,100000]))
        self.assertEqual(c.compact_endpoint(Q(3)*Q(2)**-100000),dict(dyadic=[3,-100000]))


if __name__=='__main__':unittest.main()
