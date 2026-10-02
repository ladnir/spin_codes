import unittest
from fractions import Fraction as Q

from flint import arb,ctx

import post_mix as candidate
from test_fresh_history import endpoint


class PostMixTests(unittest.TestCase):
    def test_zero_and_auxiliary_coordinates(self):
        spectrum={48:3,56:7,64:11,72:13,80:17}
        for rounds in (0,1,2,3,None):
            matrix=candidate.refresh(spectrum,rounds)
            a=Q(0) if rounds is None else Q(1,2**rounds)
            self.assertEqual(matrix[0,0],1)
            for i in range(1,11):self.assertEqual(matrix[i,0],0)
            for source in (3,9,10):
                for target in range(11):
                    self.assertEqual(endpoint(matrix[source,target]),a if source==target else 0)

    def test_uniform_refresh_matches_direct_nonzero_kernel(self):
        ctx.prec=192
        spectrum={48:1,56:2,64:3,72:4,80:5};total=sum(spectrum.values())
        labels=[v for v,n in spectrum.items() for _ in range(n)]
        for rounds in (0,1,3,None):
            matrix=candidate.refresh(spectrum,rounds)
            a=Q(0) if rounds is None else Q(1,2**rounds)
            for source,level in enumerate(spectrum,4):
                states=[s for s,v in enumerate(labels) if v==level]
                for target,v in enumerate(spectrum,4):
                    mass=sum((a*int(s==t)+(1-a)/total)
                             for s in states for t,w in enumerate(labels) if w==v)/len(states)
                    upper=endpoint(matrix[source,target])
                    self.assertGreaterEqual(upper,mass)
                    self.assertLess(upper-mass,Q(1,2**180))

    def test_mass_is_counted_once_not_per_auxiliary_bound(self):
        ctx.prec=192
        spectrum={48:1,56:2,64:3,72:4,80:5}
        matrix=candidate.refresh(spectrum,2)
        terminal=[1,1,1,0,1,1,1,1,1,0,0]
        for source in range(11):
            mass=sum((matrix[source,j]*terminal[j] for j in range(11)),arb(0))
            self.assertTrue(mass.contains(terminal[source]))


if __name__=='__main__':unittest.main()
