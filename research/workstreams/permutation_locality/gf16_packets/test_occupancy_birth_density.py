from fractions import Fraction as Q
from math import comb
import unittest
from flint import arb,ctx
import occupancy_birth_density_probe as density
from test_scalar_cover import endpoint


class FixedBirthDensityTests(unittest.TestCase):
    def test_feedback_caps_and_ordered_products(self):
        ctx.prec=192;z=Q(3,4)
        images=[sum(((a>>i)&1)*v for i,v in enumerate((1,6,120))) for a in range(8)]
        for columns in ([1,2,4,3,5,7,6,1],[0]*8,[1]*8):
            for updates in (2,3):
                data=density.base.prepare(images,columns,3,updates)
                local=density.base.outward_at_z(data,arb(3)/4)
                caps=density.feedback_caps(data,arb(3)/4)
                feedbacks=[[Q(0)]*8 for _ in range(3)]
                exact=[[[Q(0)]*8 for _ in range(8)] for _ in range(3)]
                for x in range(256):
                    j=int(bool(x&15))+int(bool(x>>4));chance=Q(1,comb(2,j)*15**j)
                    feedback=0
                    for bit,c in enumerate(columns):
                        if x>>bit&1:feedback^=c
                    feedbacks[j][feedback]+=chance*z**x.bit_count()
                    for a in range(8):
                        mass=chance*z**(images[a]^x).bit_count()
                        if not a:exact[j][a][feedback]+=mass
                        else:
                            exact[j][a][a^feedback]+=mass*Q(1,2**updates)
                            for fresh in range(1,8):
                                exact[j][a][fresh^feedback]+=mass*(1-Q(1,2**updates))/7
                for j in range(3):
                    self.assertGreaterEqual(endpoint(caps[j]),max(feedbacks[j][1:]))
                for cutoff in (1,2,3):
                    matrices=density.candidate(local,caps,3,cutoff)
                    bounds=[[[endpoint(m[i,k]) for k in range(m.ncols())]
                             for i in range(m.nrows())] for m in matrices]
                    for sequence in ((0,1,2),(2,0,1),(1,1,1),(2,2,2),(0,0,0)):
                        truth=[Q(1)]*8;upper=[Q(1)]*len(bounds[0])
                        for j in sequence:
                            truth=[sum(a*v for a,v in zip(row,truth)) for row in exact[j]]
                            upper=[sum(a*v for a,v in zip(row,upper)) for row in bounds[j]]
                            self.assertGreaterEqual(upper[0],truth[0])
                            self.assertGreaterEqual(upper[1],max(truth[1:]))
                            self.assertGreaterEqual(upper[2],sum(truth[1:])/7)
                            for i,level in enumerate(data['birth_class_levels'],3):
                                self.assertGreaterEqual(upper[i],max(truth[a] for a in range(1,8)
                                                                    if images[a].bit_count()==level))

    def test_invalid_inputs(self):
        data=density.base.prepare(list(range(8)),[1,2,4,3,5,7,6,1],3)
        for z in (arb(0),arb(2)):
            with self.assertRaises(ValueError):density.feedback_caps(data,z)
        local=density.base.outward_at_z(data,arb(3)/4)
        caps=density.feedback_caps(data,arb(3)/4)
        for cutoff in (0,4,True,1.5):
            with self.assertRaises(ValueError):density.candidate(local,caps,3,cutoff)


if __name__=='__main__':unittest.main()
