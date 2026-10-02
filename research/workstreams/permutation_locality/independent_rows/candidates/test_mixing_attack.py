import unittest
from fractions import Fraction as Q
from flint import arb,arb_mat,ctx
import mixing_attack as candidate
from test_fresh_history import endpoint


class MixingAttackTests(unittest.TestCase):
    def test_retarget_abstract_lazy_refresh_components(self):
        ctx.prec=192
        spectrum={48:3,56:7,64:11,72:13,80:17};states=sum(spectrum.values())
        a=Q(1,4);f=Q(3,4);z=Q(1,7)
        old=arb_mat(11,11)
        for source in (candidate.F,candidate.M,*range(candidate.U,candidate.U+5),9,10):
            old[source,candidate.M]=arb(3)/16
            old[source,candidate.Z]=arb(9)/(16*states)
            if source==candidate.F or candidate.U<=source<candidate.U+5:old[source,candidate.Z]+=arb(1)/28
            for i,v in enumerate(spectrum):old[source,candidate.U+i]=arb(9)*spectrum[v]/(16*states)
        old[candidate.C,candidate.C]=arb(1)/28
        old[candidate.C,candidate.Z]=arb(1)/28
        for rounds in (2,3,8,None):
            b=Q(0) if rounds is None else Q(1,2**rounds)
            new=candidate.retarget_matrix(old,1,spectrum,'.056',2,rounds)
            for source in (candidate.F,*range(candidate.U,candidate.U+5)):
                self.assertGreaterEqual(endpoint(new[source,candidate.Z]),b*z+(1-b)*f/states)
            self.assertGreaterEqual(endpoint(new[candidate.M,candidate.Z]),(1-b)*f/states)
            self.assertGreaterEqual(endpoint(new[candidate.C,candidate.C]),b*z)

    def test_reject_weaker_mixing(self):
        with self.assertRaises(ValueError):candidate.retarget_matrix(arb_mat(11,11),1,{},'.05',2,1)

    def test_empty_uniform_self_loop_is_not_scaled_as_pure_refresh(self):
        ctx.prec=192
        spectrum={48:3,56:7,64:11,72:13,80:17};states=sum(spectrum.values())
        tilt=(arb(8)/7).log();old=arb_mat(11,11)
        for i,v in enumerate(spectrum):
            f=(arb(7)/8)**v
            for j,w in enumerate(spectrum):old[candidate.U+i,candidate.U+j]=f*3*spectrum[w]/(4*states)
            old[candidate.U+i,candidate.U+i]+=f/4
        for rounds in (3,8,None):
            a=Q(0) if rounds is None else Q(1,2**rounds)
            new=candidate.retarget_matrix(old,0,spectrum,tilt,2,rounds)
            for i,v in enumerate(spectrum):
                for j,w in enumerate(spectrum):
                    expected=Q(7,8)**v*((a if i==j else 0)+(1-a)*spectrum[w]/states)
                    self.assertGreaterEqual(endpoint(new[candidate.U+i,candidate.U+j]),expected)


if __name__=='__main__':unittest.main()
