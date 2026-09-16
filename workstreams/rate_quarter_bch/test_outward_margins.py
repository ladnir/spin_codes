"""Exact toy checks and retained-certificate regressions."""
from fractions import Fraction as F
import itertools
import json
import math
import unittest
from unittest.mock import patch

import numpy as np
import certify_smaller_margins as c


def exact_epoch(t,s,spectrum,kernel,z,maximum):
    rows=[]; den=(1<<s)-1; kappa=F(den,den-1)
    for j in range(maximum+1):
        total=math.comb(t,j)
        moments={w:sum((math.comb(w,v)*math.comb(t-w,j-v)*z**(w+j-2*v)
                        for v in range(max(0,j-t+w),min(w,j)+1)),F(0))/total for w in spectrum}
        arbitrary=max(moments.values()); uniform=sum((spectrum[w]*v for w,v in moments.items()),F(0))/den
        nonkernel=F(total-kernel[j],total)
        rows.append((F(kernel[j],total)*z**j,nonkernel*z**j,F(0),
                     min(nonkernel,arbitrary)/den,F(0),arbitrary,
                     kappa*min(nonkernel,uniform)/den,F(0),kappa*uniform))
    return rows


def exact_mul(a,b):
    return tuple(sum((a[3*i+k]*b[3*k+j] for k in range(3)),F(0)) for i in range(3) for j in range(3))


class OutwardMarginTests(unittest.TestCase):
    def setUp(self):
        self.precision=c.ctx.prec; c.ctx.prec=256

    def tearDown(self):
        c.ctx.prec=self.precision

    def assert_upper(self,got,expected):
        # Serialization provides ample outward slack except at exact dyadics.
        bound=c.unpack(c.pack(got))
        self.assertTrue(bound>=c.number(expected).upper())

    def test_dyadic_serialization_rounds_up(self):
        for x in (c.arb(0),c.arb(1)/3,c.arb(2)**-1000000,(c.arb(7)/3).exp()):
            self.assertTrue(c.unpack(c.pack(x))>=x.upper())

    def test_epoch_against_exact_rationals(self):
        spectrum={2:2,4:1}; kernel=[1,0,2,0,1]; z=F(3,4)
        with patch.multiple(c,T=4,S=2):
            actual=c.epoch(spectrum,kernel,-c.number(z).log(),4)
            expected=exact_epoch(4,2,spectrum,kernel,z,4)
        for a,b in zip(actual,expected):
            for x,y in zip(a,b): self.assert_upper(x,y)

    def test_region_against_exact_two_epoch_convolution(self):
        spectrum={2:2,4:1}; kernel=[1,0,2,0,1]; z=F(3,4)
        expected=exact_epoch(4,2,spectrum,kernel,z,2)
        with patch.multiple(c,T=4,S=2,L=8):
            epoch=c.epoch(spectrum,kernel,-c.number(z).log(),2)
            binary=c.region(epoch,2); linear=c.region(epoch,2,linear=True)
        for j in range(3):
            exact=[F(0)]*9
            for v in range(j+1):
                product=exact_mul(expected[v],expected[j-v])
                for k,x in enumerate(product): exact[k]+=math.comb(4,v)*math.comb(4,j-v)*x/math.comb(8,j)
            for a,b,x in zip(binary[j],linear[j],exact):
                self.assert_upper(a,x); self.assert_upper(b,x)

    def test_q1_against_explicit_supports(self):
        z=F(3,4); spectrum={2:2,4:1}; den=3; kappa=F(3,2)
        m0=sum((n*z**w for w,n in spectrum.items()),F(0))/den
        m1=sum((n*(w*z**(w-1)+(4-w)*z**(w+1)) for w,n in spectrum.items()),F(0))/(4*den)
        zero=(F(1),F(0),F(0),F(0),F(0),z**2,F(0),F(0),kappa*m0)
        one=(F(0),z,F(0),z/den,F(0),z,kappa*m1/den,F(0),kappa*m1)
        rz=exact_mul(zero,zero)
        ra=tuple((a+b)/2 for a,b in zip(exact_mul(zero,one),exact_mul(one,zero)))
        expected=[]
        for w in range(5):
            total=F(0)
            for support in itertools.combinations(range(4),w):
                state=tuple(F(int(i==j)) for i in range(3) for j in range(3))
                for j in range(4): state=exact_mul(state,ra if j in support else rz)
                total+=sum(state[:3])
            expected.append(total/math.comb(4,w))
        with patch.multiple(c,T=4,S=2,L=8,B=4):
            binary=c.q1(spectrum,-c.number(z).log()); linear=c.q1(spectrum,-c.number(z).log(),linear=True)
        for a,b,x in zip(binary,linear,expected):
            self.assert_upper(a,x); self.assert_upper(b,x)

    def test_real_epoch_matches_numerical_reference(self):
        spectrum,kernel,_=c.fixed.load_inner()
        for z in (-8.,-.5):
            actual=c.epoch(spectrum,kernel,c.number(z).exp(),128)
            logs=np.array([[float(v.log()) if v>0 else -np.inf for v in row] for row in actual]).reshape(129,3,3)
            reference=c.fixed.general.epoch_logs(128,19,spectrum,kernel,math.exp(z),128)
            np.testing.assert_allclose(logs,reference,atol=3e-11,rtol=3e-13)

    def test_certificate_and_replay_receipts(self):
        path=c.fixed.HERE/'SMALLER_MARGIN_CERTIFICATE.json'
        certificate=json.loads(path.read_text())
        replay=json.loads((c.fixed.HERE/'SMALLER_MARGIN_CERTIFICATE_REPLAY.json').read_text())
        self.assertEqual(replay['producer_sha256'],c.fixed.sha(path))
        self.assertGreater(replay['precision_bits'],certificate['precision_bits'])
        self.assertEqual(replay['source_sha256'],certificate['source_sha256'])
        for name,digest in certificate['source_sha256'].items():
            self.assertEqual(c.fixed.sha(c.fixed.ROOT/name),digest,name)
        witnesses=json.loads((c.fixed.HERE/'SMALLER_OUTWARD_WITNESSES.json').read_text())
        for row,witness in zip(certificate['results'],witnesses['results']):
            c.check_coverage(witness['dense']['selected_boxes'],c.L,row['maximum_sparse']+1)
            total=c.unpack(row['union_upper'])
            self.assertTrue(total<c.arb(2)**(-row['target_bits']))
            self.assertTrue(sum((c.unpack(v) for v in row['group_upper'].values()),c.arb(0))<=total)
            for group,value in row['group_upper'].items():
                self.assertTrue(sum((c.unpack(v) for v in row[group]),c.arb(0))<=c.unpack(value))


if __name__=='__main__': unittest.main()
