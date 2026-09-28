from fractions import Fraction as Q
from unittest.mock import patch
import unittest

import numpy as np

from flint import arb,arb_mat,ctx

import density_tangent
from heterogeneous import density_loss,capped_density_loss
from mixture_probe import outward,counts_from_items,posterior_density_anchors
from probe import aq


class CompositionProbeTests(unittest.TestCase):
    def setUp(self):ctx.prec=256

    def test_explicit_composition_requires_all_counts(self):
        components=[('0,0',),('0,3',),('4,4',)]
        self.assertEqual(counts_from_items(components,['0,0=3695','0,3=400','4,4=1']),[3695,400,1])
        for items in (['0,3=400'],['0,0=4096','0,0=0'],['unknown=4096'],['0,0=4097','0,3=-1']):
            with self.assertRaises(ValueError):counts_from_items(components,items)

    def test_composition_caps_have_the_expected_exact_factor(self):
        components=[('zero',Q(1),Q(1),Q(0),Q(0),0),
                    ('active',Q(3),Q(1,4),Q(1,2),Q(1,4),1)]
        counts=[4095,1];witness=[Q(1,100),Q(1,4),Q(1,2)]
        with patch('mixture_probe.iid_kernel.outward',return_value=arb_mat([[1,0],[0,1]])):
            old=outward(components,None,counts,witness,threshold=0)
            new=outward(components,None,counts,witness,threshold=0,capped=True)
        expected=aq(capped_density_loss(4096,(4096,1,1))/density_loss(4096))**256
        self.assertTrue(abs(new/old/expected-1)<arb(2)**-180)
        self.assertTrue(new<old)

    def test_density_tangent_uses_unnormalized_packet_weights(self):
        components=[('zero',Q(1),Q(1),Q(0),Q(0),0),
                    ('active',Q(3),Q(1,4),Q(1,2),Q(1,4),1)]
        counts=[4095,1];witness=[Q(1,100),Q(1,4),Q(1,2)];anchors=[4094,1,1]
        with patch('mixture_probe.iid_kernel.outward',return_value=arb_mat([[1,0],[0,1]])):
            result=outward(components,None,counts,witness,threshold=0,density_anchors=anchors)
        logB,logtau=density_tangent.outward(4096,anchors)
        Z=Q(1,4)+Q(1,2)*witness[1]+Q(1,4)*witness[2]
        weights=[Q(4095,4096)+Q(1,4)/Z/4096,Q(1,2)/Z/4096,Q(1,4)/Z/4096]
        S=sum((aq(w)*t.exp() for w,t in zip(weights,logtau)),arb(0))
        expected=(4096*3)*aq(Z)**256*(256*logB).exp()*S**(4096*256)
        self.assertTrue(abs(result/expected-1)<arb(2)**-180)

    def test_posterior_anchor_proposal_for_constant_kernel(self):
        probabilities=np.array([[1.,0.,0.],[.25,.5,.25]])
        counts=[3695,401];witness=[.1,.25,.5];anchors=[3900,195,1]
        prepared=(None,probabilities,None,None)
        with patch('mixture_probe.iid_kernel.floating',return_value=np.eye(2)):
            proposed=posterior_density_anchors(prepared,counts,['1/10','1/4','1/2'],anchors)
        z=np.array([1.,witness[1],witness[2]])
        weights=(np.array(counts)/4096/(probabilities@z))@probabilities
        _,logtau=density_tangent.floating(4096,anchors)
        weights*=np.exp(logtau)
        expected=tuple(int(round(4096*w/weights.sum())) for w in weights)
        self.assertEqual(proposed,expected)
        self.assertTrue(all(type(k) is int and 0<=k<=4096 for k in proposed))

    def test_posterior_anchor_preserves_absent_category(self):
        prepared=(None,np.array([[1.,0.,0.]]),None,None)
        with patch('mixture_probe.iid_kernel.floating',return_value=np.eye(2)):
            self.assertEqual(posterior_density_anchors(prepared,[4096],[.1,1,1],[4000,50,46]),(4096,0,0))


if __name__=='__main__':unittest.main()
