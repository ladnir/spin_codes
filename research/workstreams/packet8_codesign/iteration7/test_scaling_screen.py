"""Independent small exact checks of changed-geometry argument plumbing."""
from fractions import Fraction as F
from itertools import combinations
from math import comb
from unittest import TestCase, main
from unittest.mock import patch
import numpy as np
import scaling_screen as screen


class ScalingTests(TestCase):
    def test_geometry_and_no_reset(self):
        for power, groups, physical, macro, cutoff in ((18,1024,128,32,52428),
                                                     (20,4096,512,128,209715)):
            g = screen.geometry(power)
            self.assertEqual((g['groups'],g['physical_epochs'],g['macro_epochs'],g['cutoff']),
                             (groups,physical,macro,cutoff))
            self.assertEqual(g['regions'],64)
            self.assertEqual(64*physical*64,2*(1<<power))
            self.assertTrue(g['continuous_state'])
            self.assertFalse(g['final_flush'])

    def test_parameter_forwarding(self):
        g=screen.geometry(20)
        with patch.object(screen.go,'placement',return_value='r') as call:
            self.assertEqual(screen.regional_family('m',16,g),'r')
            call.assert_called_once_with('m',16,epochs=128,windows=32)
        with patch.object(screen.go,'occupancy_upper',return_value='o') as call:
            screen.occupancy('r',3,F(3,4),F(2,5),g)
            call.assert_called_once_with('r',3,F(3,4),F(2,5),groups=4096,regions=64,cutoff=209715)
        with patch.object(screen.q1,'q1_support_upper',return_value=(F(1),F(1,8))) as call:
            with patch.object(screen.q1,'combine_supports',return_value=F(1,16)) as combine:
                screen.one_group('active',F(3,4),g)
                call.assert_called_once_with('active',F(3,4),regions=64,epochs=512,cutoff=209715)
                combine.assert_called_once_with((F(1),F(1,8)),groups=4096)

    def test_more_epochs_tiny_noncommuting_exact(self):
        local=np.array([[[1.,0.],[.125,.5]],[[.25,.5],[.125,.25]],[[.125,.25],[0.,.125]]])
        mats=[np.array([[F(float(v)) for v in row] for row in m],dtype=object) for m in local]
        result=screen.go.placement(local,3,epochs=4,windows=2)
        for q in range(4):
            exact=np.full((2,2),F(0),dtype=object)
            for locations in combinations(range(8),q):
                prod=np.eye(2,dtype=object)
                for step in range(4):
                    prod=prod@mats[sum(x//2==step for x in locations)]
                exact+=prod
            exact/=comb(8,q)
            upper=np.array([[F(float(x))*F(2)**result[q].exponent for x in row]
                            for row in result[q].value],dtype=object)
            self.assertTrue(np.all(upper>=exact))
            self.assertLess(float(np.max(upper-exact)),1e-12)

    def test_probability_endpoint_caps_and_floor(self):
        for value,expected in ((F(8),F(1)),(F(1,16),F(1,16)),(F(1,1<<1000),F(1,1<<200))):
            result=screen.probability_endpoint(screen.sp.from_fraction(value))
            self.assertGreaterEqual(result,min(value,F(1)))
            self.assertLessEqual(result,expected*F(10000000001,10000000000))

    def test_rounded_log_near_floor_cannot_decide_endpoint(self):
        value=screen.sp.Scaled(np.array([[np.nextafter(1.,np.inf)]]),-200)
        self.assertEqual(screen.go.display_margin(value),200.)
        self.assertGreaterEqual(screen.probability_endpoint(value),screen.sp.scalar_fraction(value))


if __name__=='__main__':main()
