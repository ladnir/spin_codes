import unittest
from itertools import product
from math import log2

import numpy as np

from hd_grouped_packet_screen import log_multinomial
from hd_shifted_box_screen import histogram_ranges,support_count_markers,coefficient_envelope
from hd_shifted_general_count import ProfileBox,joint_is_feasible,matching_histograms,box_count_upper,exact_joint_count
from hd_shifted_general_sweep import CHECKPOINTS
from hd_shifted_good_pair_screen import regional


class BoxScreenTests(unittest.TestCase):
    def test_histogram_ranges_cover_feasible_points(self):
        box = ProfileBox((1,3),(0,3),(0,2),(0,3))
        ranges = histogram_ranges(box,pairs=4)
        for profile in product(*(range(lo,hi+1) for lo,hi in box.ranges)):
            if joint_is_feasible(4,*profile):
                for hist,bounds in zip(matching_histograms(4,*profile),ranges):
                    self.assertTrue(all(lo<=n<=hi for n,(lo,hi) in zip(hist,bounds)))

    def test_proposed_count_tilt_remains_exact_upper_bound(self):
        box = ProfileBox((1,3),(0,3),(0,2),(0,3))
        markers = support_count_markers(box,pairs=4)
        self.assertTrue(all(marker>0 for marker in markers))
        count = sum(exact_joint_count(4,*profile)
                    for profile in product(*(range(lo,hi+1) for lo,hi in box.ranges)))
        self.assertGreaterEqual(box_count_upper(4,box,markers),count)

    def test_point_envelope_reduces_to_existing_region(self):
        profile,z,alpha,odds = CHECKPOINTS['asymmetric']
        hist = matching_histograms(4096,*profile)[0]
        actual,_ = coefficient_envelope(z,alpha,odds[0],tuple((n,n) for n in hist[1:]))
        expected,_ = regional(z,alpha,hist,odds[0])
        expected += log_multinomial(hist)+np.dot(hist[1:],odds[0])
        np.testing.assert_allclose(actual,expected,rtol=0,atol=1e-10)

    def test_neighborhood_contains_point_majorant(self):
        profile,z,alpha,odds = CHECKPOINTS['asymmetric']
        hist = matching_histograms(4096,*profile)[0]
        point,_ = coefficient_envelope(z,alpha,odds[0],tuple((n,n) for n in hist[1:]))
        box = ProfileBox(*((n-1,n+1) for n in profile))
        neighborhood,_ = coefficient_envelope(z,alpha,odds[0],histogram_ranges(box)[0][1:])
        self.assertTrue(np.all(neighborhood>=point-1e-10))


if __name__=='__main__':
    unittest.main()
