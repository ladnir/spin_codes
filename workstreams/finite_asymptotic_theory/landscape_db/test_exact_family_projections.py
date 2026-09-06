import math
import unittest

from extrapolate_exact_families import linear,project,rm_minimum


class ProjectionTest(unittest.TestCase):
    def test_rm_minimum_shell_counts(self):
        self.assertEqual(rm_minimum(8),(4,14))
        self.assertEqual(rm_minimum(32),(8,620))
        self.assertEqual(rm_minimum(2048)[0],64)
        with self.assertRaises(ValueError):rm_minimum(1024)

    def test_linear_fit_and_exact_counting_cost_restoration(self):
        self.assertAlmostEqual(linear([2,3,5],[7,9,13],9),21.)
        anchors=[]
        for block in (8,32,128,512):
            distance,count=rm_minimum(block)
            anchors.append(dict(block_bits=block,distance=distance,log_minimum_count=math.log2(count),
                margin_bits=2*distance+7-math.log2(count)-20+math.log2(block//2)))
        result=project(anchors,'rm',2048,20,4)
        distance,count=rm_minimum(2048)
        self.assertAlmostEqual(result['margin_bits'],2*distance+7-math.log2(count)-20+10)


if __name__=='__main__':unittest.main()
