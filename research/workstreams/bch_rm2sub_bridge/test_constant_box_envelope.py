import unittest
import numpy as np
from screen_larger_constant_boxes import box_envelope


class BoxEnvelopeTests(unittest.TestCase):
    def test_every_entry_and_truncated_endpoint(self):
        # Full length exercises h+j>8192, which must be excluded, not wrapped.
        index = np.arange(8193,dtype=float)
        region = np.stack([np.sin(index/(k+1)) for k in range(9)],axis=1).reshape(-1,3,3)
        for lo,hi in ((0,0),(0,7),(13,29),(8170,8192)):
            maximum = 8192-lo
            actual = box_envelope(region,lo,hi,maximum)
            self.assertEqual(actual.shape,(maximum+1,3,3))
            for j in range(maximum+1):
                expected = region[lo+j:min(hi+j,8192)+1].max(axis=0)
                np.testing.assert_array_equal(actual[j],expected)


if __name__ == '__main__':unittest.main()
