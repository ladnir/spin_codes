"""Coverage and exact-endpoint validation for the low-occupancy driver."""
from copy import deepcopy
from fractions import Fraction as F
from pathlib import Path
import unittest
from unittest.mock import patch

import low_complete as low


class LowCompleteTests(unittest.TestCase):
    def fixture(self):
        return dict(occupancy_covered=[3,4,5], occupancy_uppers={str(q): [q,-100] for q in (3,4,5)},
            whole_code_certificate=False, K=1048576, N=2097152, threshold=209715,
            groups=8192, state_bits=22, map_record=dict(map_sha256=low.MAP))

    def test_exact_positive_dyadic_and_consecutive_coverage(self):
        self.assertEqual(low.tail_endpoints(self.fixture(), 5), {q:F(q, 1 << 100) for q in (3,4,5)})
        for endpoint in ([0,-100], [-1,2], [True,-10], [1,1.5], [1]):
            with self.subTest(endpoint=endpoint), self.assertRaises(ValueError):
                low.fraction(endpoint)

    def test_missing_occupancy_or_wrong_geometry_is_rejected(self):
        changes = [dict(K=262144), dict(N=524288), dict(threshold=52428), dict(groups=2048),
                   dict(state_bits=16), dict(whole_code_certificate=True),
                   dict(occupancy_covered=[3,5]), dict(occupancy_uppers={'3':[1,-100],'5':[1,-100]}),
                   dict(map_record=dict(map_sha256='other'))]
        for change in changes:
            record = self.fixture()
            record.update(change)
            with self.subTest(change=change), self.assertRaises(ValueError):
                low.tail_endpoints(record, 5)

    def test_invalid_range_and_existing_path_fail_before_census(self):
        with patch.object(low.maps, 'prepare') as prepare:
            for maximum in (2,129,True):
                with self.subTest(maximum=maximum), self.assertRaises(ValueError):
                    low.run(Path('unused'), maximum=maximum)
            with patch.object(Path, 'exists', return_value=True), self.assertRaises(ValueError):
                low.run(Path('unused'))
            prepare.assert_not_called()


if __name__ == '__main__':
    unittest.main()
