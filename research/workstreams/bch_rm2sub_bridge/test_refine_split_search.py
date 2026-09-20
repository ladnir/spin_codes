import unittest

import refine_split_search as refine


class ProposalTests(unittest.TestCase):
    def test_complete_has_no_proposal(self):
        self.assertIsNone(refine.propose([[-80]*201], 200, 200, []))

    def test_largest_gap_and_worst_case(self):
        powers = [[-80]*(q+1) for q in range(150, 513)]
        for q in range(210, 252):
            powers[q-150][3] = 100
        powers[38][0] = 100
        tilt, q, h = refine.propose(powers, 150, 512, [])
        self.assertEqual((q, h), (230, 3))
        shards = [dict(tilt=tilt, anchor_occupation=q, anchor_all_one_rows=h)]
        self.assertNotEqual(refine.propose(powers, 150, 512, shards), (tilt, q, h))

    def test_shape_gap_rejected(self):
        with self.assertRaises(ValueError):
            refine.propose([[-80]*200], 200, 200, [])


if __name__ == '__main__':
    unittest.main()
