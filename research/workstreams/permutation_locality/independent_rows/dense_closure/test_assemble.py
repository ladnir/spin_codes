from fractions import Fraction as Q
import unittest
from assemble import sparse_upper,SPARSE_UPPER,SPARSE_EXTENSION


class HandoffTests(unittest.TestCase):
    def test_only_verified_ranges_and_cutoffs_are_accepted(self):
        self.assertEqual(sparse_upper(59,209715),SPARSE_UPPER)
        self.assertEqual(sparse_upper(65,104857),SPARSE_UPPER+SPARSE_EXTENSION)
        for q,d in ((59,209716),(65,104858),(66,10485),(64,10485),(65,-1),(65,None),(65.0,10485)):
            with self.assertRaises(ValueError):sparse_upper(q,d)

    def test_extension_decimal_sum_is_conservative(self):
        values=['5.257711e-319','1.398782e-252','6.904732e-184',
                '5.411004e-113','5.78982038e-40','1.040080e-104']
        self.assertLess(sum(map(Q,values)),SPARSE_EXTENSION)
        self.assertLess(sparse_upper(65,104857),Q(1,2**43))


if __name__=='__main__':unittest.main()
