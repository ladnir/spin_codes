import unittest
from fractions import Fraction as F
from flint import arb,ctx
import hd_lambda_shear_a12_group_range as gate


class GroupRangeTests(unittest.TestCase):
    def test_r18_formula_matches_frozen_controller(self):
        with ctx.workprec(256):
            raw,reference=arb(2)**817,arb(-7000)
            left=gate.log_bound(raw,18,reference)
            right=gate.IDEAL.method_expression(raw,4,reference)
            self.assertTrue((left-right).contains(0))

    def test_group_factor_is_linear(self):
        with ctx.workprec(256):
            difference=gate.log_bound(arb(2)**30,20,arb(-7000))-gate.log_bound(arb(2)**30,19,arb(-7000))
            self.assertTrue((difference-arb(8)).contains(0))

    def test_source_union_not_theta_scaled(self):
        with ctx.workprec(256):
            difference=gate.log_bound(arb(2)**30,18,arb(-6999))-gate.log_bound(arb(2)**30,18,arb(-7000))
            self.assertTrue((difference-arb(1)).contains(0))

    def test_invalid_group_lists(self):
        for groups in ((),(0,),(33,),(18,18),(19,18),(True,),(18.0,)):
            with self.assertRaises(ValueError): gate.checked_groups(groups)

    def test_valid_groups(self):
        self.assertEqual(gate.checked_groups([1,18,32]),(1,18,32))

    def test_bad_raw_weight(self):
        for raw in (arb(0),arb(-1)):
            with self.assertRaises(ValueError): gate.log_bound(raw,18,arb(-7000))


if __name__=='__main__': unittest.main()
