"""Regression checks for the new density checker and refresh certificates."""
from fractions import Fraction as F
import math
import unittest
from flint import arb,ctx
import frontier_dense as old
import frontier_dense_density as new
import frontier_ledger as ledger
import frontier_sparse as sparse
import bridge as base
from poisson_density_factor import density_factor


class DensityFrontierTest(unittest.TestCase):
    def test_count_cost_replaces_only_density_factor(self):
        ctx.prec=256
        for m in (20,26,28,30):
            checker=object.__new__(new.Checker);checker.rows=1<<(m-7)
            difference=256*(arb(checker.rows+1)/density_factor(checker.rows)).log()
            for lo,hi in ((2,7),(checker.rows//4,3*checker.rows//4),(checker.rows,checker.rows)):
                previous=old.Checker._count_cost(checker,lo,hi)
                current=checker._count_cost(lo,hi)
                self.assertTrue((previous-current).overlaps(difference))
                self.assertTrue(current<previous)

    def test_bound_changes_by_the_proved_constant(self):
        ctx.prec=256
        ps,_=old.discovery.row_witnesses(old.caps_module.caps(),.75)
        a=old.Checker(28,ps);b=new.Checker(28,ps)
        expected=256*(arb(a.rows+1)/density_factor(a.rows)).log()
        for q in (8192,524288):
            lo,hi=F(1,8),F(9,64);tilt=b.witness(q,q,lo,hi)
            delta=a.bound(q,q,lo,hi,tilt)-b.bound(q,q,lo,hi,tilt)
            # Both methods take outward point endpoints after summing;
            # compare to a much wider tolerance than their 256-bit ulps.
            self.assertTrue(abs(delta-expected)<arb('1e-50'))

    def test_refresh_receipt_rounding_and_partition(self):
        data=base.read(base.HERE/'generated/frontier_k28_q2_q7_v1.json')
        self.assertEqual(data['method'],'refresh_four_state_kernel_exclusion')
        self.assertEqual([r['occupation'] for r in data['rows']],list(range(2,8)))
        for row,power in zip(data['rows'],data['upper_powers']):
            self.assertLessEqual(base.decode(row['upper']),sparse.as_fraction(power))
        self.assertEqual(base.decode(data['range_upper']),sum(map(sparse.as_fraction,data['upper_powers']),F(0)))
        ledger.validate_coverage(((1,1),(2,7),(8,6143),(6144,8191),(8192,524287),
            (524288,1572864),(1572865,2097152)),2097152)


if __name__=='__main__':unittest.main()
