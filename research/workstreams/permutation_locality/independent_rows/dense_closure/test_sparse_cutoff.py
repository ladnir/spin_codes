"""The adjustable cutoff changes only the final Chernoff factor."""
from contextlib import redirect_stdout
from fractions import Fraction as Q
from io import StringIO
from types import SimpleNamespace
import unittest

import numpy as np
from flint import arb,arb_mat,ctx
import kernel  # Adds the shared proof directory to sys.path.
from occupancy_cdf_cover import cover


class SparseCutoffTests(unittest.TestCase):
    def test_replay_cutoff_factor_and_legacy_default(self):
        ctx.prec=192
        args=SimpleNamespace(groups=1,joint_witness=False,joint_top=1,
            probe_supports=[],probe_vector=[],retain_parents=True,
            target_bits=40,max_splits=0,screen_only=False,precision=192)
        # A synthetic tiny region moment lets every cover pass. There is
        # one operator/tilt, so the chosen witnesses do not change with D.
        tiny=arb(2)**-20
        operators={('0.001','1'):([arb_mat([[tiny]])]*2,
                                  [np.array([[2.**-20]])]*2)}
        counts={'1':[0]*38+[1]*219}
        def run(**kwargs):
            with redirect_stdout(StringIO()):return cover(args,operators,counts,[1],**kwargs)
        low=run(cutoff=10000);high=run(cutoff=20000)
        # The routine returns rounded-up dyadic bounds, not enclosures of
        # the exact expression, so compare logs within a tight tolerance.
        self.assertAlmostEqual(float((high/low).log()),10.,places=12)
        self.assertTrue(run().overlaps(run(cutoff=209715)))
        for bad in (-1,1<<21,Q(1,2),True):
            with self.assertRaises(ValueError):run(cutoff=bad)


if __name__=='__main__':unittest.main()
