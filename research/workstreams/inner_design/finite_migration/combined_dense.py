"""Retain both complete subset bounds: all-one-only and high-band splits."""
import argparse
from pathlib import Path
from types import SimpleNamespace

from flint import arb
import high_band_dense
import zero_constant_dense


class Checker(high_band_dense.Checker):
    def __init__(self,exponent):
        super().__init__(exponent)
        self.zero_checker = zero_constant_dense.Checker(exponent)
        assert self.zero_checker.engine.identity() == self.engine.identity()
        assert self.zero_checker.ps == self.ps and self.zero_checker.rows == self.rows
        # Only the unrestricted routing comparison is shared. The two
        # constrained variance caches have different allowed references.
        self.zero_checker.comparison = self.comparison

    def bound(self,lo,hi,a,b,witness):
        high = super().bound(lo,hi,a,b,witness)
        if high < -80*arb(2).log():
            return high
        zero = self.zero_checker.bound(lo,hi,a,b,witness)
        return min(high,zero)


def run(*args):
    driver = high_band_dense.parent.parent.budget_dense
    original = driver.short_dense
    driver.short_dense = SimpleNamespace(Checker=Checker,
        mixed_dense=high_band_dense.parent.parent.short_dense.mixed_dense)
    try:
        return driver.run(*args)
    finally:
        driver.short_dense = original


if __name__ == '__main__':
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--output',type=Path,required=True)
    p.add_argument('--seed',type=Path,required=True)
    p.add_argument('--nodes',type=int,default=300)
    p.add_argument('--seconds',type=float,default=600)
    p.add_argument('--verify',action='store_true')
    a = p.parse_args()
    run(a.output.resolve(),16,64,a.seed.resolve(),a.nodes,a.seconds,42,a.verify)
