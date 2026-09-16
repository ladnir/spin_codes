"""Complete subset bounds with the equivalent pruned routing calculation."""
import argparse
from fractions import Fraction as F
from pathlib import Path
from types import SimpleNamespace

from flint import arb,ctx
import fast_density
import high_band_dense
import zero_constant_dense
import poisson_density_factor as old_density


def unrestricted(checker,lo,hi,a,b,precision):
    assert precision == ctx.prec
    nlo,nhi = checker.density(a,b)
    old = arb(old_density.density_factor(checker.rows))
    if nhi == 1:
        return old
    pmin = min(checker.ps)
    # The mean constraint itself implies this all-one fraction bound.
    # Adding it to the variance program excludes no feasible assignment.
    maximum = (nhi-pmin)/(1-pmin)
    return min(old,fast_density.factor(checker.rows,lo,hi,nlo,nhi,checker.ps,maximum))


class ZeroChecker(zero_constant_dense.Checker):
    def _constant_comparison(self,lo,hi,a,b,precision):
        assert precision == ctx.prec
        nlo,nhi = self.density(a,b)
        return fast_density.factor(self.rows,lo,hi,nlo,nhi,self.ps,F(0))


class Checker(high_band_dense.Checker):
    def __init__(self,exponent):
        super().__init__(exponent)
        self.zero_checker = ZeroChecker(exponent)
        assert self.zero_checker.engine.identity() == self.engine.identity()
        self.zero_checker.comparison = self.comparison

    def _comparison(self,lo,hi,a,b,precision):
        return unrestricted(self,lo,hi,a,b,precision)

    def _constant_comparison(self,lo,hi,a,b,precision):
        assert precision == ctx.prec
        nlo,nhi = self.density(a,b)
        ps = tuple(p for p in self.ps if p <= high_band_dense.CUTOFF)
        return fast_density.factor(self.rows,lo,hi,nlo,nhi,ps,F(0))

    def bound(self,lo,hi,a,b,witness):
        high = super().bound(lo,hi,a,b,witness)
        if high < -80*arb(2).log():
            return high
        return min(high,self.zero_checker.bound(lo,hi,a,b,witness))


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
    p.add_argument('--nodes',type=int,default=600)
    p.add_argument('--seconds',type=float,default=900)
    p.add_argument('--verify',action='store_true')
    a = p.parse_args()
    run(a.output.resolve(),16,64,a.seed.resolve(),a.nodes,a.seconds,42,a.verify)
