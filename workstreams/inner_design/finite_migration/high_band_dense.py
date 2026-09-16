"""Separate assignments containing a reference band above 3/4."""
import argparse
from fractions import Fraction as F
from pathlib import Path
from types import SimpleNamespace

from flint import arb,ctx
import zero_constant_dense as parent

model = parent.model
CUTOFF = F(3,4)
XI = F(-5123,64)


class Checker(parent.Checker):
    def _constant_comparison(self,lo,hi,a,b,precision):
        assert precision == ctx.prec
        nlo,nhi = self.density(a,b)
        ps = tuple(p for p in self.ps if p <= CUTOFF)
        return parent.parent.constant_density.factor(self.rows,lo,hi,nlo,nhi,ps,F(0))

    def split_bound(self,lo,hi,a,b,witness):
        nlo,nhi = self.density(a,b)
        if nhi == 1:
            return None
        old = self.fixed_input_bound(lo,hi,a,b,witness)
        if old is None:
            return None
        theta = (F(lo,self.rows)*nlo+F(hi,self.rows)*nhi)/2
        r,eta = (model.base.decode(witness[k]) for k in ('fixed_r','eta'))
        x = theta*(1-r)/(r*(1-theta))
        scalar_module = parent.parent.fixed_input.scalar
        costs = scalar_module.tilted_costs(self.bands,self.ps,self.caps,model.number(x).log())
        low,high = arb(0),arb(0)
        for cost,p in zip(costs,self.ps+(F(1),)):
            term = (cost-model.number(eta*p)).exp()
            if p <= CUTOFF:
                low += term
            else:
                high += term
        def scalar(total):
            return scalar_module.scalar_bound(self.rows,lo,hi,nlo,nhi,total.log().upper(),eta,model.number(x))
        original_scalar = scalar(low+high)
        exceptional = model.up(old+scalar(low+high*(-model.number(XI)).exp())-original_scalar+model.number(XI))
        if nlo > max(p for p in self.ps if p <= CUTOFF):
            return exceptional
        ratio = self.constant_comparison(lo,hi,a,b,ctx.prec)
        original = self.comparison(lo,hi,a,b,ctx.prec)
        ordinary = model.up(old+scalar(low)-original_scalar+256*(ratio.log()-original.log()))
        return model.up((ordinary.exp()+exceptional.exp()).log())

    def bound(self,lo,hi,a,b,witness):
        # Keep the already-accepted unsplit/zero-only alternatives as whole
        # bounds; this class does not replace their matrix entries.
        old = parent.parent.short_dense.Checker.bound(self,lo,hi,a,b,witness)
        if old < -80*arb(2).log():
            return old
        improved = self.split_bound(lo,hi,a,b,witness)
        return old if improved is None else min(old,improved)


def run(*args):
    original = parent.parent.budget_dense.short_dense
    parent.parent.budget_dense.short_dense = SimpleNamespace(Checker=Checker,mixed_dense=parent.parent.short_dense.mixed_dense)
    try:
        return parent.parent.budget_dense.run(*args)
    finally:
        parent.parent.budget_dense.short_dense = original


if __name__ == '__main__':
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--output',type=Path,required=True)
    p.add_argument('--seed',type=Path,required=True)
    p.add_argument('--nodes',type=int,default=300)
    p.add_argument('--seconds',type=float,default=600)
    p.add_argument('--verify',action='store_true')
    a = p.parse_args()
    run(a.output.resolve(),16,64,a.seed.resolve(),a.nodes,a.seconds,42,a.verify)
