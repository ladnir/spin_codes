"""Split exactly zero versus at least one all-one outer row."""
import argparse
from fractions import Fraction as F
from pathlib import Path
from types import SimpleNamespace

from flint import arb, ctx
import constant_dense as parent

model = parent.model


class Checker(parent.Checker):
    def _constant_comparison(self, lo, hi, a, b, precision):
        assert precision == ctx.prec
        nlo, nhi = self.density(a, b)
        return parent.constant_density.factor(self.rows, lo, hi, nlo, nhi, self.ps, F(0))

    def split_bound(self, lo, hi, a, b, witness):
        nlo, nhi = self.density(a, b)
        if nhi == 1:
            return None
        old = self.fixed_input_bound(lo, hi, a, b, witness)
        if old is None:
            return None
        theta = (F(lo, self.rows)*nlo+F(hi, self.rows)*nhi)/2
        r, eta = (model.base.decode(witness[k]) for k in ('fixed_r', 'eta'))
        x = theta*(1-r)/(r*(1-theta))
        costs = parent.fixed_input.scalar.tilted_costs(self.bands, self.ps, self.caps, model.number(x).log())
        ordinary_sum = sum(((g-model.number(eta*p)).exp()
                            for g,p in zip(costs[:-1], self.ps)), arb(0))
        constant = (costs[-1]-model.number(eta)).exp()
        def scalar(total):
            return parent.fixed_input.scalar.scalar_bound(
                self.rows, lo, hi, nlo, nhi, total.log().upper(), eta, model.number(x))
        original_scalar = scalar(ordinary_sum+constant)
        # h is an integer and xi<=0: xi*h<=xi on h>=1.
        assert parent.XI < 0
        exceptional = model.up(old+scalar(ordinary_sum+constant*(-model.number(parent.XI)).exp())
                               -original_scalar+model.number(parent.XI))
        if nlo > max(self.ps):
            return exceptional
        comparison = self.constant_comparison(lo, hi, a, b, ctx.prec)
        original = self.comparison(lo, hi, a, b, ctx.prec)
        ordinary = model.up(old+scalar(ordinary_sum)-original_scalar
                            +256*(comparison.log()-original.log()))
        return model.up((ordinary.exp()+exceptional.exp()).log())


def run(*args):
    original = parent.budget_dense.short_dense
    parent.budget_dense.short_dense = SimpleNamespace(Checker=Checker, mixed_dense=parent.short_dense.mixed_dense)
    try:
        return parent.budget_dense.run(*args)
    finally:
        parent.budget_dense.short_dense = original


if __name__ == '__main__':
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--output', type=Path, required=True)
    p.add_argument('--seed', type=Path, required=True)
    p.add_argument('--nodes', type=int, default=400)
    p.add_argument('--seconds', type=float, default=900)
    p.add_argument('--verify', action='store_true')
    a = p.parse_args()
    run(a.output.resolve(), 16, 64, a.seed.resolve(), a.nodes, a.seconds, 42, a.verify)
