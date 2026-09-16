"""Total-budget cover with separate bounds for rare all-one-rich labels."""
import argparse
from fractions import Fraction as F
from functools import lru_cache
import math
from pathlib import Path
from types import SimpleNamespace

from flint import arb, ctx
import budget_dense
import constant_density
import fixed_input
import retune_fixed_input
import short_dense

model = fixed_input.model
MAXIMUM = F(1, 32)
XI = F(-3335, 16)


class Checker(short_dense.Checker):
    def __init__(self, exponent):
        super().__init__(exponent)
        self.constant_comparison = lru_cache(maxsize=2048)(self._constant_comparison)

    def _constant_comparison(self, lo, hi, a, b, precision):
        assert precision == ctx.prec
        nlo, nhi = self.density(a, b)
        return constant_density.factor(self.rows, lo, hi, nlo, nhi, self.ps, MAXIMUM)

    def split_bound(self, lo, hi, a, b, witness):
        nlo, nhi = self.density(a, b)
        if nhi == 1:
            return None
        old = self.fixed_input_bound(lo, hi, a, b, witness)
        if old is None:
            return None
        # fixed_input_bound uses the midpoint of the theta endpoints.
        theta = (F(lo, self.rows)*nlo+F(hi, self.rows)*nhi)/2
        r, eta = (model.base.decode(witness[k]) for k in ('fixed_r', 'eta'))
        x = theta*(1-r)/(r*(1-theta))
        costs = fixed_input.scalar.tilted_costs(self.bands, self.ps, self.caps, model.number(x).log())
        def scalar(xi):
            log_s = sum(((g-model.number(eta*p+(xi if i == len(costs)-1 else 0))).exp()
                         for i, (g,p) in enumerate(zip(costs, self.ps+(F(1),)))), arb(0)).log().upper()
            return fixed_input.scalar.scalar_bound(self.rows, lo, hi, nlo, nhi, log_s, eta, model.number(x))
        exceptional = model.up(old+scalar(XI)-scalar(F(0))+model.number(lo*XI*MAXIMUM))
        maximum_mean = MAXIMUM+(1-MAXIMUM)*max(self.ps)
        if nlo > maximum_mean:
            return exceptional
        comparison = self.constant_comparison(lo, hi, a, b, ctx.prec)
        original = self.comparison(lo, hi, a, b, ctx.prec)
        ordinary = model.up(old+256*(comparison.log()-original.log()))
        return model.up((ordinary.exp()+exceptional.exp()).log())

    def bound(self, lo, hi, a, b, witness):
        old = super().bound(lo, hi, a, b, witness)
        if old < -80*arb(2).log():
            return old
        split = self.split_bound(lo, hi, a, b, witness)
        return old if split is None else min(old, split)

    def witness(self, lo, hi, a, b):
        best = super().witness(lo, hi, a, b)
        value = self.bound(lo, hi, a, b, best)
        if value < -80*arb(2).log() or b == 1:
            return best
        q, coordinate = (lo+hi)//2, (a+b)/2
        nu = self.density(coordinate, coordinate)[0]
        theta, r = float(F(q, self.rows)*nu), float(model.base.decode(best['fixed_r']))
        u = math.log(theta/(1-theta))-math.log(r/(1-r))
        index = best['tilt']
        for tilt in (index-2, index, index+2):
            for offset in (-0.05, -0.025, 0, 0.025, 0.05):
                witness = retune_fixed_input.proposal(self, q, coordinate, tilt, u+offset)
                trial = self.bound(lo, hi, a, b, witness)
                if trial < value:
                    value, best = trial, witness
        return best


def run(*args):
    original = budget_dense.short_dense
    budget_dense.short_dense = SimpleNamespace(Checker=Checker, mixed_dense=short_dense.mixed_dense)
    try:
        return budget_dense.run(*args)
    finally:
        budget_dense.short_dense = original


if __name__ == '__main__':
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--output', type=Path, required=True)
    p.add_argument('--seed', type=Path, required=True)
    p.add_argument('--nodes', type=int, default=200)
    p.add_argument('--seconds', type=float, default=300)
    p.add_argument('--verify', action='store_true')
    a = p.parse_args()
    run(a.output.resolve(), 16, 64, a.seed.resolve(), a.nodes, a.seconds, 42, a.verify)
