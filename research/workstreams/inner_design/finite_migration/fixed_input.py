"""A fixed input tilt avoids paying for every possible band composition."""
import argparse
from fractions import Fraction as F
from functools import lru_cache
import math
from pathlib import Path

from flint import arb, ctx
import bernoulli_activation
import density_comparison
import density_probe
import ladder_dense_tilt as scalar
import poisson_density_factor as old_density

model = density_probe.model


class Checker(bernoulli_activation.Checker):
    def __init__(self, exponent):
        super().__init__(exponent)
        self.comparison = lru_cache(maxsize=2048)(self._comparison)

    def _comparison(self, lo, hi, a, b, precision):
        assert precision == ctx.prec
        nlo, nhi = self.density(a, b)
        old = arb(old_density.density_factor(self.rows))
        return old if nhi == 1 else min(old, density_comparison.factor(
            self.rows, lo, hi, nlo, nhi, self.ps))

    def fixed_input_bound(self, lo, hi, a, b, witness):
        assert 1 <= lo <= hi <= self.rows and 0 <= a <= b <= 1
        nlo, nhi = self.density(a, b)
        theta_lo, theta_hi = F(lo, self.rows) * nlo, F(hi, self.rows) * nhi
        theta = (theta_lo + theta_hi) / 2
        r = model.base.decode(witness['fixed_r'])
        eta = model.base.decode(witness['eta'])
        assert 0 < r < 1 and 0 < theta < 1 and abs(eta) <= 10000
        # This x is fixed across the entire box, not chosen per composition.
        x = theta * (1 - r) / (r * (1 - theta))
        costs = scalar.tilted_costs(self.bands, self.ps, self.caps, model.number(x).log())
        log_s = sum(((g - model.number(eta * p)).exp()
                     for g, p in zip(costs, self.ps + (F(1),))), arb(0)).log().upper()
        scalar_upper = scalar.scalar_bound(self.rows, lo, hi, nlo, nhi, log_s, eta, model.number(x))
        lower = theta_lo / (x + (1 - x) * theta_lo)
        upper = theta_hi / (x + (1 - x) * theta_hi)
        if upper == 1:
            return None
        probability = model.number(lower).union(model.number(upper))
        lam = (arb(witness['tilt']) / 40).exp()
        matrix = self.engine.bernoulli(probability, lam)
        activated = bernoulli_activation.activated(matrix, self.engine, probability, lam)
        n, epochs = self.engine.n, 1 << self.power
        moment = min(model.up(model.independent.terminal(matrix, n, epochs)),
                     model.up(model.independent.terminal(activated, n, epochs)))
        # This exact sum counts active positions, including all Q in the box.
        locations = sum(math.comb(self.rows, q) for q in range(lo, hi + 1))
        comparison = self.comparison(lo, hi, a, b, ctx.prec)
        return model.up(arb(locations).log() + 256 * comparison.log()
                        + scalar_upper + self.cutoff * lam + moment.log())

    def bound(self, lo, hi, a, b, witness):
        old = super().bound(lo, hi, a, b, witness)
        fixed = self.fixed_input_bound(lo, hi, a, b, witness)
        return old if fixed is None else min(old, fixed)


def run(output, verify=False):
    original = density_probe.Checker
    density_probe.Checker = Checker
    try:
        return density_probe.run(output, verify)
    finally:
        density_probe.Checker = original


if __name__ == '__main__':
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--output', type=Path, required=True)
    p.add_argument('--verify', action='store_true')
    a = p.parse_args()
    run(a.output.resolve(), a.verify)
