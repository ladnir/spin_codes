"""Preserve syndrome-density information on a Bernoulli zero-state exit."""
import argparse
from pathlib import Path

from flint import arb
import density_probe

model = density_probe.model


def activated(matrix, engine, theta, lam):
    z = (-lam).exp()
    g0 = 1 - theta + theta * z
    rho = min(arb(1), model.up(abs(1 - 2 * theta * z / g0)))
    cap = model.up(g0**128 * (1 + sum((count * rho**w for w, count in engine.b_spectrum.items()), arb(0)))
                   / (engine.m + 1))
    result = list(matrix)
    result[1] = arb(0)
    for k, w in enumerate(engine.levels):
        result[k + 2] = model.up(cap * engine.spectrum[w])
    return tuple(result)


class Checker(density_probe.Checker):
    def _fixed_moment(self, index, r):
        old = super()._fixed_moment(index, r)
        theta, lam = model.number(r), (arb(index) / 40).exp()
        matrix = self.engine.bernoulli(theta, lam)
        improved = activated(matrix, self.engine, theta, lam)
        moment = model.independent.terminal(improved, self.engine.n, 1 << self.power)
        return min(old, model.up(self.cutoff * lam + moment.log()))


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
