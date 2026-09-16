"""Bounded dense-cover search with composition-sensitive routing bounds."""
import argparse
from pathlib import Path
from types import SimpleNamespace

from flint import arb, ctx
import dense_ladder
import fixed_input
import mixed_dense
import poisson_density_factor as old_density


class Checker(fixed_input.Checker):
    def bound(self, lo, hi, a, b, witness):
        # Retain an already-strong old bound without computing the new
        # O(L) density envelope. This is proof-work avoidance, not a shortcut.
        old = mixed_dense.Checker.bound(self, lo, hi, a, b, witness)
        if old < -90 * arb(2).log():
            return old
        comparison = self.comparison(lo, hi, a, b, ctx.prec)
        improved = (old + 256 * (comparison.log() - arb(old_density.density_factor(self.rows)).log())).upper()
        if improved < -90 * arb(2).log():
            return improved
        fixed = self.fixed_input_bound(lo, hi, a, b, witness)
        return improved if fixed is None else min(improved, fixed)


def run(*args):
    original = dense_ladder.ladder
    dense_ladder.ladder = SimpleNamespace(Checker=Checker, candidate=mixed_dense.ladder.candidate)
    try:
        return dense_ladder.run(*args)
    finally:
        dense_ladder.ladder = original


if __name__ == '__main__':
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--output', type=Path, required=True)
    p.add_argument('--m', type=int, choices=(16, 18), default=18)
    p.add_argument('--minimum', type=int, default=512)
    p.add_argument('--seed', type=Path)
    p.add_argument('--nodes', type=int, default=100)
    p.add_argument('--seconds', type=float, default=180)
    p.add_argument('--verify', action='store_true')
    a = p.parse_args()
    run(a.output.resolve(), a.m, a.minimum, a.seed.resolve() if a.seed else None,
        a.nodes, a.seconds, a.verify)
