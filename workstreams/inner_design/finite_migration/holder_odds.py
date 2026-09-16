"""Propose a common odds shift; verify using the stored rational probabilities."""
import argparse
from fractions import Fraction as F
from pathlib import Path

import holder_split as holder


def run(output, exponent, q, factors, tilt, verify=False):
    original = holder.base.choose
    proposed = iter(factors)

    def choose(engine, region, ordinary):
        factor = next(proposed)
        assert factor > 0
        return [factor * p / (1 - p + factor * p)
                for p in original(engine, region, ordinary)]

    holder.base.choose = choose
    try:
        jobs = [dict(tilt=str(tilt), anchor=[q, 0]) for _ in factors]
        return holder.run(output, exponent, q, q, jobs, verify)
    finally:
        holder.base.choose = original


if __name__ == '__main__':
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--output', type=Path, required=True)
    p.add_argument('--m', type=int, choices=(16, 18), default=16)
    p.add_argument('--q', type=int, default=218)
    p.add_argument('--tilt', type=F, default=F(-1, 8))
    p.add_argument('--factors', nargs='+', type=F, default=list(map(F, ['4/5', '9/10', '19/20', '21/20', '11/10', '6/5'])))
    p.add_argument('--verify', action='store_true')
    a = p.parse_args()
    run(a.output.resolve(), a.m, a.q, a.factors, a.tilt, a.verify)
