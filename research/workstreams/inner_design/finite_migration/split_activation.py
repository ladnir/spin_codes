"""All-one-split bounds using bounded-density IMT activation."""
import argparse
from fractions import Fraction as F
from pathlib import Path
from types import SimpleNamespace

import activation_density
import split_imt


def run(*args):
    original = split_imt.ladder
    split_imt.ladder = SimpleNamespace(Engine=activation_density.Engine,
                                      candidate=activation_density.ladder.candidate)
    try:
        return split_imt.run(*args)
    finally:
        split_imt.ladder = original


if __name__ == '__main__':
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--output', type=Path, required=True)
    p.add_argument('--m', type=int, choices=(16, 18), default=16)
    p.add_argument('--first', type=int, default=218)
    p.add_argument('--last', type=int, default=218)
    p.add_argument('--witness', action='append')
    p.add_argument('--verify', action='store_true')
    a = p.parse_args()
    if not a.verify and not a.witness:
        p.error('Provide witness proposals')
    jobs = []
    for text in a.witness or []:
        tilt, q, h = text.split(':')
        jobs.append(dict(tilt=str(F(tilt) / 10), anchor=[int(q), int(h)]))
    run(a.output.resolve(), a.m, a.first, a.last, jobs, a.verify)
