"""Search both complete IMT moment bounds, retaining the same outward verifier."""
import argparse
from pathlib import Path
from types import SimpleNamespace

import ladder
import dense_ladder
import ladder_dense_fixed_reference as occupation_search


class Checker(ladder.Checker):
    def witness(self, lo, hi, a, b):
        # The frozen inherited search uses only the direct Fourier transfer
        # in discovery. The original fixed-reference search proposes from
        # the fixed-weight occupation transfer instead. Both proposals are
        # checked by the unchanged complete-moment bound; no matrix entries
        # are mixed. Proposal failures do not supply certificates.
        direct = super().witness(lo, hi, a, b)
        occupation = occupation_search.Checker.witness(self, lo, hi, a, b)
        return min((direct, occupation),
                   key=lambda w: self.bound(lo, hi, a, b, w))


def run(*args):
    # Inject a search backend into this driver only. Do not replace a class
    # in the frozen producer module or change any historical source bytes.
    original = dense_ladder.ladder
    dense_ladder.ladder = SimpleNamespace(Checker=Checker, candidate=ladder.candidate)
    try:
        return dense_ladder.run(*args)
    finally:
        dense_ladder.ladder = original


if __name__ == '__main__':
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--output', type=Path, required=True)
    p.add_argument('--m', type=int, choices=ladder.EXPONENTS, default=22)
    p.add_argument('--minimum', type=int, default=512)
    p.add_argument('--seed', type=Path)
    p.add_argument('--nodes', type=int, default=300)
    p.add_argument('--seconds', type=float, default=240)
    p.add_argument('--verify', action='store_true')
    a = p.parse_args()
    if not 2 <= a.minimum <= 1 << (a.m - 7) or a.nodes < 0 or a.seconds <= 0:
        p.error('Invalid range or search budget')
    run(a.output.resolve(), a.m, a.minimum, a.seed.resolve() if a.seed else None,
        a.nodes, a.seconds, a.verify)
