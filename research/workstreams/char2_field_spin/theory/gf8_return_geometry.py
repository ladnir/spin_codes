"""Exact weight-one/two return geometry for restricted symmetric GF256 maps.

On a nonzero-to-zero state transition, R=CX and Y=(I+2AC)X. Enumerating
projectively normalized weight-one/two inputs gives a deterministic output
weight bound, without consuming or reusing any random diagonal. This is a
packet fact, not a whole-code distance certificate.
"""
from __future__ import annotations

import argparse
from functools import lru_cache
from itertools import combinations
import json

import numpy as np

from state_capacity_reference import matrices, multiplication_table, mul


@lru_cache(maxsize=1)
def table():
    return np.asarray([list(row) for row in multiplication_table()], dtype=np.uint8)


@lru_cache(maxsize=32)
def return_matrix(state, phase):
    if state not in (4, 6, 8) or not 0 <= phase < 4:
        raise ValueError('checked restrictions are states4/6/8 and phases0..3')
    a, c = matrices(phase, True)
    return tuple(tuple(int(i == j) ^ mul(2, _sum(
        mul(a[i][h], c[h][j]) for h in range(state)))
        for j in range(16)) for i in range(16))


def _sum(values):
    result = 0
    for value in values:
        result ^= value
    return result


def enumerate_phase(state, phase):
    """All nonzero scalar orbits: 16 singletons and 120*255 pairs."""
    b = np.asarray(return_matrix(state, phase), dtype=np.uint8)
    products = table()
    ratios = np.arange(1, 256, dtype=np.uint8)
    singleton_weights = np.count_nonzero(b, axis=0)
    minimum_one = int(singleton_weights.min())
    one_index = int(singleton_weights.argmin())
    minimum_two, witness = 17, None
    for j, k in combinations(range(16), 2):
        output = b[:, j, None] ^ products[b[:, k, None], ratios]
        weights = np.count_nonzero(output, axis=0)
        index = int(weights.argmin())
        if int(weights[index]) < minimum_two:
            minimum_two = int(weights[index])
            witness = {'positions': [j, k], 'values': [1, index+1]}
    return dict(state=state, phase=phase, eta=2,
                minimum_return_weight=[minimum_one, minimum_two],
                witnesses=[{'positions': [one_index], 'values': [1]}, witness],
                projective_inputs_checked=16+120*255,
                exact=True, whole_code=False)


def certify_sparse_returns(state):
    phases = [enumerate_phase(state, phase) for phase in range(4)]
    return dict(state=state, eta=2,
                minimum_return_weight=[min(p['minimum_return_weight'][w]
                                           for p in phases) for w in range(2)],
                phases=phases, exact=True, whole_code=False)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--state', type=int, choices=(4, 6, 8), nargs='+', default=[4, 6, 8])
    for state in parser.parse_args().state:
        print(json.dumps(certify_sparse_returns(state)), flush=True)
