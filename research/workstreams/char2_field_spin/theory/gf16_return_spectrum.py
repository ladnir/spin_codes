"""Exact sparse return spectra and floating whole-code screens using them.

For w<=s, the s-by-w MDS feedback restriction has trivial kernel. On a
nonzero-to-zero transition R=CX, the emitted word is (I+AC)X. Uniform GLs
refresh gives each nonzero R mass 1/(16**s-1). Enumerate supports and
projective input values (first active value fixed to one) exactly.
The histograms are integer counts; whole-code searches are not certificates.
"""

import argparse
from functools import lru_cache
import itertools
import json
import math
from unittest.mock import patch

import numpy as np

from batched_transfer import multiply
import gf16_sizes
from state_screen import local as generic_local


@lru_cache(maxsize=None)
def return_matrix(state):
    if not 1 <= state <= 8:
        raise ValueError('supported monomial state is one through eight')
    a = np.empty((16, state), dtype=np.uint8)
    for j in range(16):
        power = 1
        for p in range(state):
            a[j, p] = power
            power = multiply(power, j, 0x13)
    table = np.array([[multiply(x, y, 0x13) for y in range(16)] for x in range(16)], dtype=np.uint8)
    b = np.eye(16, dtype=np.uint8)
    for p in range(state):
        b ^= table[a[:, p, None], a[None, :, p]]
    return b, table


@lru_cache(maxsize=None)
def histogram(state, weight):
    if not 1 <= weight <= min(4, state):
        raise ValueError('enumeration supports only weights one through min(4,s)')
    matrix, table = return_matrix(state)
    # product[coordinate, coefficient-1, output-coordinate].
    product = table[np.arange(1, 16)[None, :, None], matrix.T[:, None, :]]
    choices = np.array(list(itertools.product(range(15), repeat=weight-1)), dtype=np.intp)
    cases = 15**(weight-1)
    result = np.zeros(17, dtype=np.int64)
    for support in itertools.combinations(range(16), weight):
        values = np.broadcast_to(matrix[:, support[0]], (cases, 16)).copy()
        for p in range(1, weight):
            values ^= product[support[p], choices[:, p-1], :]
        result += np.bincount(np.count_nonzero(values, axis=1), minlength=17)
    counts = tuple(int(x) for x in result)
    if sum(counts) != math.comb(16, weight)*cases or counts[0]:
        raise ArithmeticError('return matrix enumeration violates invertibility or count')
    return counts


def local(z, state=4, maximum=4):
    result = gf16_sizes.local_model(z, state)
    denominator = 16**state-1
    for w in range(2, min(state, maximum)+1):
        counts = histogram(state, w)
        moment = sum(count*z**v for v, count in enumerate(counts))/sum(counts)/denominator
        # The existing nonreturn envelope subtracts the valid lower bound z^16/Q.
        # Replace both return entries by their exact joint moment here.
        result[w, 1, 1] += z**16/denominator-moment
        result[w, 1, 0] = moment
    return result


def search(k_log, active, state, maximum=4, optimize_source=True):
    # Capture original provider before patching: local() must not recurse.
    original = gf16_sizes.local_model
    def provider(z, s=4):
        result = original(z, s)
        denominator = 16**s-1
        for w in range(2, min(s, maximum)+1):
            counts = histogram(s, w)
            moment = sum(count*z**v for v, count in enumerate(counts))/sum(counts)/denominator
            result[w, 1, 1] += z**16/denominator-moment
            result[w, 1, 0] = moment
        return result
    with patch.object(gf16_sizes, 'local_model', provider):
        answer = (gf16_sizes.search_source(k_log, active, state) if optimize_source and active > 4
                  else gf16_sizes.search(k_log, active, state=state,
                                        source_rho=.871 if active > 4 else None))
    return dict(answer, exact_return_through=maximum)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--state', type=int, default=4)
    parser.add_argument('--maximum', type=int, choices=[1, 2, 3, 4], default=4)
    parser.add_argument('--k-log', type=int, choices=[16, 18, 20], default=20)
    parser.add_argument('--rows', nargs='+', type=int, default=[1, 4, 8, 16, 32, 64, 128, 256])
    parser.add_argument('--histograms', action='store_true')
    args = parser.parse_args()
    for w in range(1, min(args.state, args.maximum)+1):
        counts = histogram(args.state, w)
        if args.histograms:
            print(json.dumps(dict(state=args.state, input_weight=w, counts=counts, exact=True)), flush=True)
    if not args.histograms:
        for active in args.rows:
            print(json.dumps(search(args.k_log, active, args.state, args.maximum)), flush=True)
