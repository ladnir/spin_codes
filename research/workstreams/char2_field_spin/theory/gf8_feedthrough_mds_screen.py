"""Width-eight positive screens for full-Cauchy packets with state4/8/16.

Start with a Cauchy matrix and normalize ALL bottom-state rows so the
state-to-state block is identity. The normalized matrix is not itself
superregular. Selected output rows still have maximal numerical rank,
as do selected output rows jointly with ALL state-return equations.
These are different designs from the cubic mu4 feedthrough.
"""
from __future__ import annotations

import argparse
from functools import partial, lru_cache
import json
import math

from gf8_feedthrough_screen import live_kernel_factor, search as base_search
from state_rank_cap_transfer import rank_cap_weights
from partial_live_rank import mds_return_bound
from diagonal_refresh_rank import rank_potential_bound


@lru_cache(maxsize=8192, typed=True)
def partial_return(weight, state, cap, tau, zeros=0):
    return mds_return_bound(weight, state, cap+weight, tau, selected_zeros=zeros)


@lru_cache(maxsize=512, typed=True)
def two_state_root(tau, roots, state):
    return max(rank_potential_bound(r, roots, tau, 256, state) for r in range(1, state+1))


def weighted_matrices(y, tau, rank_caps=True, live_kernel=True, state=4, partial_live=True):
    if state not in (4, 8, 16) or not 0 < y <= 1 or not 1 <= tau <= 255:
        raise ValueError('full-MDS screen requires state4/8/16 and valid marks')
    zero = y*0+tau*0
    states = state+1 if rank_caps else 2
    cap_rows = tuple((cap, cap) for cap in range(1, state+1)) if rank_caps else ((1, state),)
    result = [[[zero for _ in range(states)] for _ in range(states)] for _ in range(9)]
    result[0][0][0] = y**0
    for w in range(9):
        supports = math.comb(8, w)
        for z in range(9):
            factor = supports*math.comb(8, z)*y**(8-z)*(1-y)**z
            # Eliminate min(z,w) input pivots, leaving z-w Cauchy state roots.
            pivots = min(z, w)
            roots = max(0, z-w)
            if roots < state:
                if rank_caps:
                    for row, cap in cap_rows:
                        for h, mass in enumerate(rank_cap_weights(cap, roots, tau, 256, state), 1):
                            result[w][row][min(state, h+w-pivots)] += factor*mass/tau**pivots
                else:
                    result[w][1][1] += factor*two_state_root(tau, roots, state)/tau**pivots
            if w and z < w:
                result[w][0][min(state, w-z) if rank_caps else 1] += factor/tau**z
            if w > state+z:
                charge = live_kernel_factor(tau, state+z, w, w) if live_kernel else tau**(-state-z)
                result[w][0][0] += factor*charge
        if not w:
            continue
        for row, cap in cap_rows:
            pure = joint = zero
            for k in range(min(w-1, cap-1)+1):
                flat = math.comb(state, k)/tau**k
                charge = live_kernel_factor(tau, state, cap-k+w, state-k+w) if live_kernel else tau**-state
                pure += supports*flat*y**(9-w)*charge
                # A return with z selected outputs needs h+w>state+z.
                # Since h<=state-k, only z<=w-k-1 survives this flat.
                for z in range(w-k):
                    charge = (live_kernel_factor(tau, state+z, cap-k+w, state-k+w)
                              if live_kernel else tau**(-state-z))
                    joint += supports*flat*math.comb(8, z)*y**(8-z)*(1-y)**z*charge
            result[w][row][0] = min(pure, joint)
            if partial_live:
                pure_new = supports*y**(9-w)*partial_return(w, state, cap, tau)
                joint_new = supports*sum(math.comb(8, z)*y**(8-z)*(1-y)**z
                                        *partial_return(w, state, cap, tau, z)
                                        for z in range(w))
                result[w][row][0] = min(result[w][row][0], pure_new, joint_new)
    return result


def search(value, sparse=False, dimension=96, k_log=20, iterations=65, start=None, live_kernel=True,
           source_tilt=False, state=4, rank_caps=True, partial_live=True):
    result = base_search(value, sparse, dimension, k_log, iterations, start, rank_caps, live_kernel,
                         partial(weighted_matrices, state=state, partial_live=partial_live), source_tilt)
    result.update(construction=f'gf256-full{8+state}-cauchy-w8-s{state}-normalized-state-one-regional-route',
                  mu=None, eta=None, birth_minima=list(range(8, 0, -1)),
                  return_minima=list(range(8, 0, -1)), state=state,
                  partial_live_return=partial_live, live_kernel_return=live_kernel)
    return result


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--dimension', type=int, default=96)
    parser.add_argument('--density', type=float, nargs='*', default=[.03, .10, .14])
    parser.add_argument('--sparse', type=int, nargs='*', default=[64, 128])
    parser.add_argument('--iterations', type=int, default=65)
    parser.add_argument('--k-log', type=int, default=20)
    parser.add_argument('--no-live-kernel', action='store_true')
    parser.add_argument('--no-partial-live', action='store_true')
    parser.add_argument('--min-source-weight', action='store_true')
    parser.add_argument('--state', type=int, choices=(4, 8, 16), default=4)
    parser.add_argument('--two-state', action='store_true')
    args = parser.parse_args()
    for sparse, values in ((True, args.sparse), (False, args.density)):
        previous = None
        for value in values:
            result = search(value, sparse, args.dimension, args.k_log, args.iterations, previous,
                            not args.no_live_kernel, args.min_source_weight, args.state,
                            not args.two_state, not args.no_partial_live)
            previous = result['log_parameters']
            print(json.dumps(result), flush=True)
