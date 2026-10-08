"""Floating diagnostic for randomized-outer SPIN, not a certificate.

All input amplitudes are averaged before the fractional moment.  The two
states mean zero/nonzero, not uniform nonzero.  Actual regional routing is
conditioned on its exact occupancy; a Bernoulli marker bounds its coefficient.
"""

from functools import lru_cache
from math import comb, expm1, inf, log, log2
import argparse
import json

from hd_randomized_outer_proof import prescribed_mds_count

NEG = -inf


def add(a, b):
    if a < b:
        a, b = b, a
    return a if b == NEG else a + log2(1 + 2**(b-a))


def total(values):
    result = NEG
    for value in values:
        result = add(result, value)
    return result


def subtract(a, b):
    """log2(2**a-2**b); used only for the known all-zero input mass."""
    if b >= a:
        return NEG
    return a + log2(-expm1((b-a)*log(2)))


def multiply(a, b):
    return tuple(tuple(add(a[i][0]+b[0][j], a[i][1]+b[1][j])
                       for j in range(2)) for i in range(2))


def power(a, n):
    result = ((0.0, NEG), (NEG, 0.0))
    while n:
        if n & 1:
            result = multiply(result, a)
        a = multiply(a, a)
        n >>= 1
    return result


def mds_moment(a, log_z, q=65536):
    if a == 0:
        return 0.0
    terms = [0.0]
    for w in range(9-a, 9):
        count = comb(8, w)*prescribed_mds_count(w, 8, a, q)
        if count:
            terms.append(log2(count) + w*log_z)
    return total(terms) - a*log2(q)


@lru_cache(None)
def local_operators(log_z, mixer='rank_one', q=65536):
    """Log2 upper operators, averaged over a uniform j-slot set.

    An affine coset's q-ary Hamming tilt is at most the zero coset's tilt.
    This follows from the nonnegative Fourier transform of z**[x != 0].
    Untouched tiles emit A Q and contain at least 64-8u-3 nonzeros.
    """
    if not log_z < 0 or mixer not in ('rank_one', 'uniform', 'diagonal'):
        raise ValueError('Invalid tilt or state mixer')
    log_q = log2(q)
    tile = [log2(comb(8, a))+mds_moment(a, log_z, q) for a in range(9)]
    dp = {(0, 0): 0.0}
    for unused in range(8):
        nxt = {}
        for (j, u), value in dp.items():
            for a in range(9):
                key = j+a, u+(a != 0)
                nxt[key] = add(nxt.get(key, NEG), value+tile[a])
        dp = nxt
    operators = []
    for j in range(65):
        normalization = log2(comb(64, j))
        zero_moment = total(value for (jj, u), value in dp.items() if jj == j)
        zero_moment -= normalization
        live_moment = total(value + max(0, 64-8*u-3)*log_z
                            for (jj, u), value in dp.items() if jj == j)
        live_moment -= normalization
        rank = min(j, 4)
        zero_zero = min(zero_moment, -rank*log_q)
        zero_live = subtract(zero_moment, -j*log_q)
        live_zero = NEG
        if j:
            if mixer == 'rank_one':
                weighted = 3*log_q-log2(q**4-1)+live_moment
                # M Q lies in im(C_S); then a uniform rank-r syndrome
                # hits it with probability q^-r. The nonuniform part is
                # upper bounded by its full mixture mass 1/q.
                unweighted = -rank*log_q + log2(
                    1/q + (1-1/q)*(q**rank-1)/(q**4-1))
                live_zero = min(weighted, unweighted)
            elif mixer == 'uniform':
                live_zero = live_moment-log2(q**4-1)
            else:
                # Only a proposed D2 H D1 mixer, NOT the native baseline.
                live_zero = live_moment-4*log2(q-1)
        operators.append(((zero_zero, zero_live), (live_zero, live_moment)))
    return tuple(operators)


@lru_cache(None)
def pair_coefficients(log_z, alpha, mixer='rank_one', q=65536):
    """Compose chronological steps BEFORE fractional power.

    Each j,k entry already averages within-packet slot placements. Jensen's
    inequality upper bounds their alpha moment. No independence of states is
    asserted. Coefficients retain both packets' unpowered routing weights.
    """
    local = local_operators(log_z, mixer, q)
    coefficients = [[[NEG]*2 for unused in range(2)] for degree in range(129)]
    chooses = [log2(comb(64, j)) for j in range(65)]
    for j in range(65):
        for k in range(65):
            product = multiply(local[j], local[k])
            weight = chooses[j]+chooses[k]
            for a in range(2):
                for b in range(2):
                    coefficients[j+k][a][b] = add(
                        coefficients[j+k][a][b], weight+alpha*product[a][b])
    return tuple(tuple(tuple(row) for row in coefficient) for coefficient in coefficients)


def route_bound(active_rows, log_z, alpha, marker, mixer='rank_one',
                rows=8192, outer_n=256, outer_k=128, limit=838860, q=65536,
                steps=1):
    if not 0 < active_rows < rows or not 0 < marker < 1 or not 0 < alpha <= 1:
        raise ValueError('Need interior occupancy/marker and 0<alpha<=1')
    if steps == 1:
        local = local_operators(log_z, mixer, q)
        weights = [log2(comb(64, j))+j*log2(marker)+(64-j)*log2(1-marker)
                   for j in range(65)]
        transfer = tuple(tuple(total(weights[j]+alpha*local[j][a][b]
                                     for j in range(65)) for b in range(2))
                         for a in range(2))
    elif steps == 2:
        coefficients = pair_coefficients(log_z, alpha, mixer, q)
        weights = [j*log2(marker)+(128-j)*log2(1-marker) for j in range(129)]
        transfer = tuple(tuple(total(weights[j]+coefficients[j][a][b]
                                     for j in range(129)) for b in range(2))
                         for a in range(2))
    else:
        raise ValueError('Only one/two-step screens are implemented')
    run = power(transfer, rows*outer_n//(64*steps))
    tail = add(run[0][0], run[0][1])
    choose = log2(comb(rows, active_rows))
    normalizer = choose+active_rows*log2(marker)+(rows-active_rows)*log2(1-marker)
    log_beta = outer_n*log2(q)-(outer_n-outer_k)*log2(q-1)
    # Projective witness division is inside alpha; active-row union is not.
    return (choose+alpha*(active_rows*log_beta-log2(q-1)-limit*log_z)
            +tail-outer_n*normalizer)


def screen(active_rows, mixer='rank_one', steps=1):
    best = (inf, None)
    for log_z in (-1/4096, -1/2048, -1/1024, -1/512, -1/256,
                  -1/128, -1/64, -.03125, -.0625, -.125, -.25, -.5, -1, -2, -4,
                  -6, -8, -10, -12, -14, -16):
        for alpha in (1/64, 1/32, 1/16, 1/8, 1/4, 1/2, 3/4, 1):
            for factor in (.5, 1, 2, 4, 8):
                marker = active_rows/8192*factor
                if not 0 < marker < 1:
                    continue
                value = route_bound(active_rows, log_z, alpha, marker, mixer,
                                    steps=steps)
                if value < best[0]:
                    best = value, dict(log_z=log_z, alpha=alpha, marker=marker)
    return dict(active_rows=active_rows, mixer=mixer, steps=steps,
                log2_upper_diagnostic=best[0],
                witness=best[1], directed=False, whole_code_certificate=False)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--rows', nargs='+', type=int,
                        default=[1, 16, 128, 512, 1024, 2048, 4096, 6144, 7547])
    parser.add_argument('--mixer', choices=['rank_one', 'uniform', 'diagonal'],
                        default='rank_one')
    parser.add_argument('--steps', type=int, choices=[1, 2], default=1)
    args = parser.parse_args()
    for rows in args.rows:
        print(json.dumps(screen(rows, args.mixer, args.steps)), flush=True)
