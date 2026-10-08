"""Bounded proof-side FFT8/FFT16 quiet comparison; no encoder or timings.

The float search chooses witnesses only. Directed replay evaluates the same
positive rank-only moment and the same marker bound with Arb.
"""

from fractions import Fraction as F
from functools import lru_cache
from math import comb, log2

from hd_randomized_outer_proof import Parameters


def novel_basis(index, x):
    from hd_blockwise_maps import mul, power
    from hd_blockwise_shared_maps import subspace_polynomial
    if not 0 <= index < 16:
        raise ValueError('basis index must be below 16')
    result = 1
    for level in range(4):
        if index & (1 << level):
            denominator = subspace_polynomial(level, 1 << level)
            result = mul(result, mul(subspace_polynomial(level, x), power(denominator, 254)))
    return result


@lru_cache(None)
def evaluation(offset, width):
    if width not in (4, 8, 16) or offset % width:
        raise ValueError('aligned width 4, 8, or 16 required')
    return tuple(tuple(novel_basis(j, offset ^ i) for j in range(width))
                 for i in range(width))


@lru_cache(None)
def matrices(width=16):
    """Return F,A,C,T for translation 64 and highest-four feedback."""
    from hd_blockwise_maps import matmul
    from hd_blockwise_shared_maps import inverse
    if width not in (8, 16):
        raise ValueError('probe supports widths 8 and 16 only')
    feedthrough, injection, feedback = [], [], [[] for _ in range(4)]
    for offset in range(0, 64, width):
        inv = inverse(evaluation(offset, width))
        local = matmul(evaluation(offset ^ 64, width), inv)
        feedthrough.extend(tuple([0]*offset+list(row)+[0]*(64-width-offset))
                           for row in local)
        injection.extend(row[:4] for row in evaluation(offset ^ 64, width))
        for j in range(4):
            feedback[j].extend(inv[width-4+j])
    translation = matmul(inverse(evaluation(0, 4)), evaluation(64, 4))
    return tuple(feedthrough), tuple(injection), tuple(map(tuple, feedback)), translation


@lru_cache(None)
def local_terms(mark, width, field=65536):
    """Positive tile expansion, scaled by q^a, keyed by a and min(4,d)."""
    mark = F(mark)
    result = [[F(0) for _ in range(5)] for _ in range(width+1)]
    for a in range(width+1):
        for v in range(width+1):
            d = max(a-v, 0)
            result[a][min(4, d)] += (comb(width, a)*comb(width, v)
                *mark**(width-v)*(1-mark)**v*field**d)
    return tuple(map(tuple, result))


def moments(mark, width, *, directed=False, field=65536):
    """K_j upper bounds and exact free-final T_j (logs in float mode).

    Uses exactly rho=max(max_tile min(4,d),min(2,number of positive d)).
    No 3+3 refinement is included, to keep the comparison fixed.
    """
    mark = F(mark)
    if width not in (8, 16) or not 0 < mark <= 1:
        raise ValueError('need width 8/16 and 0<mark<=1')
    if directed:
        from flint import arb
        zero, one = arb(0), arb(1)
        convert = lambda value: arb(value.numerator)/value.denominator
    else:
        zero, one, convert = 0.0, 1.0, float
    if mark == 1:
        if directed:
            return tuple(one/field**min(j, 4) for j in range(65)), (one,)*65
        return tuple(-min(j, 4)*log2(field) for j in range(65)), (0.0,)*65
    local = tuple(tuple(convert(x) for x in row) for row in local_terms(mark, width, field))
    dp = [[one, zero, zero, zero, zero]]
    for _ in range(64//width):
        nxt = [[zero for _ in range(5)] for _ in range(len(dp)+width)]
        for j, row in enumerate(dp):
            for rank, previous in enumerate(row):
                if not previous:
                    continue
                for a, terms in enumerate(local):
                    for residual, coefficient in enumerate(terms):
                        if coefficient:
                            successor = 2 if rank == residual == 1 else max(rank, residual)
                            nxt[j+a][successor] += previous*coefficient
        dp = nxt
    quiet, free = [], []
    for j, row in enumerate(dp):
        weighted = sum(row[k]/field**k for k in range(5))
        unweighted = sum(row)
        if directed:
            denominator = field**j*comb(64, j)
            bound, cap = weighted/denominator, one/field**min(j, 4)
            quiet.append(bound if bound < cap else cap)
            free.append(unweighted/denominator)
        else:
            denominator = j*log2(field)+log2(comb(64, j))
            quiet.append(min(log2(weighted)-denominator, -min(j, 4)*log2(field)))
            free.append(log2(unweighted)-denominator)
    return tuple(quiet), tuple(free)


def screen(occupancies=(1024, 2048, 4096)):
    """Same bounded mark/alpha grid for both widths, optimized marker odds."""
    import numpy as np
    from scipy.optimize import minimize_scalar
    from scipy.special import logsumexp
    p = Parameters()
    marks = tuple(F(1, 1 << shift) for shift in (0, 1, 2, 4, 6, 8, 10, 12, 14, 16))
    marks += (F(2, 3*(p.field-1)), F(133, 128*(p.field-1)))
    alphas = tuple(F(a, 32) for a in (1, 2, 3, 4, 6, 8, 12, 16, 24, 32))
    log_combinations = np.array([log2(comb(64, j)) for j in range(65)])
    indices = np.arange(65)
    ln2 = np.log(2)
    beta = 256*log2(p.field)-128*log2(p.field-1)
    result = []
    for width in (8, 16):
        best = {r: None for r in occupancies}
        for mark in marks:
            quiet, free = (np.asarray(x) for x in moments(mark, width))
            for alpha in alphas:
                a = float(alpha)
                k, t = log_combinations+a*quiet, log_combinations+a*free
                for rows in occupancies:
                    fixed = (-255*log2(comb(p.rows, rows))
                        +a*(rows*beta-log2(p.field-1)-p.limit*log2(float(mark))))
                    def objective(odds_log):
                        ordinary = logsumexp((k+indices*odds_log)*ln2)/ln2
                        final = logsumexp((t+indices*odds_log)*ln2)/ln2
                        return fixed+(p.packets-1)*ordinary+final-256*rows*odds_log
                    optimum = minimize_scalar(objective, bounds=(-16, 24), method='bounded')
                    marker = F(round((1/(1+2**(-optimum.x)))*(1 << 24)), 1 << 24)
                    value = objective(log2(float(marker/(1-marker))))
                    if best[rows] is None or value < best[rows]['screen_log2']:
                        best[rows] = dict(width=width, active_rows=rows, mark=str(mark),
                            marker=str(marker), alpha=str(alpha), screen_log2=float(value))
        result.extend(best[r] for r in occupancies)
    return tuple(result)


def directed_endpoint(width, rows, mark, marker, alpha):
    from flint import arb
    p = Parameters()
    mark, marker, alpha = map(F, (mark, marker, alpha))
    quiet, free = moments(mark, width, directed=True)
    convert = lambda value: arb(value.numerator)/value.denominator
    z, m, a = map(convert, (mark, marker, alpha))
    ordinary, final = arb(0), arb(0)
    for j in range(65):
        weight = comb(64, j)*m**j*(1-m)**(64-j)
        ordinary += weight*quiet[j]**a
        final += weight*free[j]**a
    logarb = lambda x: x.log()/arb(2).log()
    choose = logarb(arb(comb(p.rows, rows)))
    correction = 256*(choose+rows*logarb(m)+(p.rows-rows)*logarb(1-m))
    beta = 256*logarb(arb(p.field))-128*logarb(arb(p.field-1))
    value = (choose+a*(rows*beta-logarb(arb(p.field-1))-p.limit*logarb(z))
        +(p.packets-1)*logarb(ordinary)+logarb(final)-correction)
    return dict(width=width, active_rows=rows, mark=str(mark), marker=str(marker),
        alpha=str(alpha), failure_log2_enclosure=str(value),
        below_2_to_minus_40=bool(value < -40), whole_code_certificate=False)


def main():
    import json
    from flint import ctx
    ctx.prec = 160
    for witness in screen():
        print(json.dumps(witness), flush=True)
        print(json.dumps(directed_endpoint(witness['width'], witness['active_rows'],
            witness['mark'], witness['marker'], witness['alpha'])), flush=True)


if __name__ == '__main__':
    main()
