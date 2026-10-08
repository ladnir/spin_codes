"""Positive rank-cap packet cover for the symmetric eight-state Cauchy inner.

State 0 means actual zero incoming witness state. States 1..8 mean actual
nonzero witness state whose form rank is at most that cap. Domination is
on nonnegative continuations nondecreasing on states 1..8, not entrywise
against an exact rank-transition matrix. No whole-code certificate.
"""
from __future__ import annotations

from fractions import Fraction
from functools import lru_cache
from math import comb


def _number(value):
    return Fraction(value) if isinstance(value, int) else value


@lru_cache(maxsize=4096, typed=True)
def moment_tail_bound(rank, roots, threshold, tau, q=256, state=8):
    """Sharp scalar moment-LP bound for E[tau^(H-rank); H>=threshold].

    Only H<=min(rank,state-roots) and
    E[q^H-1]<=(q^rank-1)/(q-1)^roots are used. H=0 contributes zero.
    The actual diagonal law need not attain this relaxed LP optimum.
    """
    if not (1 <= rank <= state and 0 <= roots <= state and
            1 <= threshold <= state and 1 <= tau <= q-1):
        raise ValueError("invalid rank-moment parameters")
    tau = _number(tau)
    zero = tau * 0
    maximum = min(rank, state-roots)
    if threshold > maximum:
        return zero
    if roots == 0:
        return tau**0 if rank >= threshold else zero
    exact = not isinstance(tau, float)
    mean = (Fraction(q**rank-1, (q-1)**roots) if exact
            else (q**rank-1)/(q-1)**roots)
    maximum_moment = Fraction(q**maximum-1) if exact else float(q**maximum-1)
    mean = min(mean, maximum_moment)
    xs = tuple(q**h-1 for h in range(maximum+1))
    ys = tuple(tau**(h-rank) if h >= threshold else zero
               for h in range(maximum+1))
    # One inequality and normalization: an optimum has <=2 support points.
    # A slack optimum is a point with x<=mean; a tight one brackets mean.
    best = max(y for x, y in zip(xs, ys) if x <= mean)
    for left in range(maximum+1):
        if xs[left] > mean:
            break
        for right in range(left+1, maximum+1):
            if xs[right] < mean:
                continue
            ratio = (mean-xs[left])/(xs[right]-xs[left])
            value = ys[left] + ratio*(ys[right]-ys[left])
            best = max(best, value)
    return best


@lru_cache(maxsize=4096, typed=True)
def rank_tail_envelope(cap, roots, tau, q=256, state=8):
    """One simultaneous nonnegative upper-tail envelope for all ranks<=cap.

    W[j-1] bounds the tilted tail at H>=j. Differences W[h-1]-W[h]
    form a positive measure dominating every admissible actual law on
    nonnegative, nondecreasing continuations. This measure is an upper
    envelope, not a claimed distribution of H.
    """
    if not 1 <= cap <= state:
        raise ValueError("rank cap must be positive")
    tau = _number(tau)
    tails = [max(moment_tail_bound(r, roots, j, tau, q, state)
                 for r in range(1, cap+1)) for j in range(1, state+1)]
    # Exact tails already decrease. This only widens any floating inversion.
    for index in range(state-2, -1, -1):
        tails[index] = max(tails[index], tails[index+1])
    return tuple(tails)


def rank_cap_weights(cap, roots, tau, q=256, state=8):
    tails = rank_tail_envelope(cap, roots, tau, q, state)
    return tuple(tails[i]-(tails[i+1] if i+1 < state else 0)
                 for i in range(state))


def return_flat_factor(weight, tau, cap, state=8):
    tau = _number(tau)
    return sum(comb(state, k)/tau**k for k in range(min(weight-1, cap-1)+1))


def return_minimum_weight(weight):
    """Symmetric-map return weight from its two superregular half-blocks.

    One occupied half gives >=9 outputs. With both halves occupied, each
    outside-input half gives >=max(0,9-weight) outputs. A return cannot have
    zero output because CA=0. This improves output weight, not rank charge.
    """
    if not 1 <= weight <= 16:
        raise ValueError("a nonempty return needs input weight1..16")
    if weight == 1:
        return 9
    if weight <= 8:
        return min(9, 2*(9-weight))
    return 1


def local_matrix(weight, y, tau, state=8, width=16, q=256):
    """9x9 positive packet matrix, ordered actual-zero then rank caps1..8.

    Valid for a fixed common source quotient and fixed active input support,
    with fresh independent nonzero input/state diagonals. State-cap input
    forms may be dependent on each other and on incoming state forms.
    Selected-zero nonreturn and pure-return covers are added separately.
    """
    if (state, width, q) != (8, 16, 256):
        raise ValueError("this helper is scoped to the checked state8,width16 geometry")
    if not (0 <= weight <= width and 0 < y <= 1 and 1 <= tau <= q-1):
        raise ValueError("invalid packet mark or weight")
    y, tau = _number(y), _number(tau)
    zero = y*0 + tau*0
    out = [[zero for _ in range(state+1)] for _ in range(state+1)]
    if weight == 0:
        out[0][0] = y**0
    else:
        out[0][min(state, weight)] = y**weight
        if weight > state:
            out[0][0] = y**weight/tau**state
    for cap in range(1, state+1):
        if weight:
            minimum = return_minimum_weight(weight)
            out[cap][0] = return_flat_factor(weight, tau, cap, state)*y**minimum/tau**state
        for i in range(weight+1):
            for z in range(min(state-1, width-weight)+1):
                factor = (comb(weight, i)*comb(width-weight, z)
                          *y**(width-i-z)*(1-y)**(i+z)/tau**i)
                if not factor:
                    continue
                for h, mass in enumerate(rank_cap_weights(cap, z, tau, q, state), 1):
                    following_cap = min(state, h+weight-i)
                    out[cap][following_cap] += factor*mass
    return tuple(tuple(row) for row in out)


def monotone_continuation(matrix, values):
    """Evaluate only after checking the continuation cone in the contract."""
    if len(values) != len(matrix) or any(value < 0 for value in values):
        raise ValueError("continuation must be nonnegative with matching dimension")
    if any(values[j] > values[j+1] for j in range(1, len(values)-1)):
        raise ValueError("nonzero-cap continuation must be nondecreasing")
    return tuple(sum(a*b for a, b in zip(row, values)) for row in matrix)
