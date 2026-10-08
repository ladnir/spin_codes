"""Conservative, unconditional common-quotient packet rank operator.

This is a proof screen, not a native encoder.  State zero means that the
state equations were imposed.  States 1..s are ranks of a required-nonzero
state map.  Unknown ranks are bounded entrywise, never silently identified.
See RANK_PACKET_OPERATOR.md for the conditional theorem.
"""

from collections import defaultdict
from functools import lru_cache
from math import comb, factorial

import numpy as np


def rank_terms(f0, d0, f1, d1, output_only, empty, s=4):
    """Yield (incoming, outgoing, tau_charge, exceptional) for one support type.

    f0/d0 are fresh/dependent inputs outside J; f1/d1 are inside J.
    Each yielded entry has coefficient tau**(-charge)/(q-1)**exceptional.
    Repeated entries are ADDED.  This intentional overcount makes dominance
    entrywise, valid for arbitrary nonnegative continuation vectors.
    """
    a = f1 + d1
    w = f0 + d0 + a
    v = a + output_only
    t = w + output_only + empty
    if f0 == d0 == output_only == 0:
        if w == 0:
            yield 0, 0, 0, 0
        else:
            if w >= s + 1:
                yield 0, 0, s, 0
            lo = max(1, min(s, f1))
            hi = min(s, w)
            for k in range(lo, hi + 1):
                yield 0, k, 0, 0
    if empty >= s:
        return
    for r in range(1, s + 1):
        # Root equations reduce r to u.  Dependent zero-output equations
        # can reduce u by at most d0; fresh ones preserve the old state rank.
        for h in range(max(1, r - empty - d0), min(r, s - empty) + 1):
            largest_u = min(r, s - empty, h + d0)
            charge = f0 + d0 + r - largest_u
            if v == 0:
                yield r, h, charge, 0
                continue
            if v == 1:
                if a == 0:
                    yield r, h, charge, 0
                    continue  # invertible transvection cannot return
                if f1:
                    lo, hi = (h+1 if empty else h), min(s, h + 1)
                    for k in range(lo, hi + 1):
                        yield r, k, charge, 0
                    return_rank = h
                else:
                    yield r, h, charge, 0
                    if h < s:
                        yield r, h + 1, charge, 0
                    if h > 1 and not empty:
                        yield r, h - 1, charge, 1
                    return_rank = max(0, h - 1)
            else:
                lo = max(0, min(s, f1), h - v - d1)
                hi = min(s, h + a)
                for k in range(max(1, lo), hi + 1):
                    yield r, k, charge, 0
                return_rank = lo
            if w and v and not (v == 1 and w < t):
                # Return pivots and deterministic state rank use the SAME
                # equations/randomness.  Take max, never add these charges.
                yield r, 0, charge + max(min(a, s), return_rank), 0


@lru_cache(maxsize=None)
def dependent_updates(rank, count, s=4):
    """Entrywise rank law for fresh diagonals on arbitrary nonzero forms."""
    states = {(rank, 0): 1}
    for _ in range(count):
        next_states = defaultdict(int)
        for (r, exponent), multiplicity in states.items():
            if r == 0:
                next_states[1, exponent] += multiplicity
            else:
                next_states[r, exponent] += multiplicity
                if r < s:
                    next_states[r + 1, exponent] += multiplicity
                next_states[r - 1, exponent + 1] += multiplicity
        states = next_states
    return tuple((r, exponent, mult) for (r, exponent), mult in states.items())


def scaled_rank_terms(f0, d0, f1, d1, o, z, s=4):
    """Weighted terms for the separate eta-scaled, all-F_J-invertible law."""
    a, v = f1 + d1, f1 + d1 + o
    for r, k, charge, exceptional in rank_terms(f0, d0, f1, d1, o, z, s=s):
        if r == 0 or v <= 1:
            yield r, k, charge, exceptional, 1
    if v <= 1 or z >= s:
        return
    for r in range(1, s + 1):
        for h in range(max(1, r - z - d0), min(r, s - z) + 1):
            charge = f0 + d0 + r - min(r, s - z, h + d0)
            for middle, exceptional, multiplicity in dependent_updates(h, d1, s):
                lo, hi = max(middle, min(s, f1)), min(s, middle + f1)
                for k in range(max(1, lo), hi + 1):
                    yield r, k, charge, exceptional, multiplicity
            t = f0+d0+f1+d1+o+z
            if a and not (v == 2 and f0+d0+a < t-s+1):
                return_rank = max(min(s, f1), h - d1)
                yield r, 0, charge + max(min(a, s), return_rank), 0, 1


@lru_cache(maxsize=None)
def coefficient_table(t=16, s=4, law='unscaled'):
    """Exact integer terms indexed by (f,d,v,r,k,charge,exceptional)."""
    if law not in ('unscaled', 'eta-scaled'):
        raise ValueError('unknown construction law')
    result = defaultdict(int)
    fact = [factorial(i) for i in range(t + 1)]
    for f0 in range(t + 1):
        for d0 in range(t - f0 + 1):
            for f1 in range(t - f0 - d0 + 1):
                for d1 in range(t - f0 - d0 - f1 + 1):
                    for o in range(t - f0 - d0 - f1 - d1 + 1):
                        z = t - f0 - d0 - f1 - d1 - o
                        counts = f0, d0, f1, d1, o, z
                        multiplicity = fact[t]
                        for n in counts:
                            multiplicity //= fact[n]
                        f, d, v = f0 + f1, d0 + d1, f1 + d1 + o
                        if law == 'unscaled':
                            terms = ((*term, 1) for term in rank_terms(*counts, s=s))
                        else:
                            terms = scaled_rank_terms(*counts, s=s)
                        for r, k, charge, exceptional, weight in terms:
                            result[f, d, v, r, k, charge, exceptional] += multiplicity * weight
    return dict(result)


@lru_cache(maxsize=None)
def _arrays(t, s, law):
    table = coefficient_table(t, s, law)
    keys = np.asarray(list(table), dtype=np.int64)
    values = np.asarray(list(table.values()), dtype=float)
    f, d, v, r, k, charge, exceptional = keys.T
    index = (((f * (t + 1) + d) * (s + 1) + r) * (s + 1) + k)
    normalizer = np.zeros((t + 1, t + 1))
    for fresh in range(t + 1):
        for dep in range(t - fresh + 1):
            normalizer[fresh, dep] = comb(t, fresh) * comb(t - fresh, dep)
    return index, values, v, charge, exceptional, normalizer


def packet_matrices(tau, y, q=256, t=16, s=4, normalized=False, law='unscaled'):
    """Return K[f,d,r,k], summed over placements and marked outputs.

    If normalized=True divide by multinomial(t;f,d,t-f-d), giving the
    uniform-placement average N_fd used by outer_freshness_transfer.
    Floating evaluation is exploratory; coefficient_table is exact.
    """
    if not 1 <= tau <= q - 1 or not 0 < y <= 1:
        raise ValueError('require 1 <= tau <= q-1 and 0 < y <= 1')
    index, values, v, charge, exceptional, normalizer = _arrays(t, s, law)
    weights = values * y**v * tau**(-charge.astype(float)) / float(q - 1)**exceptional
    result = np.bincount(index, weights=weights,
                        minlength=(t + 1)**2 * (s + 1)**2)
    result = result.reshape(t + 1, t + 1, s + 1, s + 1)
    if normalized:
        np.divide(result, normalizer[:, :, None, None], out=result,
                  where=normalizer[:, :, None, None] != 0)
    return result


def dependent_action(values, q=256):
    """Sharp one-diagonal upper action on any nonnegative rank payoff."""
    result = [values[1]]
    for r in range(1, len(values)):
        ordinary = max(values[r:r+2])
        result.append(ordinary + max(0, values[r-1] - ordinary)/(q-1))
    return np.asarray(result, dtype=float)


@lru_cache(maxsize=None)
def _support_types(t):
    counts, multiplicities = [], []
    fact = [factorial(i) for i in range(t + 1)]
    for f0 in range(t + 1):
        for d0 in range(t - f0 + 1):
            for f1 in range(t - f0 - d0 + 1):
                for d1 in range(t - f0 - d0 - f1 + 1):
                    for o in range(t - f0 - d0 - f1 - d1 + 1):
                        z = t - f0 - d0 - f1 - d1 - o
                        count = f0, d0, f1, d1, o, z
                        multiplicity = fact[t]
                        for n in count:
                            multiplicity //= fact[n]
                        counts.append(count)
                        multiplicities.append(multiplicity)
    return np.asarray(counts, dtype=int).T, np.asarray(multiplicities, dtype=float)


def gaussian_binomial(n, k, q):
    if k < 0 or k > n:
        return 0
    numerator, denominator = 1, 1
    for i in range(k):
        numerator *= q**(n-i)-1
        denominator *= q**(k-i)-1
    return numerator//denominator


@lru_cache(maxsize=None)
def root_intersection_probability(q, s, r, z, h):
    """Exact dimension law for U intersect ker(A), U uniform in Gr(r,s)."""
    from fractions import Fraction
    if not 0 <= z <= s or not max(0, r-z) <= h <= min(r, s-z):
        return Fraction(0)
    count = (gaussian_binomial(s-z, h, q)*gaussian_binomial(z, r-h, q)
             *q**((s-z-h)*(r-h)))
    return Fraction(count, gaussian_binomial(s, r, q))


def packet_action(tau, y, continuation, q=256, t=16, s=4, normalized=False,
                  law='eta-scaled', fresh_orientation=False):
    """Uniform Bellman action K_fd(V), without duplicating unknown-rank mass.

    Returns shape(t+1,t+1,s+1).  This is NONLINEAR in V.  Matrix powers of
    values returned here are invalid; use a positive supersolution or compose
    the actions.  Currently requires the separate all-F_J-invertible law.
    """
    if law != 'eta-scaled':
        raise ValueError('Bellman action currently requires eta-scaled law')
    values = np.asarray(continuation, dtype=float)
    if (values.shape != (s+1,) or np.any(values < 0) or
            not 1 <= tau <= q-1 or not 0 < y <= 1):
        raise ValueError('invalid continuation or tilts')
    nonzero = values.copy()
    nonzero[0] = 0
    both = nonzero + values[0] * tau**(-np.arange(s+1, dtype=float))
    allowed = np.zeros((t+1, t+1, s+1))
    allowed_nonzero = np.zeros_like(allowed)
    for fresh in range(t+1):
        g = np.asarray([max(nonzero[max(r, min(s, fresh)):min(s, r+fresh)+1])
                        for r in range(s+1)])
        b = np.asarray([max(both[max(r, min(s, fresh)):min(s, r+fresh)+1])
                        for r in range(s+1)])
        for dep in range(t-fresh+1):
            a = fresh + dep
            allowed_nonzero[fresh, dep] = g
            if a:
                allowed[fresh, dep] = np.minimum(b, g + values[0]/tau**min(a, s))
            else:
                # No overlap means the state undergoes an invertible map,
                # so an incoming nonzero state cannot return.
                allowed[fresh, dep] = nonzero
            g, b = dependent_action(g, q), dependent_action(b, q)
    counts, multiplicities = _support_types(t)
    f0, d0, f1, d1, o, z = counts
    f, d, v = f0+f1, d0+d1, f1+d1+o
    w = f+d
    result_by_type = np.zeros((s+1, len(f0)))
    # The zero-state row uses exact I=J and the MDS return threshold.
    for i in np.nonzero((f0 == 0) & (d0 == 0) & (o == 0))[0]:
        if w[i] == 0:
            result_by_type[0, i] = values[0]
            continue
        low, high = max(1, min(s, f1[i])), min(s, w[i])
        active = max(nonzero[low:high+1])
        if w[i] <= s:
            result_by_type[0, i] = active
        else:
            formal_low = min(s, f1[i])
            result_by_type[0, i] = min(
                max(both[formal_low:high+1]), active + values[0]/tau**s)
    for r in range(1, s+1):
        for u in range(1, r+1):
            valid = (z < s) & (u >= r-z) & (u <= s-z)
            continuation_bound = np.zeros(len(f0))
            for h in range(1, u+1):
                forbidden_return = ((v == 1) & (w < t)) | ((v == 2) & (w < t-s+1))
                candidate = np.where(forbidden_return, allowed_nonzero[f1, d1, h],
                                     allowed[f1, d1, h])
                if h < s:
                    single_root = (v == 1) & (z > 0)
                    single_value = np.where(f1 == 1, values[h+1],
                                            np.where(d1 == 1, max(values[h:h+2]), values[h]))
                    candidate = np.where(single_root, single_value, candidate)
                continuation_bound = np.maximum(
                    continuation_bound, np.where(h >= u-d0, candidate, 0))
            candidate = tau**(-(f0+d0+r-u).astype(float))*continuation_bound
            if fresh_orientation:
                probabilities = np.asarray([float(root_intersection_probability(q, s, r, zz, u))
                                            for zz in range(s)])
                selected = probabilities[np.minimum(z, s-1)]
                result_by_type[r] += np.where(valid, candidate*selected, 0)
            else:
                result_by_type[r] = np.maximum(result_by_type[r], np.where(valid, candidate, 0))
    index = f*(t+1)+d
    weights = multiplicities*y**v
    result = np.empty((t+1, t+1, s+1))
    for r in range(s+1):
        result[:, :, r] = np.bincount(index, weights=weights*result_by_type[r],
                                     minlength=(t+1)**2).reshape(t+1, t+1)
    if normalized:
        for fresh in range(t+1):
            for dep in range(t-fresh+1):
                result[fresh, dep] /= comb(t, fresh)*comb(t-fresh, dep)
    return result


def packet_action_exact(tau, y, continuation, q=256, t=16, s=4, normalized=False,
                        fresh_orientation=False):
    """Independent rational reference for the eta-scaled Bellman action."""
    from fractions import Fraction
    tau, y = Fraction(tau), Fraction(y)
    values = list(map(Fraction, continuation))
    if len(values) != s+1 or min(values) < 0 or not 1 <= tau <= q-1 or not 0 < y <= 1:
        raise ValueError('invalid continuation or tilts')
    zero = Fraction(0)
    nonzero = [zero] + values[1:]
    both = [nonzero[k]+values[0]/tau**k for k in range(s+1)]

    def step(payoff):
        out = [payoff[1]]
        for k in range(1, s+1):
            ordinary = max(payoff[k:k+2])
            out.append(ordinary+max(zero, payoff[k-1]-ordinary)/(q-1))
        return out

    allowed, allowed_nonzero = {}, {}
    for fresh in range(t+1):
        g = [max(nonzero[max(k, min(s, fresh)):min(s, k+fresh)+1]) for k in range(s+1)]
        b = [max(both[max(k, min(s, fresh)):min(s, k+fresh)+1]) for k in range(s+1)]
        for dep in range(t-fresh+1):
            a = fresh+dep
            allowed_nonzero[fresh, dep] = g
            allowed[fresh, dep] = ([min(b[k], g[k]+values[0]/tau**min(a, s))
                                   for k in range(s+1)] if a else nonzero)
            g, b = step(g), step(b)
    out = [[[zero for _ in range(s+1)] for _ in range(t+1)] for _ in range(t+1)]
    for raw in _support_types(t)[0].T:
        f0, d0, f1, d1, o, z = map(int, raw)
        f, d, v = f0+f1, d0+d1, f1+d1+o
        w = f+d
        factor = factorial(t)
        for count in (f0, d0, f1, d1, o, z):
            factor //= factorial(count)
        factor *= y**v
        local = [zero]*(s+1)
        if f0 == d0 == o == 0:
            if not w:
                local[0] = values[0]
            else:
                active = max(nonzero[max(1, min(s, f1)):min(s, w)+1])
                local[0] = (active if w <= s else
                    min(max(both[min(s, f1):min(s, w)+1]), active+values[0]/tau**s))
        if z < s:
            for r in range(1, s+1):
                for u in range(max(1, r-z), min(r, s-z)+1):
                    forbidden_return = (v == 1 and w < t) or (v == 2 and w < t-s+1)
                    selected = allowed_nonzero if forbidden_return else allowed
                    if v == 1 and z:
                        payoff = max(values[h+1] if f1 else
                                     (max(values[h:h+2]) if d1 else values[h])
                                     for h in range(max(1, u-d0), u+1))
                    else:
                        payoff = max(selected[f1, d1][max(1, u-d0):u+1])
                    candidate = payoff/tau**(r-u+f0+d0)
                    if fresh_orientation:
                        local[r] += root_intersection_probability(q, s, r, z, u)*candidate
                    else:
                        local[r] = max(local[r], candidate)
        for r in range(s+1):
            out[f][d][r] += factor*local[r]
    if normalized:
        for f in range(t+1):
            for d in range(t-f+1):
                out[f][d] = [value/(comb(t, f)*comb(t-f, d)) for value in out[f][d]]
    return out


def fixed_witness_diagnostics(iterations=120, refresh_packets=0, witnesses=None):
    """Fixed-witness, all-dependent comparison; no broad parameter search."""
    import math
    from scipy.special import logsumexp
    from return_cover import return_screen_objective
    if witnesses is None:
        witnesses = (
        (.01, [-3.0258570984074526, -.5857953682093398, 2.5381883850844997,
               2.014679343059534, 0]),
        (.10, [.70708872410493, -1.9100611711286484, 2.6167343763275186,
               -.2683299405157173, 0]),
        (.18, [2.1865170506895493, -3.3950117622093154, 3.492359740395451,
               -1.3396605412597715, 0]))
    results = []
    for alpha, logs in witnesses:
        lx, ly, lt, ll, le = logs
        x, y, tau = map(math.exp, (lx, ly, lt))
        powers = x**np.arange(17)

        def single(vector, fresh_orientation=False):
            packet = packet_action(tau, y, vector, fresh_orientation=fresh_orientation)
            return powers @ packet[0]

        def action(vector):
            if not refresh_packets:
                return single(vector)
            # Backward continuation: the three later plain packets are
            # applied before the first packet's fresh-orientation action.
            for _ in range(refresh_packets-1):
                vector = single(vector)
            return single(vector, fresh_orientation=True)

        vector = np.ones(5)
        for _ in range(iterations):
            updated = action(vector)
            updated /= updated[0]
            if np.max(np.abs(np.log(updated/vector))) < 1e-12:
                vector = updated
                break
            vector = (vector+updated)/2
        rho = max(action(vector)/vector)
        row = logsumexp([32*math.log1p(-alpha), -10*lt-11*le
            +32*np.logaddexp(math.log1p(-alpha), lt+ll+le+math.log(alpha))])
        combined = (math.log(rho)/(16*max(1, refresh_packets))
                    -alpha*lx-.1*ly-alpha*ll+row/32)/math.log(2)
        results.append(dict(alpha=alpha, law='eta-scaled',
            source_freshness_used=False, directed_certificate=False,
            refresh_symbols=16*refresh_packets,
            finite_boundary_factors_included=False,
            return_only_per_N=float(return_screen_objective(logs, alpha)),
            bellman_per_N=float(combined), rho=float(rho),
            continuation=vector.tolist(), logarithmic_parameters=logs))
    return results


if __name__ == '__main__':
    import argparse
    import json
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--diagnostics', action='store_true')
    parser.add_argument('--refresh-packets', type=int, nargs='+', default=(0,),
                        choices=(0, 1, 4))
    parser.add_argument('--refresh-tuned', action='store_true',
                        help='use root two-state refresh witnesses, without retuning their scalars')
    args = parser.parse_args()
    if args.diagnostics:
        for refresh in args.refresh_packets:
            tuned = None
            if args.refresh_tuned:
                tuned = {1: ((.01, [-2.591153, -.1251895, 3.6061255, 1.526953, 0]),
                             (.10, [1.033507, -1.664052, 2.697537, -.304106, 0]),
                             (.18, [2.7141905, -3.3970089, 3.566338, -1.379263, 0])),
                         4: ((.10, [.8524831, -1.8988915, 2.7512284, -.3283413, 0]),
                             (.18, [2.4570803, -3.4735652, 3.7083624, -1.4563104, 0]))}[refresh]
            print(json.dumps(fixed_witness_diagnostics(refresh_packets=refresh, witnesses=tuned),
                             indent=2), flush=True)
    else:
        print(json.dumps({'coefficient_terms': len(coefficient_table()),
                          'states': 5, 'q': 256, 't': 16,
                          'interface': 'packet_matrices(tau,y,normalized=False)'}))
