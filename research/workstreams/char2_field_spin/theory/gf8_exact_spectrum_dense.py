"""Exact-MDS-spectrum dense bounds with checked Bernoulli homogenization.

All support counts and convexity checks are exact. Floating point only
proposes tilts; Arb evaluates acceptance. This is a separate refinement,
leaving the old two-product dense bound available as a fallback.
"""
from __future__ import annotations

from fractions import Fraction as F
from functools import lru_cache
from hashlib import sha256
from math import comb
import math

import numpy as np
from scipy.optimize import minimize_scalar
from scipy.special import logsumexp
from flint import arb, ctx

from gf8_outer_rate_dense import rational, lower_integer


@lru_cache(maxsize=64)
def support_counts(length=128, dimension=74, q=256):
    """Exact number on each prescribed support of a linear MDS code."""
    if not 1 <= dimension <= length <= q or q < 2:
        raise ValueError('invalid MDS geometry')
    parity = length-dimension
    result = [1]+[0]*length
    for weight in range(parity+1, length+1):
        rank = weight-parity
        result[weight] = sum((-1)**j*comb(weight, j)*(q**(rank-j)-1)
                             for j in range(rank))
        if result[weight] <= 0:
            raise ArithmeticError('nonpositive MDS support count')
    if sum(comb(length, w)*value for w, value in enumerate(result)) != q**dimension:
        raise ArithmeticError('MDS spectrum does not sum to code size')
    return tuple(result)


def marked_coefficients(tilt, length=128, dimension=74, q=256):
    tilt = F(tilt)
    if tilt < 1:
        raise ValueError('tilt below one')
    ratio = tilt/F(q-1)
    return tuple(count*ratio**w for w, count in enumerate(support_counts(length, dimension, q)))


def convexity_slacks(tilt, length=128, dimension=74, q=256):
    a = marked_coefficients(tilt, length, dimension, q)
    return tuple(a[w+2]-2*a[w+1]+a[w] for w in range(length-1))


def convex_tilt_floor(length=128, dimension=74, q=256):
    """Rational proposal above every upper quadratic root; checked exactly."""
    counts = support_counts(length, dimension, q)
    threshold = 1.
    for w in range(length-1):
        if not counts[w+2]:
            continue
        # Divide a[w+2]*lambda^2-2*a[w+1]*lambda+a[w] by
        # its positive leading coefficient before converting to float.
        b = float(F(counts[w+1]*(q-1), counts[w+2]))
        c = float(F(counts[w]*(q-1)**2, counts[w+2]))
        discriminant = b*b-c
        if discriminant >= 0:
            threshold = max(threshold, b+math.sqrt(discriminant))
    candidate = F(math.ceil((threshold+1e-10)*10**6), 10**6)
    while min(convexity_slacks(candidate, length, dimension, q)) < 0:
        candidate += F(1, 10**6)
    return candidate


def bracket(lower, upper, tilt, length=128, dimension=74, q=256):
    """Uniform positive Bernstein-term upper bound on one density cell."""
    lower, upper, tilt = map(F, (lower, upper, tilt))
    if not 0 <= lower <= upper <= 1:
        raise ValueError('invalid exact density cell')
    coefficients = marked_coefficients(tilt, length, dimension, q)
    if min(convexity_slacks(tilt, length, dimension, q)) < 0:
        raise ValueError('marked support sequence is not discretely convex')
    result = F(0)
    for w, coefficient in enumerate(coefficients):
        if coefficient:
            point = min(upper, max(lower, F(w, length)))
            result += comb(length, w)*coefficient*point**w*(1-point)**(length-w)
    return result


@lru_cache(maxsize=64)
def positive_excess(length=128, dimension=74, q=256):
    """Nonnegative excess above the constant random-subspace membership."""
    baseline = F(1, q**(length-dimension))
    return tuple(max(F(0), F(count, (q-1)**w)-baseline)
                 for w, count in enumerate(support_counts(length, dimension, q)))


def excess_bracket(lower, upper, tilt, length=128, dimension=74, q=256):
    """Exact endpoint bound using Maclaurin only on forced-zero sets.

    For each excess support term, discard the factors alpha_j on live
    positions. Summing the remaining products gives an elementary
    symmetric polynomial in 1-alpha_j. Maclaurin bounds it by its value
    at the common mean. The constant baseline uses ordinary AM-GM.
    """
    lower, upper, tilt = map(F, (lower, upper, tilt))
    if not 0 <= lower <= upper <= 1 or tilt < 1:
        raise ValueError('invalid exact density cell or tilt')
    baseline = F(1, q**(length-dimension))
    result = baseline*(1+(tilt-1)*upper)**length
    for w, excess in enumerate(positive_excess(length, dimension, q)):
        if excess:
            result += comb(length, w)*excess*tilt**w*(1-lower)**(length-w)
    return result


def propose_tilt(lower, upper, length=128, dimension=74, q=256):
    """Floating point proposes only; returned witness is checked exactly."""
    lower, upper = float(lower), float(upper)
    floor = convex_tilt_floor(length, dimension, q)
    counts = support_counts(length, dimension, q)
    terms, weights = [], []
    for w, count in enumerate(counts):
        if not count:
            continue
        point = min(upper, max(lower, w/length))
        if (point == 0 and w) or (point == 1 and w < length):
            continue
        logterm = math.log(comb(length, w))+math.log(count)-w*math.log(q-1)
        if w:
            logterm += w*math.log(point)
        if w < length:
            logterm += (length-w)*math.log1p(-point)
        terms.append(logterm)
        weights.append(w)
    terms, weights = np.asarray(terms), np.asarray(weights)
    objective = lambda x: float(logsumexp(terms+weights*x))-length*lower*x
    result = minimize_scalar(objective, bounds=(math.log(float(floor)), 20), method='bounded')
    candidates = [floor, max(floor, F(math.exp(float(result.x))).limit_denominator(10**6))]
    proposal = min(candidates, key=lambda x: objective(math.log(float(x))))
    if min(convexity_slacks(proposal, length, dimension, q)) < 0:
        raise ArithmeticError('unsafe proposed tilt')
    return proposal


def propose_excess_tilt(lower, upper, length=128, dimension=74, q=256):
    lower, upper = float(lower), float(upper)
    terms, weights = [], []
    for w, excess in enumerate(positive_excess(length, dimension, q)):
        if not excess or (lower == 1 and w < length):
            continue
        term = math.log(comb(length, w))+math.log(excess.numerator)-math.log(excess.denominator)
        if w < length:
            term += (length-w)*math.log1p(-lower)
        terms.append(term)
        weights.append(w)
    terms, weights = np.asarray(terms), np.asarray(weights)
    def objective(logtilt):
        baseline = -(length-dimension)*math.log(q)+length*math.log1p(math.expm1(logtilt)*upper)
        tail = float(logsumexp(terms+weights*logtilt)) if len(terms) else -math.inf
        return float(np.logaddexp(baseline, tail))-length*lower*logtilt
    result = minimize_scalar(objective, bounds=(0, 20), method='bounded')
    candidates = [F(1), max(F(1), F(math.exp(float(result.x))).limit_denominator(10**6))]
    return min(candidates, key=lambda x: objective(math.log(float(x))))


@lru_cache(maxsize=16)
def check_ball_exponent(exponent):
    exponent = F(exponent)
    a, b = exponent.numerator, exponent.denominator
    divisor = math.gcd(10, b)
    if exponent <= 0 or not F(10, 9)**(10*b//divisor)*2295**(b//divisor) < 2**(10*a//divisor):
        raise ValueError('proposed output-ball exponent is not certified')
    return exponent


def interval_margin(g, first, last, tilt, ball_exponent=F(127, 100), method='excess'):
    """Directed margin for a fixed or periodic full-row-rank puncture."""
    if not 1 <= first <= last <= g['input_n']:
        raise ValueError('invalid inverse-input interval')
    ball_exponent = check_ball_exponent(ball_exponent)
    n = g['input_n']
    length, dimension = g['outer']
    if method not in ('excess', 'convex'):
        raise ValueError('unknown exact-spectrum method')
    upper_bound = excess_bracket if method == 'excess' else bracket
    base = upper_bound(F(first, n), F(last, n), tilt, length, dimension)
    return ((first*rational(tilt).log()-g['rows']*rational(base).log()
             -length*arb(g['rows']+1).log())/arb(2).log()
            -rational(ball_exponent)*g['output_n']-8*(n-g['output_n']))


def point_screen(g, density, ball_exponent=F(127, 100), method='excess'):
    density = F(density)
    n = g['input_n']
    weight = min(n, max(1, (n*density.numerator+density.denominator-1)//density.denominator))
    proposer = propose_excess_tilt if method == 'excess' else propose_tilt
    tilt = proposer(F(weight, n), F(weight, n), *g['outer'])
    previous = ctx.prec
    try:
        ctx.prec = 192
        bits = lower_integer(interval_margin(g, weight, weight, tilt, ball_exponent, method))
    finally:
        ctx.prec = previous
    return dict(outer=g['outer'], input_n=n, output_n=g['output_n'], weight=weight,
                density=float(F(weight, n)), tilt=str(tilt), margin_bits_lower=bits,
                method=method, exact_convexity_checked=method == 'convex',
                rigorous_point=True, complete=False, whole_code=False)


def certify_profile(g, first_density=F(9, 10), cell_bits=128, initial_grid=200,
                    precision=192, include_witnesses=False):
    """Cover every dense integer weight; caller supplies construction identity.

    The receipt assumes the stated independent full-row-rank inner and
    GF256 regional ensemble. It does not itself prove those hypotheses.
    """
    if (not isinstance(first_density, F) or not 0 < first_density <= 1
            or cell_bits < 40 or initial_grid < 1 or precision < 80):
        raise ValueError('invalid exact replay parameters')
    if (g['outer'][0] != 128 or g.get('q', 256) != 256
            or g['input_n'] != 128*g['rows']
            or not 0 < g['output_n'] <= g['input_n']):
        raise ValueError('invalid GF256 dense geometry')
    if not F(10, 9)**100*2295**10 < 2**127:
        raise ArithmeticError('output-ball rational comparison failed')
    n = g['input_n']
    first = (n*first_density.numerator+first_density.denominator-1)//first_density.denominator
    previous = ctx.prec
    try:
        ctx.prec = precision
        span = (n+initial_grid-1)//initial_grid
        pending = [(lo, min(n, lo+span-1)) for lo in range(first, n+1, span)][::-1]
        entries = []
        while pending:
            lo, hi = pending.pop()
            tilt = propose_excess_tilt(F(lo, n), F(hi, n), *g['outer'])
            bits = lower_integer(interval_margin(g, lo, hi, tilt))
            if bits >= cell_bits:
                entries.append((lo, hi, tilt, bits))
            elif lo == hi:
                raise ArithmeticError(f'exact-spectrum dense weight {lo} has only {bits} bits')
            else:
                mid = (lo+hi)//2
                pending.extend(((mid+1, hi), (lo, mid)))
        if (entries[0][0] != first or entries[-1][1] != n
                or any(a[1]+1 != b[0] for a, b in zip(entries, entries[1:]))):
            raise ArithmeticError('dense integer coverage gap')
        bits = min(entry[3] for entry in entries)
        witness_digest = sha256(''.join(f'{lo}:{hi}:{tilt}:{margin}\n'
                                       for lo, hi, tilt, margin in entries).encode('ascii')).hexdigest()
        spectrum_digest = sha256(','.join(map(str, support_counts(*g['outer']))).encode('ascii')).hexdigest()
        result = dict(g)
        result.update(weights=[first, n], covered_weights=n-first+1, first_density=str(first_density),
            checked_cells=len(entries), precision_bits=precision, initial_grid=initial_grid,
            target_cell_bits=cell_bits, ball_exponent='127/100', conditioning='(rows+1)^128',
            omitted_symbols=n-g['output_n'], completion_factor_log2=8*(n-g['output_n']),
            margin_bits_lower=bits, failure_upper_dyadic_exponent=bits,
            failure_upper_dyadic_hex=['0x1', hex(1 << bits)],
            witness_sha256=witness_digest, exact_support_spectrum_sha256=spectrum_digest,
            dense_method='exact-MDS-positive-excess-Maclaurin',
            complete=True, rigorous=True, whole_code=False,
            scope='complete dense inverse-input interval including omitted-output completions',
            dense_inner_hypothesis='invertible full causal inner; independent nonzero input diagonals',
            abstract_dense_inner_hypothesis='retained full-row-rank inner independent of regional route and nonzero input diagonals')
        if include_witnesses:
            result['witnesses'] = [[lo, hi, str(tilt), margin] for lo, hi, tilt, margin in entries]
        return result
    finally:
        ctx.prec = previous


if __name__ == '__main__':
    import argparse, json
    from gf8_puncture_frontier import profile
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--dimension', type=int, default=74)
    parser.add_argument('--schedule', type=int, nargs='+', default=[5, 6])
    parser.add_argument('--density', type=F, nargs='+', default=[F(99, 100), F(999, 1000), F(1)])
    parser.add_argument('--method', choices=('excess', 'convex'), default='excess')
    parser.add_argument('--certify', action='store_true')
    parser.add_argument('--first', type=F, default=F(9, 10))
    args = parser.parse_args()
    if args.certify:
        from gf8_periodic_climb_profile import metadata
        result = certify_profile(metadata(args.dimension, tuple(args.schedule)), args.first)
        result.pop('failure_upper_dyadic_hex')
        print(json.dumps(result), flush=True)
    else:
        g = profile(args.dimension, tuple(args.schedule))
        for density in args.density:
            print(json.dumps(point_screen(g, density, method=args.method)), flush=True)
