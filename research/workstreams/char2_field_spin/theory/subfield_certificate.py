"""Partial certificates and whole-code bound development for GF256 SPIN.

--certify rigorously covers exactly one and two active outer rows. It is
NOT a whole-code certificate. Screen modes use floating-point arithmetic.
"""

import argparse
import json
import math
from dataclasses import dataclass
from decimal import Context, Decimal, ROUND_CEILING, ROUND_FLOOR
from fractions import Fraction

import numpy as np

from q1_batched_screen import fixed_support_multiplicity, return_weights


LOW = Context(prec=80, rounding=ROUND_FLOOR)
HIGH = Context(prec=80, rounding=ROUND_CEILING)


@dataclass(frozen=True)
class Interval:
    lo: Decimal
    hi: Decimal

    @staticmethod
    def integer(value):
        return Interval(Decimal(value), Decimal(value))

    @staticmethod
    def rational(numerator, denominator=1):
        return Interval(LOW.divide(Decimal(numerator), Decimal(denominator)),
                        HIGH.divide(Decimal(numerator), Decimal(denominator)))

    def __add__(self, other):
        if not isinstance(other, Interval):
            other = Interval.integer(other)
        return Interval(LOW.add(self.lo, other.lo), HIGH.add(self.hi, other.hi))

    __radd__ = __add__

    def __sub__(self, other):
        if not isinstance(other, Interval):
            other = Interval.integer(other)
        return Interval(LOW.subtract(self.lo, other.hi), HIGH.subtract(self.hi, other.lo))

    def __rsub__(self, other):
        return Interval.integer(other)-self

    def __mul__(self, other):
        if not isinstance(other, Interval):
            other = Interval.integer(other)
        low_products = [LOW.multiply(a, b) for a in (self.lo, self.hi) for b in (other.lo, other.hi)]
        high_products = [HIGH.multiply(a, b) for a in (self.lo, self.hi) for b in (other.lo, other.hi)]
        return Interval(min(low_products), max(high_products))

    __rmul__ = __mul__

    def __truediv__(self, other):
        if not isinstance(other, Interval):
            other = Interval.integer(other)
        if other.lo <= 0 <= other.hi:
            raise ValueError('division interval contains zero')
        lows = [LOW.divide(a, b) for a in (self.lo, self.hi) for b in (other.lo, other.hi)]
        highs = [HIGH.divide(a, b) for a in (self.lo, self.hi) for b in (other.lo, other.hi)]
        return Interval(min(lows), max(highs))

    def __pow__(self, exponent):
        if exponent < 0:
            return Interval.integer(1)/(self**(-exponent))
        out, power = Interval.integer(1), self
        while exponent:
            if exponent & 1:
                out = out*power
            exponent >>= 1
            if exponent:
                power = power*power
        return out


ZERO, ONE = Interval.integer(0), Interval.integer(1)


def interval_local_envelope(z, width=16, q=256):
    """Entrywise upper envelope, with outward intervals around its entries."""
    matrices = []
    f = Interval.rational(1, q-1)+Interval.rational(q-2, q-1)*z
    for w in range(width+1):
        if w == 0:
            p0 = ONE
        else:
            p0 = Interval.rational((q-1)**w+(-1)**w*(q*q-1-w*q*(q-1)),
                                   q*q*(q-1)**w)
        g = Interval.rational(q+1-width, q+1)*z**(width-w)*f**w
        if w:
            g += Interval.rational(w, q+1)*z**(width-w+1)*f**(w-1)
        if w < width:
            g += Interval.rational(width-w, q+1)*z**(width-w-1)*f**w
        ret = (1-p0)/(q*q-1)
        jlo, jhi = ret*z**width, ret*z**max(0, width-1-w)
        if w == 1:
            jlo = jhi = sum((z**d for d in return_weights(width, q)), ZERO)/(width*(q*q-1))
        matrices.append([[z**w*p0, z**w*(1-p0)], [jhi, g-jlo]])
    return matrices


def interval_matrix_product(left, right):
    return [[sum((left[i][k]*right[k][j] for k in range(2)), ZERO)
             for j in range(2)] for i in range(2)]


def interval_polynomial_product(left, right, limit):
    size = min(limit+1, len(left)+len(right)-1)
    out = [[[ZERO, ZERO], [ZERO, ZERO]] for _ in range(size)]
    for i, li in enumerate(left):
        for j, rj in enumerate(right[:max(0, size-i)]):
            product = interval_matrix_product(li, rj)
            for row in range(2):
                for col in range(2):
                    out[i+j][row][col] += product[row][col]
    return out


def interval_regions(z, maximum_h, rows=8192, width=16):
    local = interval_local_envelope(z, width)
    local = [[[entry*math.comb(width, w) for entry in row] for row in matrix]
             for w, matrix in enumerate(local[:maximum_h+1])]
    out = [[[ONE, ZERO], [ZERO, ONE]]]
    exponent = rows//width
    while exponent:
        if exponent & 1:
            out = interval_polynomial_product(out, local, maximum_h)
        exponent >>= 1
        if exponent:
            local = interval_polynomial_product(local, local, maximum_h)
    return [[[entry/math.comb(rows, h) for entry in row] for row in matrix]
            for h, matrix in enumerate(out)]


def interval_q1(z):
    """Exact MDS supports; interval-round the one-row expectation."""
    regions = interval_regions(z, 1)
    states = [[ONE, ZERO]]+[[ZERO, ZERO] for _ in range(16)]
    for region in range(16):
        next_states = [[ZERO, ZERO] for _ in range(17)]
        for w in range(region+1):
            for active in range(2):
                for end in range(2):
                    next_states[w+active][end] += sum(
                        (states[w][start]*regions[active][start][end] for start in range(2)), ZERO)
        states = next_states
    total = ZERO
    for w in range(9, 17):
        count = fixed_support_multiplicity(16, 8, w, 256)//255
        total += (states[w][0]+states[w][1])*count
    return 8192*total/(z**13107)


def interval_coupled(z, active_rows):
    """Coupled regional bound with n_w <= q^(w-8) source envelope."""
    regions = interval_regions(z, active_rows)
    matrix = [[sum((regions[h][i][j]*math.comb(active_rows, h)*256**h
                    for h in range(active_rows+1)), ZERO)
               for j in range(2)] for i in range(2)]
    vector = [ONE, ZERO]
    for region in range(16):
        vector = [sum((vector[i]*matrix[i][j] for i in range(2)), ZERO) for j in range(2)]
    return ((vector[0]+vector[1])*math.comb(8192, active_rows)
            / (255*256**(8*active_rows)*z**13107))


def integer_margin(upper):
    exact = Fraction(upper)
    if exact <= 0:
        raise ValueError('a positive upper bound is required')
    bits = exact.denominator.bit_length()-exact.numerator.bit_length()
    while exact > Fraction(2)**(-bits):
        bits -= 1
    while exact <= Fraction(2)**(-bits-1):
        bits += 1
    return bits


def window_obstruction(active_rows=32, state_size=2, window=816):
    """Exact first-moment lower bound; NOT a failure-probability bound."""
    if not 0 < active_rows <= window <= 8192 or window % 16:
        raise ValueError('invalid aligned window')
    dimension = 8*active_rows-15*state_size
    if dimension <= 0 or 16*window > 13107:
        raise ValueError('the proposed window does not give a nonzero bad subspace')
    event_probability = Fraction(math.comb(window, active_rows), math.comb(8192, active_rows))**16
    bad_lines = (256**dimension-1)//255
    expected_lower = event_probability*bad_lines
    return dict(active_rows=active_rows, state_size=state_size, window=window,
                maximum_output_weight=16*window, subspace_dimension_lower=dimension,
                route_event_probability=event_probability,
                expected_bad_lines_lower=expected_lower)


def certify_sparse():
    entries = []
    total = ZERO
    for active, numerator, denominator in [(1, 8183, 8192), (2, 511, 512)]:
        z = Interval.rational(numerator, denominator)
        bound = interval_q1(z) if active == 1 else interval_coupled(z, active)
        total += bound
        entries.append(dict(active_rows=active, z=[numerator, denominator],
                            failure_upper=str(bound.hi), certified_integer_bits=integer_margin(bound.hi)))
    return dict(format='gf256-batched-sparse-v1', complete=False,
                config=dict(q=256, k=65536, n=131072, outer=[16, 8], width=16, state=2,
                            bad_weight_threshold=13107, refresh='uniform-sl2-pre-output',
                            diagonal='independent-uniform-nonzero'),
                covered_active_rows=[1, 2], uncovered_active_rows=[3, 8192], entries=entries,
                sparse_failure_upper=str(total.hi), certified_integer_bits=integer_margin(total.hi))


def verify_sparse(record):
    expected = certify_sparse()
    for key in ['format', 'complete', 'config', 'covered_active_rows', 'uncovered_active_rows']:
        if record.get(key) != expected[key]:
            raise ValueError(f'incorrect scope or configuration: {key}')
    if len(record.get('entries', [])) != 2:
        raise ValueError('incomplete sparse coverage')
    for supplied, recomputed in zip(record['entries'], expected['entries']):
        for key in ['active_rows', 'z']:
            if supplied.get(key) != recomputed[key]:
                raise ValueError(f'incorrect witness: {key}')
        upper = Decimal(supplied['failure_upper'])
        if not upper.is_finite() or upper < Decimal(recomputed['failure_upper']):
            raise ValueError('claimed class bound is too small')
        if supplied['certified_integer_bits'] > integer_margin(upper):
            raise ValueError('claimed class margin is too large')
    upper = Decimal(record['sparse_failure_upper'])
    supplied_sum = sum((Fraction(Decimal(e['failure_upper'])) for e in record['entries']), Fraction(0))
    if not upper.is_finite() or Fraction(upper) < supplied_sum:
        raise ValueError('claimed sparse total is too small')
    if record['certified_integer_bits'] > integer_margin(upper):
        raise ValueError('claimed sparse margin is too large')
    return True


def logsumexp(values):
    maximum = max(values)
    return maximum + math.log(sum(math.exp(v-maximum) for v in values))


def local_envelope(z, width=16, q=256):
    result = np.zeros((width+1, 2, 2))
    f = 1/(q-1) + (q-2)/(q-1)*z
    for w in range(width+1):
        if w == 0:
            p0 = 1.0
        else:
            p0 = (1 + (-1)**w * (q*q-1-w*q*(q-1))/(q-1)**w)/(q*q)
        g = (q+1-width)/(q+1)*z**(width-w)*f**w
        if w:
            g += w/(q+1)*z**(width-w+1)*f**(w-1)
        if w < width:
            g += (width-w)/(q+1)*z**(width-w-1)*f**w
        ret = (1-p0)/(q*q-1)
        jlo = ret*z**width
        jhi = ret*z**max(0, width-1-w)
        if w == 1:
            jlo = jhi = sum(z**d for d in return_weights(width, q))/(width*(q*q-1))
        result[w] = [[z**w*p0, z**w*(1-p0)], [jhi, g-jlo]]
    return result


def polynomial_product(left, right, limit):
    size = min(limit+1, len(left)+len(right)-1)
    out = np.zeros((size, 2, 2))
    for i in range(2):
        for j in range(2):
            for k in range(2):
                out[:, i, j] += np.convolve(left[:, i, k], right[:, k, j])[:size]
    return out


def regional_envelopes(z, maximum_h, rows=8192, width=16, q=256):
    assert rows % width == 0
    p = maximum_h/rows
    local = local_envelope(z, width, q)
    for w in range(width+1):
        local[w] *= math.comb(width, w)*p**w*(1-p)**(width-w)
    result = np.eye(2)[None, :, :]
    exponent = rows//width
    while exponent:
        if exponent & 1:
            result = polynomial_product(result, local, maximum_h)
        exponent >>= 1
        if exponent:
            local = polynomial_product(local, local, maximum_h)
    for h in range(maximum_h+1):
        logprob = (math.lgamma(rows+1)-math.lgamma(h+1)-math.lgamma(rows-h+1)
                   + h*math.log(p)+(rows-h)*math.log1p(-p))
        result[h] /= math.exp(logprob)
    return result


def norm_screen(active_rows, k=2**16, width=16, outer_length=16,
                outer_dimension=8, delta=0.1, search_steps=70):
    rows = k//outer_dimension
    n = rows*outer_length
    threshold = math.ceil(delta*n)-1
    logchoose = math.lgamma(rows+1)-math.lgamma(active_rows+1)-math.lgamma(rows-active_rows+1)
    outer_counts = [(w, fixed_support_multiplicity(outer_length, outer_dimension, w, 256)
                    * math.comb(outer_length, w))
                    for w in range(outer_length-outer_dimension+1, outer_length+1)]
    best = None
    for iz in range(search_steps+1):
        tilt = math.exp(2+5*iz/search_steps)
        logz = -tilt*active_rows/n
        z = math.exp(logz)
        regions = regional_envelopes(z, active_rows, rows, width)
        for ie in range(81):
            eta = math.exp(-10+10*ie/80)
            row0 = regions[1:, 0, 0]+eta*regions[1:, 0, 1]
            row1 = regions[1:, 1, 0]/eta+regions[1:, 1, 1]
            logrhos = np.log(np.maximum(row0, row1))/np.arange(1, active_rows+1)
            logrho = float(np.max(logrhos))
            logp = logsumexp([math.log(count)+w*logrho for w, count in outer_counts])
            logbound = logchoose+active_rows*logp-math.log(255)-math.log(eta)-threshold*logz
            bits = -logbound/math.log(2)
            if best is None or bits > best['margin_bits']:
                best = dict(active_rows=active_rows, margin_bits=bits, z=z,
                            eta=eta, rho=math.exp(logrho), worst_h=int(np.argmax(logrhos))+1)
    return best


def coupled_screen(active_rows, k=2**16, width=16, outer_length=16,
                   outer_dimension=8, delta=0.1, search_steps=70):
    rows = k//outer_dimension
    n = rows*outer_length
    threshold = math.ceil(delta*n)-1
    logchoose = math.lgamma(rows+1)-math.lgamma(active_rows+1)-math.lgamma(rows-active_rows+1)
    best = None
    for iz in range(search_steps+1):
        tilt = math.exp(1+6*iz/search_steps)
        logz = -tilt*active_rows/n
        z = math.exp(logz)
        regions = regional_envelopes(z, active_rows, rows, width)
        logs = np.array([math.log(math.comb(active_rows, h))+(h-active_rows/2)*math.log(256)
                         for h in range(active_rows+1)])
        # Common normalization controls intermediate magnitudes; add it
        # back after taking the product across outer_length regions.
        scale = float(np.max(logs))
        matrix = np.sum(regions*np.exp(logs-scale)[:, None, None], axis=0)
        vector = np.array([1.0, 0.0])
        cumulative = outer_length*scale
        for region in range(outer_length):
            vector = vector @ matrix
            normalization = float(np.max(vector))
            if normalization == 0:
                break
            vector /= normalization
            cumulative += math.log(normalization)
        if normalization == 0:
            continue
        logbound = (logchoose-math.log(255)-threshold*logz
                    +cumulative+math.log(float(np.sum(vector))))
        bits = -logbound/math.log(2)
        if best is None or bits > best['margin_bits']:
            best = dict(active_rows=active_rows, margin_bits=bits, z=z)
    return best


def exact_sparse_screen(active_rows, k=2**16, width=16, outer_length=16,
                        outer_dimension=8, delta=0.1, search_steps=30):
    """Keep each selected outer row's support weight (exponential in a)."""
    assert 1 <= active_rows <= 4
    rows = k//outer_dimension
    n = rows*outer_length
    threshold = math.ceil(delta*n)-1
    logchoose = math.lgamma(rows+1)-math.lgamma(active_rows+1)-math.lgamma(rows-active_rows+1)
    factors = np.array([fixed_support_multiplicity(outer_length, outer_dimension, w, 256)
                        / 256**(w-outer_dimension) if w >= outer_length-outer_dimension+1 else 0
                        for w in range(outer_length+1)])
    masks = [tuple((mask >> j) & 1 for j in range(active_rows)) for mask in range(1 << active_rows)]
    best = None
    for iz in range(search_steps+1):
        tilt = math.exp(2+4*iz/search_steps)
        logz = -tilt*active_rows/n
        z = math.exp(logz)
        regions = regional_envelopes(z, active_rows, rows, width)
        matrices = [regions[h]*256**(h-active_rows/2) for h in range(active_rows+1)]
        current = np.zeros((1,)*active_rows+(2,))
        current[(0,)*active_rows+(0,)] = 1
        cumulative = 0.0
        for region in range(outer_length):
            next_state = np.zeros((region+2,)*active_rows+(2,))
            for mask in masks:
                target = tuple(slice(bit, region+1+bit) for bit in mask)+(slice(None),)
                next_state[target] += current @ matrices[sum(mask)]
            scale = float(np.max(next_state))
            next_state /= scale
            cumulative += math.log(scale)
            current = next_state
        for axis in range(active_rows):
            shape = [1]*(active_rows+1)
            shape[axis] = outer_length+1
            current *= factors.reshape(shape)
        total = float(np.sum(current))
        logbound = logchoose-math.log(255)-threshold*logz+cumulative+math.log(total)
        bits = -logbound/math.log(2)
        if best is None or bits > best['margin_bits']:
            best = dict(active_rows=active_rows, margin_bits=bits, z=z)
    return best


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--screen', type=int, nargs='+', default=[2, 3, 4, 8, 16])
    parser.add_argument('--search-steps', type=int, default=70)
    parser.add_argument('--coupled', action='store_true')
    parser.add_argument('--exact-sparse', action='store_true')
    parser.add_argument('--certify', action='store_true')
    args = parser.parse_args()
    if args.certify:
        print(json.dumps(certify_sparse(), indent=2))
        raise SystemExit(0)
    for active in args.screen:
        fn = exact_sparse_screen if args.exact_sparse else coupled_screen if args.coupled else norm_screen
        print(json.dumps(fn(active, search_steps=args.search_steps)), flush=True)
