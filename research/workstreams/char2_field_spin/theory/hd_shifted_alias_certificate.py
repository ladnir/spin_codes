"""Directed targeted bounds for shifted pairs, using an interval radix-2 FFT.

Exact rational packet caps are raised to a rational alpha with Arb. The
three-variable polynomial is transformed using explicit binary64 interval
butterflies; neither numpy.fft nor a guessed numerical cushion is used.
Receipts cover one class full under the even matching, never the whole code.
"""
from fractions import Fraction as F
from functools import lru_cache
from itertools import combinations_with_replacement
from math import comb, factorial, isfinite, sqrt

import numpy as np
from flint import arb, arb_poly, ctx

from hd_common_refresh_packets import local_moments, matrix
from hd_grouped_quiet import quiet_cap
from hd_shifted_pairs_count import circular_run_count
from hd_sumrank_paired_certificate import poly_product, round_up, scalar_power


LIMIT = (1 << 21) * 7 // 20
DEFAULT_MARK = F('0.48180002666992056')
DEFAULT_ALPHA = F('0.09468173632710614')
# Positive rational markers are witnesses; their choice needs no proof.
DEFAULT_SINGLETON = F('0.049459')
DEFAULT_FULL = F('0.18676')


def rational(value):
    value = F(value)
    return arb(value.numerator) / value.denominator


def binary64_enclosure(value):
    """Return endpoints, checking containment with Arb after conversion."""
    if not value.is_finite():
        raise ArithmeticError('A nonfinite Arb input has no finite enclosure')
    if value.is_zero():
        return 0., 0.
    lo = float(np.nextafter(float(value.lower()), -np.inf))
    hi = float(np.nextafter(float(value.upper()), np.inf))
    if not isfinite(lo) or not isfinite(hi):
        raise ArithmeticError('Binary64 enclosure overflow')
    if not arb(lo) <= value.lower() or not value.upper() <= arb(hi):
        raise ArithmeticError('Arb-to-binary64 enclosure check failed')
    return lo, hi


def _down(value):
    return np.nextafter(value, -np.inf)


def _up(value):
    return np.nextafter(value, np.inf)


def _add(a, b):
    return _down(a[0] + b[0]), _up(a[1] + b[1])


def _sub(a, b):
    return _down(a[0] - b[1]), _up(a[1] - b[0])


def _mul(a, b):
    products = (a[0] * b[0], a[0] * b[1],
                a[1] * b[0], a[1] * b[1])
    return (_down(np.minimum.reduce(products)),
            _up(np.maximum.reduce(products)))


def _cadd(a, b):
    return _add(a[0], b[0]), _add(a[1], b[1])


def _csub(a, b):
    return _sub(a[0], b[0]), _sub(a[1], b[1])


def _cmul(a, b):
    return (_sub(_mul(a[0], b[0]), _mul(a[1], b[1])),
            _add(_mul(a[0], b[1]), _mul(a[1], b[0])))


def _map(function, value):
    return tuple(tuple(function(endpoint) for endpoint in part) for part in value)


def _slice(value, index):
    return _map(lambda endpoint: endpoint[index], value)


def _assign(target, index, value):
    for target_part, value_part in zip(target, value):
        for target_endpoint, value_endpoint in zip(target_part, value_part):
            target_endpoint[index] = value_endpoint


def _finite(value):
    if not all(np.all(np.isfinite(endpoint)) for part in value for endpoint in part):
        raise ArithmeticError('Interval transform overflow or invalid arithmetic')


def _roots(grid, inverse=False):
    sign = 1 if inverse else -1
    real = [binary64_enclosure((arb(sign * 2 * k) / grid).cos_pi()) for k in range(grid)]
    imag = [binary64_enclosure((arb(sign * 2 * k) / grid).sin_pi()) for k in range(grid)]
    return (tuple(np.array(row) for row in zip(*real)),
            tuple(np.array(row) for row in zip(*imag)))


def _check_binary64_environment():
    if np.dtype(float).itemsize != 8 or np.finfo(np.float64).nmant != 52:
        raise ArithmeticError('Expected IEEE-754 binary64 arrays')
    smallest = np.float64(float.fromhex('0x0.0000000000001p-1022'))
    normal = np.float64(float.fromhex('0x1p-1022'))
    if (F(float(np.multiply(normal, np.float64(.5)))) != F(1, 2**1023)
            or F(float(np.add(smallest, smallest))) != F(1, 2**1073)):
        raise ArithmeticError('Gradual underflow is required; FTZ/DAZ is unsupported')


def interval_fft(value, axes, inverse=False):
    """Unnormalized radix-2 DFT with rectangular complex interval outputs.

    Each real operation is a separate NumPy float64 ufunc, expanded by one
    representable number. This assumes IEEE-754 binary64 round-to-nearest
    with gradual underflow; invalid values and overflow fail closed.
    """
    _check_binary64_environment()
    if any(endpoint.dtype != np.float64 for part in value for endpoint in part):
        raise ValueError('Interval FFT endpoints must be binary64 arrays')
    for axis in axes:
        size = value[0][0].shape[axis]
        if size < 1 or size & (size - 1):
            raise ValueError('FFT axes must have power-of-two lengths')
        roots = _roots(size, inverse)
        bits = size.bit_length() - 1
        reverse = np.array([int(format(i, f'0{bits}b')[::-1], 2) for i in range(size)])
        value = _map(lambda endpoint: np.take(np.moveaxis(endpoint, axis, -1), reverse, axis=-1), value)
        original_shape = value[0][0].shape
        span = 2
        while span <= size:
            blocks = _map(lambda endpoint: endpoint.reshape(-1, size // span, span), value)
            left = _map(np.copy, _slice(blocks, (..., slice(0, span // 2))))
            right = _slice(blocks, (..., slice(span // 2, span)))
            twiddle = _slice(roots, slice(0, size // 2, size // span))
            rotated = _cmul(right, twiddle)
            _assign(blocks, (..., slice(0, span // 2)), _cadd(left, rotated))
            _assign(blocks, (..., slice(span // 2, span)), _csub(left, rotated))
            value = _map(lambda endpoint: endpoint.reshape(original_shape), blocks)
            span *= 2
        value = _map(lambda endpoint: np.moveaxis(endpoint, -1, axis), value)
        _finite(value)
    return value


def interval_matrix_power(value, steps):
    if steps < 1 or steps & (steps - 1) or value[0][0].shape[0] != 4:
        raise ValueError('Expected four matrix entries and power-of-two steps')
    for _ in range(steps.bit_length() - 1):
        a, b, c, d = (_slice(value, i) for i in range(4))
        trace, bc = _cadd(a, d), _cmul(b, c)
        entries = (_cadd(_cmul(a, a), bc), _cmul(b, trace),
                   _cmul(c, trace), _cadd(_cmul(d, d), bc))
        value = tuple(tuple(np.stack([entry[part][end] for entry in entries])
                            for end in range(2)) for part in range(2))
        _finite(value)
    return value


def selected_inverse(value, target):
    """Enclose one coefficient of the normalized inverse DFT."""
    shape = value[0][0].shape
    if len(shape) != 4 or len(target) != 3 or len(set(shape[1:])) != 1:
        raise ValueError('Expected a cubic three-variable transform')
    grid = shape[1]
    roots = _roots(grid, True)
    for axis, degree in enumerate(target, 1):
        if degree % grid:
            indices = np.arange(grid) * (degree % grid) % grid
            broadcast = [1, 1, 1, 1]
            broadcast[axis] = grid
            phases = _map(lambda endpoint: endpoint[indices].reshape(broadcast), roots)
            value = _cmul(value, phases)
    value = _map(lambda endpoint: endpoint.reshape(4, -1), value)
    while value[0][0].shape[1] > 1:
        value = _cadd(_slice(value, (slice(None), slice(0, None, 2))),
                      _slice(value, (slice(None), slice(1, None, 2))))
    factor = 1. / (grid ** 3)  # Exact binary power; multiplication is enclosed.
    value = tuple(_mul(part, (factor, factor)) for part in value)
    _finite(value)
    return _map(lambda endpoint: endpoint[:, 0], value)


@lru_cache(2)
def packet_table(z, bits=60):
    """Exact upward dyadic h4 caps, including diagonal and quiet factors."""
    z = F(z)
    if not 0 < z < 1 or bits < 16:
        raise ValueError('Invalid exact packet witness')
    local = local_moments(z, q=256, h=4)
    density = F(256 ** 16 - 1, 255 ** 16)
    quiet = quiet_cap(z, q=256, h=4)
    result = {}
    for occupancy in combinations_with_replacement(range(17), 4):
        k, birth, ret, live = matrix(occupancy, z, local, q=256, h=4)
        entries = (k, density * birth, ret,
                   min(density * live, quiet) if not sum(occupancy) else density * live)
        result[occupancy] = tuple(round_up(entry, bits) for entry in entries)
    return result


def _multinomial(counts):
    result = factorial(sum(counts))
    for count in counts:
        result //= factorial(count)
    return result


def odd_coefficients(z, alpha, singleton, full, *, singleton01=None, bits=60, progress=None):
    """Arb coefficients, with singleton orientations kept as separate axes."""
    alpha, singleton, full = F(alpha), F(singleton), F(full)
    singleton01 = singleton if singleton01 is None else F(singleton01)
    if not 0 < alpha <= 1 or min(singleton, singleton01, full) <= 0:
        raise ValueError('Invalid alpha or positive marker')
    entries = {occupancy: tuple(rational(x) ** rational(alpha) if x else arb(0)
                                for x in values)
               for occupancy, values in packet_table(F(z), bits).items()}
    half = []
    for a in range(17):
        for b in range(17 - a):
            for c in range(17 - a - b):
                half.append((a, b, c, _multinomial((16-a-b-c, a, b, c))))
    coeff = [[arb(0) for _ in range(33 ** 3)] for _ in range(4)]
    for index, (a, b, c, n) in enumerate(half):
        # Exchanging halves preserves the sorted four lane occupancies.
        for j in range(index, len(half)):
            d, e, f, m = half[j]
            degree = (a+d) * 33 ** 2 + (b+e) * 33 + c+f
            multiplicity = n * m * (1 if j == index else 2)
            values = entries[tuple(sorted((a+c, b+c, d+f, e+f)))]
            for entry in range(4):
                coeff[entry][degree] += multiplicity * values[entry]
        if progress and (index + 1) % 200 == 0:
            progress(dict(stage='odd-coefficients', half_profiles=index+1, total=len(half)))
    singleton_powers = [rational(singleton) ** i for i in range(33)]
    other_powers = None if singleton01 == singleton else [rational(singleton01)**i for i in range(33)]
    full_powers = [rational(full) ** i for i in range(33)]
    for a in range(33):
        for b in range(33-a):
            for c in range(33-a-b):
                degree = a * 33 ** 2 + b * 33 + c
                marker = ((singleton_powers[a+b] if other_powers is None
                           else singleton_powers[a]*other_powers[b]) * full_powers[c])
                for entry in range(4):
                    coeff[entry][degree] *= marker
    return coeff


def positive_alias_power(coefficients, target, grid=64, *, steps=128, progress=None):
    """Return four Arb upper majorants of a positive polynomial coefficient."""
    width = round(len(coefficients[0]) ** (1 / 3))
    if len(coefficients) != 4 or any(len(row) != width ** 3 for row in coefficients):
        raise ValueError('Expected four cubic coefficient arrays')
    if grid < width or grid < 2 or grid & (grid - 1):
        raise ValueError('Grid must be a power of two at least the input width')
    if steps < 1 or steps & (steps - 1) or len(target) != 3 or any(t < 0 for t in target):
        raise ValueError('Invalid target or power-of-two step count')
    if any(value.lower() < 0 for row in coefficients for value in row):
        raise ValueError('Input must enclose known nonnegative coefficients')
    positive = [sum(row, arb(0)) for row in coefficients]
    a, b, c, d = map(float, positive)
    eigenvalue = (a+d+sqrt((a-d)**2+4*b*c))/2
    v = [b, eigenvalue-a] if eigenvalue-a >= eigenvalue-d else [eigenvalue-d, c]
    if not all(isfinite(x) and x > 0 for x in v):
        raise ArithmeticError('Could not propose a positive similarity vector')
    # A proposed vector need not be an exact eigenvector. These are exact
    # dyadic witnesses, and the subsequent row-sum upper bound is directed.
    vector = [arb(x / max(v)) for x in v]
    row_sums = [(positive[2*r]*vector[0]+positive[2*r+1]*vector[1])/vector[r]
                for r in range(2)]
    scale = max(value.upper() for value in row_sums)
    factors = [vector[s] / vector[r] / scale for r in range(2) for s in range(2)]
    real_lo = np.zeros((4, grid, grid, grid))
    real_hi = np.zeros_like(real_lo)
    for entry, row in enumerate(coefficients):
        for degree, coefficient in enumerate(row):
            if not coefficient.is_zero():
                a, remainder = divmod(degree, width ** 2)
                b, c = divmod(remainder, width)
                real_lo[entry, a, b, c], real_hi[entry, a, b, c] = binary64_enclosure(coefficient*factors[entry])
    value = ((real_lo, real_hi), (np.zeros_like(real_lo), np.zeros_like(real_lo)))
    if progress:
        progress(dict(stage='interval-fft', grid=grid, steps=steps))
    with np.errstate(over='raise', invalid='raise', divide='raise', under='ignore'):
        transformed = interval_fft(value, (1, 2, 3))
        del value, real_lo, real_hi
        powered = interval_matrix_power(transformed, steps)
        del transformed
        selected = selected_inverse(powered, target)
    if np.any(selected[1][0] > 0) or np.any(selected[1][1] < 0):
        raise ArithmeticError('Imaginary coefficient enclosure excludes zero')
    region = tuple(arb(max(0., selected[0][1][entry])) * scale ** steps
                   * vector[entry // 2] / vector[entry % 2] for entry in range(4))
    return region, dict(grid=grid, steps=steps, input_width=width,
        similarity_vector=[str(x) for x in vector], scale=str(scale),
        normalized_alias_real_lower=selected[0][0].tolist(),
        normalized_alias_real_upper=selected[0][1].tolist(),
        normalized_alias_imag_lower=selected[1][0].tolist(),
        normalized_alias_imag_upper=selected[1][1].tolist(),
        arithmetic='explicit-binary64-outward-interval-radix2', directed=True)


def even_region(active_pairs, z, alpha, *, bits=60):
    """Arb two-step composition followed by the exact regional coefficient."""
    alpha = rational(alpha)
    table = packet_table(F(z), bits)
    groups = [(a+b, comb(16, a)*comb(16, b)*(1 if a == b else 2),
               tuple(rational(x) for x in table[a, a, b, b]))
              for a, b in combinations_with_replacement(range(17), 2)]
    coeff = [[arb(0) for _ in range(65)] for _ in range(4)]
    for degree, count, left in groups:
        for other_degree, other_count, right in groups:
            if degree + other_degree > active_pairs:
                continue
            for r in range(2):
                for s in range(2):
                    entry = left[2*r]*right[s]+left[2*r+1]*right[2+s]
                    if not entry.is_zero():
                        coeff[2*r+s][degree+other_degree] += count*other_count*entry**alpha
    # This rational tilt is for conditioning only; coefficient extraction
    # cancels it exactly through the normalizer below.
    marker = F(active_pairs, 4096-active_pairs)
    marker_arb = rational(marker)
    base = tuple(arb_poly([value*marker_arb**degree for degree, value in enumerate(row)]) for row in coeff)
    for _ in range(6):
        base = poly_product(base, base, active_pairs)
    normalizer = arb(comb(4096, active_pairs))*marker_arb**active_pairs
    return tuple(arb(max(arb(0), polynomial[active_pairs].upper()))/normalizer for polynomial in base)


def receipt(active_pairs=512, runs=128, z=DEFAULT_MARK, alpha=DEFAULT_ALPHA,
            singleton=DEFAULT_SINGLETON, full=DEFAULT_FULL, *, grid=64,
            precision=256, bits=60, threshold=-80, progress=None):
    z, alpha, singleton, full = map(F, (z, alpha, singleton, full))
    if not 0 < active_pairs < 4096 or not 1 <= runs <= min(active_pairs, 4096-active_pairs):
        raise ValueError('Invalid contracted-cycle support class')
    if not 0 < z < 1 or not 0 < alpha <= 1 or min(singleton, full) <= 0 or precision < 128:
        raise ValueError('Invalid rational witness or precision')
    with ctx.workprec(precision):
        if progress:
            progress(dict(stage='exact-packets-and-even-region'))
        even = even_region(active_pairs, z, alpha, bits=bits)
        coeff = odd_coefficients(z, alpha, singleton, full, bits=bits, progress=progress)
        target = (runs, runs, active_pairs-runs)
        odd, diagnostics = positive_alias_power(coeff, target, grid, progress=progress)
        normalizer = (arb(_multinomial((4096-active_pairs-runs,)+target))
                      * rational(singleton)**(2*runs)*rational(full)**(active_pairs-runs))
        odd = tuple(entry/normalizer for entry in odd)
        pair = tuple(even[2*r]*odd[s]+even[2*r+1]*odd[2+s] for r in range(2) for s in range(2))
        product = scalar_power(pair, 128)
        log2 = lambda value: value.log()/arb(2).log()
        count = circular_run_count(4096, active_pairs, runs)
        exponent = rational(alpha)*(2*active_pairs*(2048-128*log2(arb(255)))
                                    -log2(arb(255))-LIMIT*log2(rational(z)))
        bound = log2(arb(count))+exponent+log2(product[0]+product[1])
        # An exact rational ceiling is convenient for machine consumers;
        # Arb's printed decimal endpoints themselves carry display radii.
        decimal_scale = 10**12
        upper_rational = F(int((bound.upper()*decimal_scale).upper().ceil().fmpz()), decimal_scale)
        return dict(construction='shifted-pairs-alternating-h4',
            active_pairs=active_pairs, active_groups=2*active_pairs, contracted_runs=runs,
            mark=str(z), alpha=str(alpha), singleton_marker=str(singleton), full_marker=str(full),
            precision=precision, packet_mantissa_bits=bits, even_composition_steps=2,
            class_support_count=str(count), forbidden_weight=LIMIT, field_order=256,
            log2_upper_enclosure=str(bound), log2_upper_endpoint=str(bound.upper()),
            log2_upper_rational=str(upper_rational),
            accepted=bool(bound < threshold), threshold=threshold,
            odd_entry_log2_upper=[str(log2(entry)) for entry in odd],
            even_entry_log2_upper=[str(log2(entry)) for entry in even],
            alias_diagnostics=diagnostics, full_under_even_matching_only=True,
            directed=True, whole_code_certificate=False)


if __name__ == '__main__':
    import argparse
    import json
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--active-pairs', type=int, default=512)
    parser.add_argument('--runs', type=int, default=128)
    parser.add_argument('--mark', type=F, default=DEFAULT_MARK)
    parser.add_argument('--alpha', type=F, default=DEFAULT_ALPHA)
    parser.add_argument('--singleton', type=F, default=DEFAULT_SINGLETON)
    parser.add_argument('--full', type=F, default=DEFAULT_FULL)
    parser.add_argument('--grid', type=int, default=64)
    parser.add_argument('--precision', type=int, default=256)
    parser.add_argument('--bits', type=int, default=60)
    parser.add_argument('--threshold', type=int, default=-80)
    args = parser.parse_args()
    print(json.dumps(receipt(args.active_pairs, args.runs, args.mark, args.alpha,
        args.singleton, args.full, grid=args.grid, precision=args.precision,
        bits=args.bits, threshold=args.threshold,
        progress=lambda item: print(json.dumps(item), flush=True))), flush=True)
