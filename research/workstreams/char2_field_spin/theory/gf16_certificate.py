"""Directed upper bounds for GF16 GRS[16,8], t16, uniform-GL4 SPIN.

This file fixes A_j=(1,j,j^2,j^3), C=A^T, independent GF16* input
diagonals, and independent uniform GL(4,16) pre-output refreshes.
All certification arithmetic is rational or directed Decimal arithmetic.
Floating searches only propose rational tilts for independent replay.
"""

import argparse
from decimal import Decimal, localcontext
from fractions import Fraction
import json
import math

import numpy as np

from q1_batched_screen import fixed_support_multiplicity
from state_screen import birth_zero_count, zero_patterns
from subfield_certificate import Interval, HIGH, LOW, ZERO, ONE, integer_margin


def interval_local_envelope(z):
    """Intervals enclosing the positive local upper-envelope entries.

    The singleton joint return is exact for the fixed monomial A.
    Other returns use the smaller of the MDS and invertible-map bounds.
    """
    q, width, state = 16, 16, 4
    if z.lo <= 0 or z.hi > 1:
        raise ValueError('require 0 < z <= 1')
    f = Interval.rational(1, 15)+Interval.rational(14, 15)*z
    zs = [z**i for i in range(17)]
    fs = [f**i for i in range(17)]
    complements = [(1-z)**i for i in range(17)]
    patterns = zero_patterns(width, state, q)
    denominator = q**state-1
    matrices = []
    for w in range(17):
        p0 = Interval.rational(birth_zero_count(w, state, q), (q-1)**w)
        g = ZERO
        for zeros, count in enumerate(patterns):
            for active_zeros in range(max(0, zeros-(width-w)), min(w, zeros)+1):
                inactive_zeros = zeros-active_zeros
                coefficient = count*math.comb(w, active_zeros)*math.comb(width-w, inactive_zeros)
                g += (Interval.rational(coefficient, denominator)
                      *zs[width-w-inactive_zeros+active_zeros]*fs[w-active_zeros])
        ret = (1-p0)/denominator
        jlo = ret*zs[width]
        jhi = ret*zs[max(0, width-state+1-w)]
        if w == 1:
            jlo = jhi = (13*zs[15]+zs[16]+2*zs[14])/(16*denominator)
        elif w >= 2:
            invertible = sum((Interval.rational(math.comb(width, h),
                                                denominator*(q-1)**max(0, h-width+w))
                              *zs[width-h]*complements[h] for h in range(17)), ZERO)
            jhi = Interval(min(jhi.lo, invertible.lo), min(jhi.hi, invertible.hi))
        matrices.append([[zs[w]*p0, zs[w]*(1-p0)], [jhi, g-jlo]])
    return matrices


def upper_polynomial_product(left, right, limit):
    """Positive 2x2 coefficients, every arithmetic operation rounded up."""
    size = min(limit+1, len(left)+len(right)-1)
    zero = Decimal(0)
    out = [[zero, zero, zero, zero] for _ in range(size)]
    add, mul = HIGH.add, HIGH.multiply
    for i, a in enumerate(left):
        for j, b in enumerate(right[:size-i]):
            c = out[i+j]
            c[0] = add(c[0], add(mul(a[0], b[0]), mul(a[1], b[2])))
            c[1] = add(c[1], add(mul(a[0], b[1]), mul(a[1], b[3])))
            c[2] = add(c[2], add(mul(a[2], b[0]), mul(a[3], b[2])))
            c[3] = add(c[3], add(mul(a[2], b[1]), mul(a[3], b[3])))
    return out


def upper_regions(z, maximum_h):
    local = interval_local_envelope(z)
    factors = [[HIGH.multiply(entry.hi, Decimal(math.comb(16, w)))
                for row in matrix for entry in row]
               for w, matrix in enumerate(local[:maximum_h+1])]
    out = [[Decimal(1), Decimal(0), Decimal(0), Decimal(1)]]
    power = 512
    while power:
        if power & 1:
            out = upper_polynomial_product(out, factors, maximum_h)
        power >>= 1
        if power:
            factors = upper_polynomial_product(factors, factors, maximum_h)
    return [[HIGH.divide(entry, Decimal(math.comb(8192, h))) for entry in matrix]
            for h, matrix in enumerate(out)]


def upper_exact_source(regions, active_rows):
    """Positive support-mask dynamic program, with exact MDS multiplicities."""
    if not 1 <= active_rows <= 4:
        raise ValueError('exact source computation is limited to four rows')
    matrices = np.array(regions, dtype=object).reshape((-1, 2, 2))
    masks = [tuple((mask >> j) & 1 for j in range(active_rows)) for mask in range(1 << active_rows)]
    with localcontext(HIGH):
        current = np.full((1,)*active_rows+(2,), Decimal(0), dtype=object)
        current[(0,)*active_rows+(0,)] = Decimal(1)
        for region in range(16):
            next_state = np.full((region+2,)*active_rows+(2,), Decimal(0), dtype=object)
            for mask in masks:
                target = tuple(slice(bit, region+1+bit) for bit in mask)+(slice(None),)
                next_state[target] += current @ matrices[sum(mask)]
            current = next_state
        counts = [Decimal(fixed_support_multiplicity(16, 8, w, 16)) if w >= 9 else Decimal(0)
                  for w in range(17)]
        for axis in range(active_rows):
            shape = [1]*(active_rows+1)
            shape[axis] = 17
            current *= np.array(counts, dtype=object).reshape(shape)
        return np.sum(current)


def upper_coupled(regions, active_rows):
    add, mul = HIGH.add, HIGH.multiply
    matrix = [Decimal(0)]*4
    for h, regional in enumerate(regions):
        weight = Decimal(math.comb(active_rows, h)*16**h)
        matrix = [add(a, mul(weight, b)) for a, b in zip(matrix, regional)]
    vector = [Decimal(1), Decimal(0)]
    for _ in range(16):
        vector = [add(mul(vector[0], matrix[0]), mul(vector[1], matrix[2])),
                  add(mul(vector[0], matrix[1]), mul(vector[1], matrix[3]))]
    return HIGH.divide(add(*vector), Decimal(16**(8*active_rows)))


def class_upper(active_rows, numerator, denominator, exact_source=None):
    if not 1 <= active_rows < 8192:
        raise ValueError('unsupported active-row count')
    if exact_source is None:
        exact_source = active_rows <= 4
    z = Interval.rational(numerator, denominator)
    regions = upper_regions(z, active_rows)
    moment = upper_exact_source(regions, active_rows) if exact_source else upper_coupled(regions, active_rows)
    numerator_upper = HIGH.multiply(moment, Decimal(math.comb(8192, active_rows)))
    denominator_lower = LOW.multiply(Decimal(15), (z**13107).lo)
    bound = HIGH.divide(numerator_upper, denominator_lower)
    return dict(active_rows=active_rows, z=[numerator, denominator],
                method='exact-source' if exact_source else 'coupled-source-envelope',
                failure_upper=str(bound), integer_bits=integer_margin(bound))


CONFIG = dict(q=16, k=65536, n=131072, outer=[16, 8], width=16, state=4,
              bad_weight_threshold=13107, refresh='uniform-gl4-pre-output',
              diagonal='independent-uniform-nonzero', maps='monomial-vandermonde-transpose')


# Rational witnesses only. Their selection used floating searches; all claimed
# bounds are recomputed below with directed arithmetic and no SciPy dependency.
# Entry a-1 is the numerator of z for active-row class a; denominator is 2^32.
SPARSE_TILT_NUMERATORS = (
    4290114385, 4288363716, 4285983836, 4283082122, 4284857262, 4281376731,
    4277645861, 4273809296, 4269925272, 4266036137, 4262174043, 4258362174,
    4254615470, 4250941498, 4247341609, 4243812231, 4240346292, 4236934787,
    4233568168, 4230237417, 4226934668, 4223653438, 4220388600, 4217136226,
    4213893368, 4210657854, 4207428100, 4204202962, 4200981618, 4197763482,
    4194548136, 4191335286, 4188124724, 4184916302, 4181709919, 4178505502,
    4175303004, 4172102391, 4168903641, 4165706740, 4162511678, 4159318453,
    4156127062, 4152937506, 4149749787, 4146563908, 4143379874, 4140197690,
    4137017359, 4133838889, 4130662285, 4127487554, 4124314700, 4121143733,
    4117974657, 4114807479, 4111642208, 4108478850, 4105317411, 4102157900,
    4099000324, 4095844690, 4092691006, 4089539279, 4086389518, 4083241729,
    4080095920, 4076952101, 4073810277, 4070670458, 4067532651, 4064396864,
    4061263106, 4058131383, 4055001706, 4051874081, 4048748517, 4045625021,
    4042503603, 4039384270, 4036267031, 4033151894, 4030038866, 4026927957,
    4023819175, 4020712527, 4017608023, 4014505670, 4011405476, 4008307451,
    4005211602, 4002117938, 3999026466, 3995937196, 3992850136, 3989765293,
    3986682676, 3983602293, 3980524153, 3977448264, 3974374634, 3971303272,
    3968234185, 3965167382, 3962102871, 3959040660, 3955980758, 3952923173,
    3949867912, 3946814985, 3943764399, 3940716163, 3937670284, 3934626772,
    3931585633, 3928546876, 3925510509, 3922476541, 3919444979, 3916415832,
    3913389107, 3910364813, 3907342957, 3904323548, 3901306594, 3898292103,
    3895280083,
)


def certify_sparse(maximum_rows=127):
    if not 1 <= maximum_rows <= len(SPARSE_TILT_NUMERATORS):
        raise ValueError('fixed sparse witnesses cover exactly rows 1 through 127')
    entries = []
    total = Decimal(0)
    for active in range(1, maximum_rows+1):
        denominator = 2**32
        numerator = SPARSE_TILT_NUMERATORS[active-1]
        entry = class_upper(active, numerator, denominator)
        entries.append(entry)
        total = HIGH.add(total, Decimal(entry['failure_upper']))
    return dict(format='gf16-gl4-sparse-v1', complete=False,
                config=CONFIG,
                covered_active_rows=[1, maximum_rows], uncovered_active_rows=[maximum_rows+1, 8192],
                entries=entries, failure_upper=str(total), integer_bits=integer_margin(total))


def verify_sparse(record):
    """Replay supplied rational witnesses without invoking floating searches."""
    if record.get('format') != 'gf16-gl4-sparse-v1' or record.get('complete') is not False:
        raise ValueError('incorrect format or claim of whole-code coverage')
    if record.get('config') != CONFIG:
        raise ValueError('incorrect construction parameters')
    entries = record.get('entries', [])
    maximum = len(entries)
    if not 1 <= maximum < 8192:
        raise ValueError('empty or unsupported sparse coverage')
    if record.get('covered_active_rows') != [1, maximum] or record.get('uncovered_active_rows') != [maximum+1, 8192]:
        raise ValueError('incorrect coverage boundaries')
    total = Fraction(0)
    for active, entry in enumerate(entries, 1):
        method = 'exact-source' if active <= 4 else 'coupled-source-envelope'
        if entry.get('active_rows') != active or entry.get('method') != method:
            raise ValueError('gap, overlap, or wrong source bound')
        numerator, denominator = entry['z']
        recomputed = class_upper(active, numerator, denominator)
        upper = Decimal(entry['failure_upper'])
        if not upper.is_finite() or upper < Decimal(recomputed['failure_upper']):
            raise ValueError('understated class bound')
        if entry['integer_bits'] > integer_margin(upper):
            raise ValueError('overstated class margin')
        total += Fraction(upper)
    upper = Decimal(record['failure_upper'])
    if not upper.is_finite() or Fraction(upper) < total:
        raise ValueError('understated total bound')
    if record['integer_bits'] > integer_margin(upper):
        raise ValueError('overstated total margin')
    return True


def certify_full():
    from coupled_chernoff import GF16_K16_MIDDLE_WITNESSES, interval_cover
    from subfield_bound_audit import exact_bulk_tail_check
    sparse = certify_sparse(127)
    middle, details = interval_cover(GF16_K16_MIDDLE_WITNESSES, interval_local_envelope)
    dense = exact_bulk_tail_check(first=5998, z=Fraction(1, 52), q=16,
                                  distance=Fraction(1, 10), margin_bits=64)
    if not dense['rigorous']:
        raise ValueError('the exact dense-tail comparison failed')
    dense_upper = Interval.rational(1, 1 << 64).hi
    total = HIGH.add(HIGH.add(Decimal(sparse['failure_upper']), middle.hi), dense_upper)
    return dict(format='gf16-gl4-full-v1', complete=True, config=CONFIG,
                covered_active_rows=[1, 8192], sparse=sparse,
                middle=dict(covered_active_rows=[128, 5997], entries=details,
                            failure_upper=str(middle.hi)),
                dense=dense, failure_upper=str(total), integer_bits=integer_margin(total))


def verify_full(record):
    from coupled_chernoff import interval_cover
    from subfield_bound_audit import exact_bulk_tail_check
    if (record.get('format') != 'gf16-gl4-full-v1' or record.get('complete') is not True
            or record.get('config') != CONFIG or record.get('covered_active_rows') != [1, 8192]):
        raise ValueError('incorrect full certificate scope')
    sparse = record['sparse']
    if sparse.get('covered_active_rows') != [1, 127]:
        raise ValueError('the sparse block does not cover exactly 1 through 127')
    verify_sparse(sparse)
    middle = record['middle']
    if middle.get('covered_active_rows') != [128, 5997]:
        raise ValueError('the middle block has incorrect boundaries')
    witnesses = [{key: entry[key] for key in ['first', 'last', 'x', 'z']} for entry in middle['entries']]
    recalculated, details = interval_cover(witnesses, interval_local_envelope)
    claimed_middle_sum = Fraction(0)
    for supplied, actual in zip(middle['entries'], details):
        contribution = Decimal(supplied['contribution_upper'])
        if not contribution.is_finite() or contribution < Decimal(actual['contribution_upper']):
            raise ValueError('understated middle-interval contribution')
        claimed_middle_sum += Fraction(contribution)
    middle_upper = Decimal(middle['failure_upper'])
    if (not middle_upper.is_finite() or middle_upper < recalculated.hi
            or Fraction(middle_upper) < claimed_middle_sum):
        raise ValueError('understated middle contribution')
    expected_dense = exact_bulk_tail_check(first=5998, z=Fraction(1, 52), q=16,
                                           distance=Fraction(1, 10), margin_bits=64)
    if not expected_dense['rigorous'] or record['dense'] != expected_dense:
        raise ValueError('incorrect dense-tail comparison')
    exact_sum = Fraction(Decimal(sparse['failure_upper']))+Fraction(middle_upper)+Fraction(1, 1 << 64)
    upper = Decimal(record['failure_upper'])
    if not upper.is_finite() or Fraction(upper) < exact_sum:
        raise ValueError('understated full contribution')
    if record['integer_bits'] > integer_margin(upper):
        raise ValueError('overstated full margin')
    return True


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--max-rows', type=int, default=127)
    parser.add_argument('--verify')
    parser.add_argument('--certify-full', action='store_true')
    args = parser.parse_args()
    if args.verify:
        with open(args.verify, encoding='utf-8-sig') as source:
            record = json.load(source)
        verifier = verify_full if record.get('complete') else verify_sparse
        verifier(record)
        print('Verified every claimed class and its rational witness.')
    elif args.certify_full:
        print(json.dumps(certify_full(), indent=2))
    else:
        print(json.dumps(certify_sparse(args.max_rows), indent=2))
