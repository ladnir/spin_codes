"""Positive actual-state moment with exact forced-birth degrees zero through four.

For one quartet, P=P0+[(1+p*s*nu*Y)^4-1], where P0 has empty-old row1
and all old-live rows of the accepted-MDS polynomial.  With g=width/4,

 K_f(u,v)=sum_b C(g,b)c_(b,f)(p*s)^f
           [X^u Y^(v-f)]P0^(g-b) / C(width,u),
 c_(b,f)=[w^f]((1+w)^4-1)^b.

One positive convolution pass saves only powers g-fmax,...,g.  Regional
polynomial squaring is truncated at fmax<=4, then each coefficient is
quota-conditioned and minimized across p in the common population basis.
Write the resulting region as Z+sum_f B_f nu^f.  B0 has row0 EXACTLY zero;
positive-degree (0,0) entries contain real regional excursions and remain.
Global propagation uses polynomial ROW VECTORS, not large polynomial
matrix powers.  Its degree-zero vector is exactly e0 and moment exactly1.

log_t is the actual FREE death marker.  No nu evaluation is performed:
outputs are bare [nu^f] coefficients.  Optional positive floors enlarge
only supported entries and never populate structural zeros.  All numerical
moments remain FLOAT diagnostics, not directed or distance certificates.
The default profile report covers only D in [32769,196607].  Other D
classes require their separately proved bounds.  No optimizer, subprocess,
raw data files, or new dependencies.  Numerical jobs must be serial.
"""

import argparse
from fractions import Fraction
from itertools import combinations, product
from math import comb
from operator import index
import os
import shutil
import threading
import time

import numpy as np


BC_LD = np.longdouble
BC_PROFILE = (4096, 786432, 786432, 4, 1048576)
BC_CENTRAL = ((32769, 196607),)


def _bc_dependency(name):
    if name in globals():
        return globals()[name]
    import quartet_actual_mds_moment as module
    return getattr(module, name)


def _bc_floor(bits):
    if bits is None:
        return None
    bits = index(bits)
    if not 1 <= bits <= 4000:
        raise ValueError('Floor precision must be in 1..4000 bits')
    return bits


def _bc_counts(b, fmax):
    value = [1] + [0] * fmax
    for _ in range(b):
        after = [0] * (fmax + 1)
        for f in range(fmax + 1):
            for j in range(1, min(4, f) + 1):
                after[f] += comb(4, j) * value[f - j]
        value = after
    return value


def _bc_positive_range(array, stage):
    if np.any(array < 0) or not np.all(np.isfinite(array)):
        raise ArithmeticError('Invalid nonnegative data at ' + stage)


def _bc_normalize(array, scale, floor_bits=None, protect_constant=False):
    """Normalize supported positive entries; structural zeros stay zero."""
    _bc_positive_range(array, 'normalization input')
    supported = array > 0
    if not np.any(supported):
        return None
    norm = array.max()
    scale = BC_LD(scale) + np.log2(norm)
    array = array / norm
    if floor_bits is not None:
        array[supported] = np.maximum(array[supported], np.exp2(-BC_LD(floor_bits)))
    elif np.any(array[supported] <= 0):
        raise ArithmeticError('Positive normalization underflow; enable the explicit floor majorant')
    if protect_constant:
        constant = np.exp2(-scale)
        if constant <= 0:
            raise ArithmeticError('The exact one-step constant is outside representable range')
        array[0, 0] = constant
    return array, scale


def _bc_from_logs(logs, floor_bits=None):
    if np.any(np.isnan(logs)) or np.any(np.isposinf(logs)):
        raise ArithmeticError('Invalid coefficient log sum')
    supported = np.isfinite(logs)
    if not np.any(supported):
        return None
    scale = logs[supported].max()
    result = np.zeros_like(logs)
    differences = logs[supported] - scale
    if floor_bits is not None:
        differences = np.maximum(differences, -BC_LD(floor_bits))
    result[supported] = np.exp2(differences)
    if np.any(result[supported] <= 0):
        raise ArithmeticError('Mixed coefficient scales underflow; enable the explicit floor majorant')
    return result, scale


def _bc_accumulate(logs, array, scale):
    _bc_positive_range(array, 'positive coefficient sum')
    positive = array > 0
    term = np.full_like(array, BC_LD('-inf'))
    np.log2(array, out=term, where=positive)
    term += scale
    np.logaddexp2(logs, term, out=logs)


def _bc_product(left, right):
    """Zero-preserving product; reject disappearing positive scalar products."""
    if left is None or right is None:
        return None
    a, sa = left
    b, sb = right
    apos, bpos = a[a > 0], b[b > 0]
    if not len(apos) or not len(bpos):
        return None
    if apos.min() * bpos.min() <= 0:
        raise ArithmeticError('Positive coefficient product underflow; enable the floor majorant')
    value = a @ b
    _bc_positive_range(value, 'coefficient matrix/vector product')
    return value, sa + sb


def actual_birth_one_step(width=512, fmax=4, log_s=1, log_t=-1,
                          log_z=0, log_eta=-2, p=1/3,
                          floor_bits=None, progress=False):
    """Return B[0..fmax], with K(nu)=Z+sum_f B_f nu^f.

    Each nonzero block is (normalized matrix, log2 scale).  None denotes
    an identically zero coefficient.  B0's entire source row0 is zero.
    """
    width, fmax = index(width), index(fmax)
    if width < 4 or width % 4 or not 0 <= fmax <= 4:
        raise ValueError('Require quartet width and fmax in 0..4')
    floor_bits = _bc_floor(floor_bits)
    groups = width // 4
    p, log_s = BC_LD(p), BC_LD(log_s)
    local, local_scale = _bc_dependency('actual_mds_local_scaled')(
        log_s, log_t, log_z, 0, log_eta, p)
    local = local.copy()
    local[0, 1:] = 0
    local, local_scale = _bc_normalize(local, local_scale, floor_bits, protect_constant=True)
    needed = set(range(max(0, groups - fmax), groups + 1))
    snapshots = {0: (np.ones((1, 1), dtype=BC_LD), BC_LD(0))} if 0 in needed else {}
    shape = width + 1, width + 1
    current = np.zeros(shape, dtype=BC_LD)
    after = np.empty(shape, dtype=BC_LD)
    scratch = np.empty(shape, dtype=BC_LD)
    current[0, 0] = 1
    scale = BC_LD(0)
    started = time.monotonic()
    shifts = [(r, v, local[r, v]) for r in range(5) for v in range(5) if local[r, v] > 0]
    local_min = min(weight for _, _, weight in shifts)
    for group in range(groups):
        size, target_size = 4 * group + 1, 4 * group + 5
        source = current[:size, :size]
        target = after[:target_size, :target_size]
        target.fill(0)
        if source[source > 0].min() * local_min <= 0:
            raise ArithmeticError('Snapshot convolution would lose a positive product')
        work = scratch[:size, :size]
        for r, v, weight in shifts:
            np.multiply(source, weight, out=work)
            np.add(after[r:r + size, v:v + size], work, out=after[r:r + size, v:v + size])
        u = np.arange(target_size)[:, None]
        v = np.arange(target_size)[None, :]
        support = v <= 4 * u
        if np.any(target[~support] != 0) or np.any(target[support] <= 0):
            raise ArithmeticError('P0 power violated its exact support or lost a positive coefficient')
        normalized, scale = _bc_normalize(target, scale + local_scale, floor_bits, protect_constant=True)
        target[:] = normalized
        current, after = after, current
        power = group + 1
        if power in needed:
            snapshots[power] = (current[:target_size, :target_size].copy(), scale)
        if progress and (power % 32 == 0 or power == groups):
            print(dict(stage='birth_P0_convolution', quartets=power,
                       seconds=round(time.monotonic() - started, 3)), flush=True)
    log_choose = np.array([np.log2(BC_LD(comb(width, u))) for u in range(width + 1)])
    counts = [_bc_counts(b, fmax) for b in range(min(groups, fmax) + 1)]
    coefficients = []
    for f in range(fmax + 1):
        logs = np.full(shape, BC_LD('-inf'))
        for b in range(min(groups, f) + 1):
            if not counts[b][f]:
                continue
            polynomial, polynomial_scale = snapshots[groups - b]
            length = len(polynomial)
            term = np.full(polynomial.shape, BC_LD('-inf'))
            np.log2(polynomial, out=term, where=polynomial > 0)
            term += (polynomial_scale + np.log2(BC_LD(comb(groups, b) * counts[b][f]))
                     + f * (np.log2(p) + log_s))
            np.logaddexp2(logs[:length, f:f + length], term, out=logs[:length, f:f + length])
        logs -= log_choose[:, None]
        if f == 0:
            logs[0, :] = BC_LD('-inf')  # The exact constant is represented by Z.
        else:
            # For old support empty, exact f births simply choose f slots.
            logs[0, :] = BC_LD('-inf')
            if f <= width:
                logs[0, f] = np.log2(BC_LD(comb(width, f))) + f * (np.log2(p) + log_s)
        block = _bc_from_logs(logs, floor_bits)
        if block is not None and f == 0 and np.any(block[0][0] != 0):
            raise ArithmeticError('Degree zero acquired a nonempty path from empty state')
        coefficients.append(block)
    return coefficients


def actual_birth_polynomial_square(coefficients, floor_bits=None):
    """Positive (Z+B(nu))^2-Z, retaining every supported birth degree."""
    fmax = len(coefficients) - 1
    shape = next(block[0].shape for block in coefficients if block is not None)
    result = []
    for f in range(fmax + 1):
        logs = np.full(shape, BC_LD('-inf'))
        for a in range(f + 1):
            term = _bc_product(coefficients[a], coefficients[f - a])
            if term is not None:
                _bc_accumulate(logs, *term)
        if coefficients[f] is not None:
            matrix, scale = coefficients[f]
            edges = np.zeros(shape, dtype=BC_LD)
            edges[0, :] += matrix[0, :]
            edges[:, 0] += matrix[:, 0]
            _bc_accumulate(logs, edges, scale)
        block = _bc_from_logs(logs, floor_bits)
        if f == 0 and block is not None and np.any(block[0][0] != 0):
            raise ArithmeticError('Regional degree zero lost its exact zero source row')
        result.append(block)
    return result


def actual_birth_regional(width=512, region_updates=16, eligible=4096,
                          fmax=4, p_values=(1/3,), log_s=1, log_t=-1,
                          log_z=0, log_eta=-2, floor_bits=None, progress=False):
    """Coefficientwise common-basis minimum AFTER whole-region quota conditioning."""
    width, region_updates, eligible, fmax = map(index, (width, region_updates, eligible, fmax))
    if region_updates < 1 or region_updates & (region_updates - 1):
        raise ValueError('Region update count must be a positive power of two')
    slots = width * region_updates
    if not 0 <= eligible <= slots:
        raise ValueError('Eligible quota is outside its support')
    floor_bits = _bc_floor(floor_bits)
    values = tuple(BC_LD(p) for p in p_values)
    if not values or any(not 0 < p < 1 for p in values):
        raise ValueError('Quota tilts must lie in (0,1)')
    common = None
    details = []
    for p in values:
        started = time.monotonic()
        coefficients = actual_birth_one_step(width, fmax, log_s, log_t, log_z, log_eta, p, floor_bits, progress)
        for step in range(region_updates.bit_length() - 1):
            coefficients = actual_birth_polynomial_square(coefficients, floor_bits)
            if progress:
                print(dict(stage='birth_region_square', p=float(p), updates=1 << (step + 1),
                           seconds=round(time.monotonic() - started, 3)), flush=True)
        quota = (np.log2(BC_LD(comb(slots, eligible))) + eligible * np.log2(p)
                 + (slots - eligible) * np.log2(1 - p))
        conditioned = [None if block is None else (block[0], block[1] - quota) for block in coefficients]
        if common is None:
            common = conditioned
        else:
            for f in range(fmax + 1):
                if common[f] is None or conditioned[f] is None:
                    common[f] = None
                    continue
                a, sa = common[f]
                b, sb = conditioned[f]
                la, lb = np.full_like(a, BC_LD('-inf')), np.full_like(b, BC_LD('-inf'))
                np.log2(a, out=la, where=a > 0)
                np.log2(b, out=lb, where=b > 0)
                common[f] = _bc_from_logs(np.minimum(la + sa, lb + sb), floor_bits)
        details.append(dict(p=float(p), log2_quota=float(quota), seconds=round(time.monotonic() - started, 3)))
    if common[0] is not None and np.any(common[0][0][0] != 0):
        raise ArithmeticError('Quota conditioning charged the exact empty source row')
    return common, details


def actual_birth_global_vectors(coefficients, regions=256, floor_bits=None, progress=False):
    """Forward polynomial ROW DP; v0=e0 exactly, at most 14 products per region.

    The e0*B_f terms are row lookups, reducing the actual vector-matrix
    products to ten per region when fmax=4.  No polynomial matrix powers.
    """
    regions = index(regions)
    if regions < 0:
        raise ValueError('Region count must be nonnegative')
    fmax = len(coefficients) - 1
    width1 = next(block[0].shape[0] for block in coefficients if block is not None)
    if coefficients[0] is not None and np.any(coefficients[0][0][0] != 0):
        raise ValueError('Degree-zero positive part must have zero source row0')
    vectors = [None] * (fmax + 1)
    empty = np.zeros(width1, dtype=BC_LD)
    empty[0] = 1
    vectors[0] = empty, BC_LD(0)
    for region in range(regions):
        after = [vectors[0]]  # This is an identity, not a rounded matrix product.
        for f in range(1, fmax + 1):
            logs = np.full(width1, BC_LD('-inf'))
            if vectors[f] is not None and vectors[f][0][0] > 0:
                logs[0] = vectors[f][1] + np.log2(vectors[f][0][0])  # v_f Z
            for b in range(f + 1):
                if coefficients[b] is None or vectors[f - b] is None:
                    continue
                if f == b:
                    matrix, scale = coefficients[b]
                    _bc_accumulate(logs, matrix[0], scale)
                else:
                    term = _bc_product(vectors[f - b], coefficients[b])
                    if term is not None:
                        _bc_accumulate(logs, *term)
            after.append(_bc_from_logs(logs, floor_bits))
        vectors = after
        if progress and ((region + 1) % 64 == 0 or region + 1 == regions):
            print(dict(stage='birth_global_vector_DP', regions=region + 1), flush=True)
    moments = [0.0]
    for f in range(1, fmax + 1):
        moments.append(float('-inf') if vectors[f] is None else
                       float(vectors[f][1] + np.log2(vectors[f][0].sum())))
    if not (vectors[0][1] == 0 and vectors[0][0][0] == 1 and np.count_nonzero(vectors[0][0]) == 1):
        raise ArithmeticError('The exact degree-zero global control failed')
    return moments


def actual_birth_probe(profile=BC_PROFILE, fmax=4, exact_f=4, sum_births=False,
                       intervals=BC_CENTRAL, p_values=(1/3,), log_s=1, log_t=-1,
                       log_z=0, log_eta=-2, exponent=16, theta=1/16,
                       floor_bits=None, progress=False):
    """Single bounded evaluation; exact f or explicit sum of f=1..fmax."""
    h, R, Ls, _, F = map(index, profile)
    fmax, exact_f = index(fmax), index(exact_f)
    if not 1 <= exact_f <= fmax <= 4:
        raise ValueError('Require 1<=exact_f<=fmax<=4')
    coefficients, tilts = actual_birth_regional(512, 16, h, fmax, p_values,
        log_s, log_t, log_z, log_eta, floor_bits, progress)
    moments = actual_birth_global_vectors(coefficients, 256, floor_bits, progress)
    if 'actual_death_profile_bound' in globals():
        evaluate = globals()['actual_death_profile_bound']
    else:
        from quartet_actual_death_profiles import actual_death_profile_bound as evaluate
    selected = range(1, fmax + 1) if sum_births else (exact_f,)
    bounds = []
    for f in selected:
        if (F - f) % 4:
            bounds.append(dict(factual=f, geometrically_empty=True, reason='F-factual is not divisible by four'))
            continue
        if not np.isfinite(moments[f]):
            bounds.append(dict(factual=f, geometrically_empty=True, reason='zero exact-birth coefficient'))
            continue
        bound = evaluate(moments[f], profile=(h, R, Ls, f, F), intervals=intervals,
            log_death_marker=log_t, exponent=exponent, theta=theta,
            log_s=log_s, log_z=log_z, log_nu=0, log_eta=log_eta,
            birth_extracted=True)
        bounds.append(dict(factual=f, **bound))
    finite = [BC_LD(item['log2_sum_bound']) for item in bounds
              if item.get('log2_sum_bound') is not None and np.isfinite(item['log2_sum_bound'])]
    combined = float('-inf')
    if finite:
        maximum = max(finite)
        combined = float(maximum + np.log2(sum(np.exp2(value - maximum) for value in finite)))
    selected_logs = [BC_LD(moments[f]) for f in selected if np.isfinite(moments[f])]
    raw_sum = float('-inf')
    if selected_logs:
        maximum = max(selected_logs)
        raw_sum = float(maximum + np.log2(sum(np.exp2(value - maximum) for value in selected_logs)))
    return dict(profile=(h, R, Ls, exact_f, F), birth_mode='sum_1_to_fmax' if sum_births else 'exact_f',
                fmax=fmax, exact_f=None if sum_births else exact_f,
                coefficient_log2_moments=moments, selected_geometry_log2_moment=raw_sum,
                degree_zero_control_exactly_one=True, regional_tilts=tilts,
                fixed_profile_death_bounds=bounds, log2_selected_profile_bound=combined,
                log_s=float(log_s), log_t=float(log_t), log_z=float(log_z), log_eta=float(log_eta),
                birth_marker_evaluated=False, positive_floor_bits=floor_bits,
                floating_diagnostic_only=True, full_code_certificate=False)


def _bc_exact_one(width, s, t, z, eta, p, fmax):
    local = _bc_dependency('actual_mds_local_exact')(s, t, z, 0, eta, p)
    groups = width // 4
    powers = [[[Fraction(1)]]]
    for _ in range(groups):
        current = powers[-1]
        size = len(current)
        after = [[Fraction(0) for _ in range(size + 4)] for _ in range(size + 4)]
        for u, v in product(range(size), repeat=2):
            for r, w in product(range(5), repeat=2):
                after[u + r][v + w] += current[u][v] * local[r][w]
        powers.append(after)
    result = [[[Fraction(0) for _ in range(width + 1)] for _ in range(width + 1)] for _ in range(fmax + 1)]
    for f in range(fmax + 1):
        for b in range(min(f, groups) + 1):
            multiplier = comb(groups, b) * _bc_counts(b, fmax)[f] * (p * s) ** f
            if not multiplier:
                continue
            polynomial = powers[groups - b]
            for u, v in product(range(len(polynomial)), repeat=2):
                result[f][u][v + f] += multiplier * polynomial[u][v] / comb(width, u)
    result[0][0][0] = 0
    return result


def _bc_exact_square(blocks):
    fmax, size = len(blocks) - 1, len(blocks[0])
    after = [[[Fraction(0) for _ in range(size)] for _ in range(size)] for _ in blocks]
    for f in range(fmax + 1):
        for u, v in product(range(size), repeat=2):
            after[f][u][v] = sum(blocks[a][u][k] * blocks[f - a][k][v]
                                  for a in range(f + 1) for k in range(size))
            after[f][u][v] += (blocks[f][u][0] if v == 0 else 0) + (blocks[f][0][v] if u == 0 else 0)
    return after


def actual_birth_checks():
    F = Fraction
    s, t, z, eta, p = F(3, 4), F(2, 3), F(5, 4), F(7, 8), F(1, 3)
    mask_weights, _ = _bc_dependency('_am_physical_local')(s, t, z, 1, eta, p)
    physical_entries = regional_entries = 0
    for width in (4, 8):
        actual = _bc_exact_one(width, s, t, z, eta, p, 4)
        reference = [[[Fraction(0) for _ in range(width + 1)] for _ in range(width + 1)] for _ in range(5)]
        partitions = [(15,)] if width == 4 else [(first, 255 ^ first)
            for rest in combinations(range(1, 8), 3) for first in [1 + sum(1 << j for j in rest)]]
        for partition in partitions:
            maps = []
            for quartet in partition:
                positions = [j for j in range(width) if quartet >> j & 1]
                maps.append([sum(((mask >> j) & 1) << k for k, j in enumerate(positions)) for mask in range(1 << width)])
            for u in range(width + 1):
                old = (1 << u) - 1
                for new in range(1 << width):
                    births, weight = 0, Fraction(1)
                    for mapping in maps:
                        o, n = mapping[old], mapping[new]
                        births += n.bit_count() if o == 0 else 0
                        weight *= mask_weights[o][n]
                    if births <= 4:
                        reference[births][u][new.bit_count()] += weight / len(partitions)
        reference[0][0][0] = 0
        assert actual == reference
        physical_entries += 5 * (width + 1) ** 2
    one = _bc_exact_one(4, s, t, z, eta, p, 4)
    region = _bc_exact_square(one)
    for f in range(5):
        direct = sum(mask_weights[0][first] * mask_weights[first][last]
            for first, last in product(range(16), repeat=2)
            if first.bit_count() + (last.bit_count() if first == 0 else 0) == f)
        assert sum(region[f][0]) + int(f == 0) == direct
        regional_entries += 1
    assert region[1][0][0] > 0 and all(value == 0 for value in region[0][0])
    # Birth and eligible counters together; conditioning removes p exactly.
    conditioned = None
    uppers = []
    for tilt in (F(1, 3), F(2, 3)):
        _, counter = _bc_dependency('_am_physical_local')(s, t, z, 1, eta, tilt)
        polynomial = [[[[Fraction(0) for _ in range(5)] for _ in range(5)] for _ in range(5)] for _ in range(5)]
        for u, v, h in product(range(5), repeat=3):
            f = v if u == 0 else 0
            polynomial[f][u][v][h] = counter[u][v][h] / comb(4, u)
        quota = comb(8, 5) * tilt ** 5 * (1 - tilt) ** 3
        exact = [[[Fraction(0) for _ in range(5)] for _ in range(5)] for _ in range(5)]
        for f, u, v in product(range(5), repeat=3):
            exact[f][u][v] = sum(polynomial[a][u][k][h1] * polynomial[f - a][k][v][5 - h1]
                for a in range(f + 1) for k in range(5) for h1 in range(1, 5)) / quota
        exact[0][0][0] -= 1
        if conditioned is not None:
            assert exact == conditioned
        conditioned = exact
        unconditioned = _bc_exact_square(_bc_exact_one(4, s, t, z, eta, tilt, 4))
        upper = [[[value / quota for value in row] for row in block] for block in unconditioned]
        assert all(exact[f][u][v] <= upper[f][u][v] for f, u, v in product(range(5), repeat=3))
        uppers.append(upper)
        regional_entries += 125
    minimum = [[[min(upper[f][u][v] for upper in uppers) for v in range(5)] for u in range(5)] for f in range(5)]
    exact_vectors = [[Fraction(int(f == 0 and u == 0)) for u in range(5)] for f in range(5)]
    for _ in range(4):
        exact_vectors = [[(exact_vectors[f][0] if v == 0 else 0)
                          + sum(exact_vectors[f - b][u] * minimum[b][u][v]
                                for b in range(f + 1) for u in range(5))
                          for v in range(5)] for f in range(5)]
    assert exact_vectors[0] == [1, 0, 0, 0, 0]
    logs = [np.log2(BC_LD(value.numerator) / value.denominator) for value in (s, t, z, eta)]
    blocks = actual_birth_one_step(8, 4, *logs, BC_LD(p.numerator) / p.denominator)
    exact8 = _bc_exact_one(8, s, t, z, eta, p, 4)
    for f, block in enumerate(blocks):
        values = block[0] * np.exp2(block[1]) if block is not None else np.zeros((9, 9))
        expected = np.array([[BC_LD(v.numerator) / v.denominator for v in row] for row in exact8[f]])
        assert np.array_equal(values == 0, expected == 0)
        assert np.max(np.abs(values - expected) / np.maximum(1, expected)) < BC_LD('1e-11')
    regional, _ = actual_birth_regional(4, 2, 5, 4, (BC_LD(1)/3, BC_LD(2)/3), *logs)
    moments = actual_birth_global_vectors(regional, 4)
    for f in range(5):
        expected = sum(exact_vectors[f])
        logarithm = np.log2(BC_LD(expected.numerator) / expected.denominator)
        assert abs(moments[f] - logarithm) < BC_LD('1e-10')
    assert moments[0] == 0.0
    floored, _ = actual_birth_regional(4, 2, 5, 4, (BC_LD(1)/3, BC_LD(2)/3), *logs, floor_bits=2)
    assert all(floored[0][0][0] == 0)
    floor_moments = actual_birth_global_vectors(floored, 4, floor_bits=2)
    assert floor_moments[0] == 0.0 and all(a >= b - 1e-10 for a, b in zip(floor_moments, moments))
    print(dict(exact_physical_birth_entries=physical_entries, exact_regional_birth_checks=regional_entries,
               floating_one_step_entries=405, floating_global_birth_degrees=5,
               exact_degree_zero_control=True, positive_birth_return_excursion=True,
               structural_zero_floors_checked=True, full_code_certificate=False), flush=True)


def actual_birth_main():
    parser = argparse.ArgumentParser(description=__doc__)
    mode = parser.add_mutually_exclusive_group(required=True)
    mode.add_argument('--checks', action='store_true')
    mode.add_argument('--probe', action='store_true')
    parser.add_argument('--fmax', type=int, default=4)
    parser.add_argument('--exact-f', type=int, default=4)
    parser.add_argument('--sum-births', action='store_true')
    parser.add_argument('--p-values', type=float, nargs='+', default=(1/3,))
    parser.add_argument('--log-s', type=float, default=1)
    parser.add_argument('--log-t', '--log-death-marker', dest='log_t', type=float, default=-1)
    parser.add_argument('--log-z', type=float, default=0)
    parser.add_argument('--log-eta', type=float, default=-2)
    parser.add_argument('--exponent', type=float, default=16)
    parser.add_argument('--theta', type=float, default=1/16)
    parser.add_argument('--interval', nargs=2, type=int, action='append')
    parser.add_argument('--floor-majorant', action='store_true')
    args = parser.parse_args()
    import resource
    if shutil.disk_usage('/tmp').free < 20 * (1 << 30):
        raise RuntimeError('20GiB /tmp free-space guard')
    resource.setrlimit(resource.RLIMIT_AS, (2 * (1 << 30), 2 * (1 << 30)))
    timer = threading.Timer(175, lambda: os._exit(124))
    timer.daemon = True
    timer.start()
    started = time.monotonic()
    try:
        if args.checks:
            actual_birth_checks()
        else:
            result = actual_birth_probe(fmax=args.fmax, exact_f=args.exact_f, sum_births=args.sum_births,
                intervals=args.interval or BC_CENTRAL, p_values=tuple(args.p_values),
                log_s=args.log_s, log_t=args.log_t, log_z=args.log_z, log_eta=args.log_eta,
                exponent=args.exponent, theta=args.theta,
                floor_bits=4000 if args.floor_majorant else None, progress=True)
            print(dict(seconds=round(time.monotonic() - started, 3), **result), flush=True)
    finally:
        timer.cancel()


if __name__ == '__main__':
    actual_birth_main()
