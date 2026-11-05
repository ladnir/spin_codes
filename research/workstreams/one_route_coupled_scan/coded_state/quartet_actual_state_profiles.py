"""Exact profile-class interface for ACTUAL state-support histories.

These statistics are distinct from canonical availability statistics:
  Ls=sum|S_t|, D=sum|S_(t-1) minus S_t|, a=|S_T|,
  V_t=Gamma_t(S_(t-1)) union S_t,
  factual=sum|S_t minus Gamma_t(S_(t-1))|, F=sum|V_t|.
Then G=F-factual, total coordinate births b=D+a, and survivors i=Ls-b.
The actual-S field weight is (qx)^b*x^D*(1+qx)^i.  In particular,
factual must NOT replace b in that formula.

The existing exact run/occupancy witness also counts arbitrary S masks.
At its frozen markers s=75/64, z^4=1/2, d0=d_availability(1/256),
this helper uses six field/geometry branches dominated by the same
d0^D moment.  Availability-only feasibility screens and census class
labels are not imported.  Separately proved actual-S weighted censuses
can remove prior classes and refine the actual death interval.  Such
outputs bound the profile REMAINDER, with prior budgets kept separate.
The --whole-profile control omits those cuts.  The forced-birth routing
ratio is valid for these actual forced births as well.

Inject quartet_run_occupancy.py, quartet_run_certificate.py, and
quartet_profile_certificate.py definitions with their main blocks removed.
--certify freshly recomputes and verifies the moment before using its bit
bound.  The callable evaluator requires an independently verified bound
as an explicit argument.  Integer arithmetic proves the resulting class
criterion; no complete distance or all-profile coverage claim is made.
No local imports, subprocesses, witness files or data dumps.  Run serially.
"""

import argparse
from collections import Counter
from fractions import Fraction
from heapq import heappush, heapreplace
from itertools import combinations, product
from math import factorial, lgamma, log
from operator import index
import os
import shutil
import threading
import time


AS_N = 1 << 21
AS_K = 1 << 20
AS_Q = 1 << 32
AS_W = 943718
AS_WIDTH = 512
AS_UPDATES = 4096
AS_ROWS = 8192
AS_UNION_BITS = 96
AS_THRESHOLD = -172
AS_CENSUS_J = tuple(range(17))


def actual_state_profile_constants():
    """Finite rational checks for every constant and closure comparison."""
    q, N = AS_Q, AS_N
    d0 = Fraction(q + 65536, 256 * (q + 257))
    assert Fraction(1, 256) <= d0 < 1
    assert Fraction(1024, 989) < 2
    assert 75 ** 1024 > 1 << 6378
    assert Fraction(75, 64) ** 1024 > 1 << 234

    x8 = Fraction(1, 256)
    lambda8 = 1 + q * x8
    t8 = q * x8 * x8 / lambda8
    assert lambda8 == (1 << 24) * Fraction((1 << 24) + 1, 1 << 24)
    # (1+2^-24)^Ls <= 1/(1-N*2^-24)=8/7<2.
    assert Fraction(1, 1) / (1 - Fraction(N, 1 << 24)) == Fraction(8, 7) < 2
    assert t8 < Fraction(1, 256) <= d0

    x16 = Fraction(1, 1 << 16)
    lambda16 = 1 + q * x16
    t16 = q * x16 * x16 / lambda16
    assert lambda16 == (1 << 16) * (1 + Fraction(1, 1 << 16))
    assert t16 == Fraction(1, 65537) < Fraction(1, 1 << 16)
    # For m=65536, (1+1/m)^m <= sum_j 1/j!.  Starting with j=5,
    # the tail ratio is <=1/6, so its sum is <=1/100.  Since N=32m,
    # (1+2^-16)^N <(11/4)^32<2^47, using rational checks only.
    series_upper = sum(Fraction(1, factorial(j)) for j in range(5)) + Fraction(1, 100)
    assert series_upper == Fraction(1631, 600) < Fraction(11, 4)
    assert N == 32 * (1 << 16)
    assert 11 ** 32 < 1 << 111
    assert t16 <= d0 * d0                 # theta=1/2, without square roots
    assert t16 / d0 <= Fraction(1, 256)  # plain branch saves 8 bits per death

    x32 = Fraction(1, q)
    lambda32 = 1 + q * x32
    t32 = q * x32 * x32 / lambda32
    assert lambda32 == 2 and t32 == Fraction(1, 1 << 33)
    assert t32 ** 8 <= d0 ** 33           # theta=8/33
    assert t32 / d0 <= Fraction(1, 1 << 25)
    assert 1 / d0 <= 256                 # geometry-only branch
    # Actual x=1 rank bound: gamma*(1+1/q)^N <2.
    assert Fraction(1024, 989) / (1 - Fraction(N, q)) < 2
    assert 9 ** 64 < 1 << 203             # the sharper j=3 binary-run census
    assert AS_ROWS * AS_K * AS_N ** 3 == 1 << AS_UNION_BITS
    assert AS_UNION_BITS + AS_THRESHOLD == -76
    return dict(statistics='actual_state_support', exact_constants_verified=True,
                closure_numerator=d0.numerator, closure_denominator=d0.denominator,
                profile_union_bits=AS_UNION_BITS, qualifying_threshold=AS_THRESHOLD,
                class_failure_exponent=-76, availability_screens_used=False)


def actual_state_profile_screen(h, R, Ls, factual, F):
    """Necessary actual-S conditions only; passing does not prove realizability."""
    h, R, Ls, factual, F = map(index, (h, R, Ls, factual, F))
    C = F - R + 128 * h
    outside = R - 128 * h
    G = F - factual
    dmin = max(0, factual - AS_WIDTH)
    dmax = min(F - Ls, AS_W // 2)
    failures = []
    if not (1 <= h <= AS_ROWS and 128 * h <= R <= AS_K):
        failures.append('rank_and_row_ranges')
    if not (1 <= factual <= Ls <= F <= AS_N):
        failures.append('nonzero_actual_support_ranges')
    if not (129 * h <= C <= 256 * h and 0 <= outside <= 128 * (AS_ROWS - h)):
        failures.append('covered_row_capacity')
    if not (factual <= C and 0 <= 256 * h - C <= AS_N - F):
        failures.append('forced_birth_routing_support')
    if not (G >= 0 and G % 4 == 0 and G <= AS_N - AS_WIDTH
            and G // 4 <= Ls <= G + AS_WIDTH):
        failures.append('actual_state_occupancy_support')
    if dmin > dmax or 2 * factual - AS_WIDTH > AS_W:
        failures.append('actual_birth_death_output_budget')
    return dict(statistics='actual_state_support', h=h, R=R, Ls=Ls,
                factual=factual, F=F, G=G, C_H=C, C_outside=outside,
                dmin=dmin, dmax=dmax, empty=bool(failures),
                failed_necessary_conditions=failures,
                availability_screens_used=False, realizability_asserted=False)


def actual_state_census_constant(j):
    """Exact exponent in sum_S 2^(-2jD)<=2^c_j, with free endpoints.

    The two-state positive vector gives 2^(width*j)*(1+2^-j)^N.
    ln(2)>2/3 follows by integrating 1/t>4/(1+t)^2 from 1 to 2.
    Thus c_j=width*j+ceil(3N/2^(j+1)); j=3 uses 9^64<2^203.
    """
    j = index(j)
    if j < 0:
        raise ValueError('Census index must be nonnegative')
    if j == 0:
        return AS_N
    if j == 3:
        return 361984
    denominator = 1 << (j + 1)
    return AS_WIDTH * j + (3 * AS_N + denominator - 1) // denominator


def actual_state_census_interval(h, R, Ls, factual, F, j_values=AS_CENSUS_J):
    """Intersect complements of independently proved ACTUAL-S prior classes.

    Each rank class 32(R-Ls)-2jD>=c_j+45 and closure class
    32R-16Ls+(16-2j)D>=16W+c_j+92 costs <2^-76, without a
    profile union.  The separate R-Ls>=65537 class costs <2^-63.
    A covered profile is not declared geometrically empty.
    """
    screen = actual_state_profile_screen(h, R, Ls, factual, F)
    if screen['empty']:
        return dict(**screen, status='empty', refined_death_interval=None)
    h, R, Ls, factual, F = (screen[key] for key in ('h', 'R', 'Ls', 'factual', 'F'))
    indices = tuple(sorted(set(index(j) for j in j_values)))
    if any(j < 0 for j in indices):
        raise ValueError('Census indices must be nonnegative')
    prior = dict(actual_census_indices=indices, actual_census_cut_count=2 * len(indices),
                 per_census_cut_failure_exponent=-76,
                 separate_actual_rank_gap_failure_exponent=-63)
    if R - Ls >= 65537:
        return dict(**screen, **prior, status='actual_rank_gap_class',
                    refined_death_interval=None)
    lower, upper = screen['dmin'], screen['dmax']
    rank_left = 32 * (R - Ls)
    closure_base = 32 * R - 16 * Ls
    for j in indices:
        c = actual_state_census_constant(j)
        if j == 0:
            if rank_left > c + 44:
                return dict(**screen, **prior, status='actual_census_classes',
                            refined_death_interval=None)
        else:
            lower = max(lower, -((c + 44 - rank_left) // (2 * j)))
        coefficient = 16 - 2 * j
        right = 16 * AS_W + c + 91 - closure_base
        if coefficient > 0:
            upper = min(upper, right // coefficient)
        elif coefficient < 0:
            lower = max(lower, -((-right) // coefficient))
        elif right < 0:
            return dict(**screen, **prior, status='actual_census_classes',
                        refined_death_interval=None)
        if lower > upper:
            return dict(**screen, **prior, status='actual_census_classes',
                        refined_death_interval=None)
    return dict(**screen, **prior, status='remaining', refined_death_interval=(lower, upper))


def actual_state_field_branches(R, Ls, factual, F, death_interval=None):
    """Integer residual exponents after extracting the common d0^D moment.

    For x>0 put lambda_x=1+qx and t_x=qx^2/lambda_x.  The capped
    per-S field bound is at most min(1,B_x*t_x^D), where
      B_x=(gamma/q)*x^-W*q^-R*lambda_x^Ls.
    The discarded terminal factor (qx/lambda_x)^a is <=1.
    Raising the conditional probability to theta in [0,1] preserves an
    upper bound.  Routing is kept separate and is never raised to theta.
    """
    R, Ls, factual, F = map(index, (R, Ls, factual, F))
    dmin = max(0, factual - AS_WIDTH)
    dmax = min(F - Ls, AS_W // 2)
    if death_interval is not None:
        lower, upper = map(index, death_interval)
        if not dmin <= lower <= upper <= dmax:
            raise ValueError('Refined death interval must lie inside the physical interval')
        dmin, dmax = lower, upper
    b32 = -31 + 32 * (AS_W - R) + Ls
    return dict(
        x8_plain=-30 + 8 * AS_W - 32 * R + 24 * Ls,
        x16_half=8 + 8 * AS_W - 16 * R + 8 * Ls,
        x16_plain=16 + 16 * AS_W - 32 * R + 16 * Ls - 8 * dmin,
        x32_eight_over_33=-((-8 * b32) // 33),
        x32_plain=b32 - 25 * dmin,
        geometry_only=8 * dmax)


def actual_state_profile_certificate(h, R, Ls, factual, F, moment_bits, cache=None,
                                     use_census=True):
    """Exact class calculation conditional on a valid supplied moment bound.

    The caller must supply a verified log2 integer upper bound for the
    width-512, 4096-update moment at s=75/64, z^4=1/2 and the stated d0.
    --certify obtains it by a fresh all-row integer verification.
    No availability profile certificate or availability class is reused.
    With use_census=True (default), the exponent bounds only the remainder
    outside the separately budgeted actual rank-gap and census classes.
    """
    moment_bits = index(moment_bits)
    if moment_bits < 0:
        raise ValueError('A bound for the moment including the empty history must be nonnegative')
    screen = actual_state_profile_screen(h, R, Ls, factual, F)
    result = dict(**screen, moment_log2_integer_upper=moment_bits,
                  conditional_on_supplied_moment_bound=True,
                  global_profile_union_bits=AS_UNION_BITS,
                  all_actual_death_counts=True, free_endpoint=True,
                  full_code_certificate=False)
    if screen['empty']:
        return dict(**result, status='empty', exact_profile_exponent_upper=None,
                    qualifies_global_profile_class=False)
    interval = None
    if use_census:
        census = actual_state_census_interval(h, R, Ls, factual, F)
        result.update(status=census['status'],
                      scope='profile_remainder_outside_actual_prior_classes',
                      actual_census_indices=census['actual_census_indices'],
                      actual_census_cut_count=census['actual_census_cut_count'],
                      per_census_cut_failure_exponent=-76,
                      separate_actual_rank_gap_failure_exponent=-63)
        if census['status'] != 'remaining':
            return dict(**result, exact_profile_remainder_exponent_upper=None,
                        qualifies_global_profile_class=False)
        interval = census['refined_death_interval']
    else:
        result.update(status='remaining', scope='whole_actual_state_profile')
    if 'profile_certificate_route' not in globals():
        raise RuntimeError('Inject the generic exact forced-birth routing function')
    h, R, Ls, factual, F = (screen[key] for key in ('h', 'R', 'Ls', 'factual', 'F'))
    routing = profile_certificate_route(h, R, factual, F, cache)
    base = (moment_bits - (117 * Ls) // 512 + (F - factual) // 4
            + routing['probability_capped_bits'])
    branches = actual_state_field_branches(R, Ls, factual, F, interval)
    branch, residual = min(branches.items(), key=lambda entry: entry[1])
    exponent = base + residual
    # At most 2^96 profiles exist.  The pointwise minimum is chosen inside
    # each profile; it costs no branch union.  D is already summed.
    exponent_key = 'exact_profile_remainder_exponent_upper' if use_census else 'exact_profile_exponent_upper'
    result[exponent_key] = exponent
    return dict(**result, routebits=routing['routebits'],
                route_probability_capped_bits=routing['probability_capped_bits'],
                common_base_exponent=base, field_branch_exponents=branches,
                selected_branch=branch, refined_death_interval=interval,
                qualifies_global_profile_class=exponent <= AS_THRESHOLD,
                qualifying_threshold=AS_THRESHOLD, class_failure_exponent=-76,
                no_branch_or_death_count_union=True)


def _as_logchoose_float(n, k):
    if k < 0 or k > n:
        return float('-inf')
    if k == 0 or k == n:
        return 0.0
    return (lgamma(n + 1) - lgamma(k + 1) - lgamma(n - k + 1)) / log(2.0)


def actual_state_float_value(h, R, Ls, factual, F, moment_bits, use_census=True):
    """FLOAT scouting value; no availability screen or class label is used.

    Besides the actual-S classes, the independent covered-footprint and
    absolute-rank classes F>=1736704 and R>=1009256 are removed.  No
    availability-defined rank gap, census or feasibility condition is used.
    """
    screen = actual_state_profile_screen(h, R, Ls, factual, F)
    if screen['empty']:
        return dict(status='empty', statistics='actual_state_support')
    if F >= 1736704 or R >= 1009256:
        return dict(status='independent_prior_class', statistics='actual_state_support')
    interval = None
    if use_census:
        census = actual_state_census_interval(h, R, Ls, factual, F)
        if census['status'] != 'remaining':
            return dict(status=census['status'], statistics='actual_state_support')
        interval = census['refined_death_interval']
    C = screen['C_H']
    route = min(0.0, _as_logchoose_float(F - factual, C - factual)
                + _as_logchoose_float(AS_N - F, 256 * h - C)
                - 255 * _as_logchoose_float(AS_ROWS, h))
    branches = actual_state_field_branches(R, Ls, factual, F, interval)
    branch, residual = min(branches.items(), key=lambda entry: entry[1])
    base = moment_bits - (117 * Ls) // 512 + (F - factual) / 4 + route
    fixed_value = base + residual
    occupancy_only = None
    value = fixed_value
    if 'occupancy_geometry_logbound' in globals():
        # The t=1 helper counts arbitrary actual S masks as well.  This is
        # an optional FLOAT comparison, never an input to the exact evaluator.
        occupancy_only = occupancy_geometry_logbound(Ls, F - factual) + route
        value = min(value, occupancy_only)
    return dict(status='remaining', statistics='actual_state_support',
                h=h, R=R, Ls=Ls, factual=factual, F=F,
                actual_death_interval=(screen['dmin'], screen['dmax']),
                refined_death_interval=interval,
                route=route, common_base=base, selected_branch=branch,
                residual_exponent=residual, fixed_moment_value=fixed_value,
                occupancy_geometry_value=occupancy_only, value=value,
                float_diagnostic_only=True, availability_screens_used=False,
                actual_census_used=use_census,
                full_code_certificate=False)


def actual_state_profile_probe(moment_bits, top=8, use_census=True):
    """Bounded coarse FLOAT scan; a finite grid cannot certify its complement."""
    moment_bits, top = index(moment_bits), index(top)
    if moment_bits < 0 or not 1 <= top <= 32:
        raise ValueError('Require a nonnegative moment bound and top in 1..32')
    counts = Counter()
    worst = []
    sequence = 0
    hvalues = (1, 4, 16, 64, 128, 256, 512, 768, 1024, 1536, 2048,
               2560, 3072, 3584, 4096, 5120, 6144, 7168)
    for h in hvalues:
        local = []
        candidates = 0
        ranks = sorted({128 * h, 1009255,
                        *range(max(128 * h, 32768), 1009256, 32768)})
        for R in ranks:
            if R > 1009255:
                continue
            for extra in (1, 16, 32, 64, 96, 112, 120, 124, 128):
                F = R + h * extra
                if F >= 1736704:
                    counts['independent_prior_F_class'] += 1
                    continue
                # F=G+factual<=4Ls+Ls gives this actual-S lower limit.
                lower_L = max(1, (F + 4) // 5)
                states = sorted({lower_L, R, AS_W, F, F // 2,
                                 R + (F - AS_W) // 4,
                                 *range(lower_L, F + 1, 32768)})
                for Ls in states:
                    if not 1 <= Ls <= F:
                        continue
                    minimum_f = max(1, F - 4 * Ls, F - AS_N + AS_WIDTH)
                    maximum_f = min(Ls, F - R + 128 * h,
                                    F - Ls + AS_WIDTH, (AS_W + AS_WIDTH) // 2)
                    for raw_f in sorted({minimum_f, maximum_f, 512, 4096,
                                         16384, 33280, 65536, 131072, 262144}):
                        factual = raw_f + (F - raw_f) % 4
                        if not minimum_f <= factual <= maximum_f:
                            continue
                        entry = actual_state_float_value(h, R, Ls, factual, F, moment_bits, use_census)
                        counts[entry['status']] += 1
                        if entry['status'] != 'remaining':
                            continue
                        candidates += 1
                        counts['float_negative' if entry['value'] <= AS_THRESHOLD else 'float_open'] += 1
                        sequence += 1
                        item = entry['value'], sequence, entry
                        if len(local) < top:
                            heappush(local, item)
                        elif item[0] > local[0][0]:
                            heapreplace(local, item)
        worst.extend(local)
        print(dict(h=h, candidates=candidates,
                   worst=max(local, key=lambda item: item[0])[2] if local else None), flush=True)
    print(dict(statistics='actual_state_support', counts=dict(counts),
               finite_grid_only=True, availability_screens_used=False,
               actual_census_used=use_census,
               full_code_certificate=False), flush=True)
    for _, _, entry in sorted(worst, reverse=True)[:top]:
        print(entry, flush=True)


def _as_partitions(width):
    if width == 4:
        return [(15,)]
    if width == 8:
        return [(first, 255 ^ first) for rest in combinations(range(1, 8), 3)
                for first in [1 + sum(1 << j for j in rest)]]
    raise ValueError('Physical identity checks use widths four and eight')


def actual_state_profile_checks():
    constants = actual_state_profile_constants()
    histories = 0
    field_statistics = set()
    violates_availability_screen = False
    for width, updates in ((4, 0), (4, 1), (4, 2), (4, 3), (8, 1), (8, 2)):
        for partition in _as_partitions(width):
            full = (1 << width) - 1
            gamma = [sum(quartet for quartet in partition if quartet & mask)
                     for mask in range(1 << width)]
            sizes = [mask.bit_count() for mask in range(1 << width)]
            for masks in product(range(1 << width), repeat=updates):
                old = Ls = D = b = i = factual = F = G = 0
                for new in masks:
                    covered = gamma[old]
                    Ls += sizes[new]
                    D += sizes[old & (full ^ new)]
                    b += sizes[new & (full ^ old)]
                    i += sizes[new & old]
                    factual += sizes[new & (full ^ covered)]
                    F += sizes[covered | new]
                    G += sizes[covered]
                    old = new
                terminal = sizes[old]
                assert b == D + terminal and i == Ls - b
                assert G == F - factual and G % 4 == 0
                assert factual <= b and max(0, factual - width) <= D <= F - Ls
                assert 2 * D + terminal == b + D
                assert 2 * factual - width <= b + D
                if updates:
                    assert G <= width * (updates - 1)
                    assert G // 4 <= Ls <= G + width
                    violates_availability_screen |= 2 * F - Ls - factual > width * (updates - 1)
                else:
                    assert (Ls, D, b, factual, F, G, terminal) == (0, 0, 0, 0, 0, 0, 0)
                field_statistics.add((Ls, D, b, i, terminal))
                histories += 1
    # Actual S histories can violate this canonical-availability inequality.
    # This assertion guards the distinction; no such screen is imported.
    assert violates_availability_screen
    rational_fields = 0
    for q in (4, 16, AS_Q):
        for x in (Fraction(1), Fraction(1, 4), Fraction(1, 256)):
            lam = 1 + q * x
            death = q * x * x / lam
            for Ls, D, b, i, terminal in field_statistics:
                direct = (q * x) ** b * x ** D * lam ** i
                rewritten = lam ** Ls * death ** D * (q * x / lam) ** terminal
                assert direct == rewritten and (q * x / lam) ** terminal <= 1
                rational_fields += 1
    rounding = 0
    for value in range(-256, 257):
        exponent = -((-8 * value) // 33)
        assert Fraction(8 * value, 33) <= exponent < Fraction(8 * value, 33) + 1
        rounding += 1
    # Branch evaluation needs no route or witness and uses actual factual.
    sample = actual_state_field_branches(900000, 900000, 20000, 1150000)
    assert sample['geometry_only'] == 8 * min(250000, AS_W // 2)
    assert sample['x16_plain'] == 16 + 16 * AS_W - 32 * 900000 + 16 * 900000 - 8 * 19488
    # Exact signed complement arithmetic, including the j=8 zero coefficient.
    complements = 0
    for coefficient in range(-18, 19, 2):
        for right in range(-51, 52):
            for death in range(31):
                if coefficient > 0:
                    retained = death <= right // coefficient
                elif coefficient < 0:
                    retained = death >= -((-right) // coefficient)
                else:
                    retained = right >= 0
                assert retained == (coefficient * death <= right)
                complements += 1
    assert actual_state_census_constant(0) == AS_N
    assert actual_state_census_constant(3) == 361984
    refined = actual_state_field_branches(900000, 900000, 20000, 1150000, (30000, 60000))
    assert refined['x16_plain'] == sample['x16_plain'] - 8 * (30000 - 19488)
    assert refined['x32_plain'] == sample['x32_plain'] - 25 * (30000 - 19488)
    assert refined['geometry_only'] == 8 * 60000
    print(dict(exact_actual_support_histories=histories,
               exact_actual_field_identities=rational_fields,
               signed_fractional_ceiling_checks=rounding,
               signed_census_complement_checks=complements,
               availability_counterexample_checked=True,
               moment_recomputed=False, full_code_certificate=False, **constants), flush=True)


def actual_state_profile_main():
    parser = argparse.ArgumentParser(description=__doc__)
    mode = parser.add_mutually_exclusive_group(required=True)
    mode.add_argument('--checks', action='store_true')
    mode.add_argument('--certify', action='store_true')
    mode.add_argument('--probe', action='store_true')
    parser.add_argument('--profile', type=int, nargs=5, action='append',
                        metavar=('H', 'R', 'LS', 'FACTUAL', 'F'),
                        help='Repeatable actual-state profile; default certifies the original and dense test points')
    parser.add_argument('--top', type=int, default=8)
    parser.add_argument('--whole-profile', action='store_true',
                        help='Control: omit actual-S census and rank-gap cuts')
    parser.add_argument('--iterations', type=int, default=256)
    parser.add_argument('--mantissa-bits', type=int, default=48)
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
            actual_state_profile_checks()
            return
        actual_state_profile_constants()
        if 'profile_certificate_recompute_moment' not in globals():
            raise RuntimeError('Inject the existing exact moment recomputation function')
        witness = profile_certificate_recompute_moment(args.iterations, args.mantissa_bits, progress=True)
        print(dict(stage='actual_state_moment_verified', **witness), flush=True)
        moment_bits = witness['moment_log2_integer_upper']
        if args.probe:
            actual_state_profile_probe(moment_bits, args.top, not args.whole_profile)
            return
        profiles = args.profile or [(2048, 900000, 900000, 20000, 1150000),
                                    (7168, 917504, 1310720, 33280, 1720320)]
        cache = {}
        for profile in profiles:
            result = actual_state_profile_certificate(*profile, moment_bits, cache,
                                                       use_census=not args.whole_profile)
            print(dict(seconds=round(time.monotonic() - started, 3),
                       moment_recomputed=True, **result), flush=True)
    finally:
        timer.cancel()


if __name__ == '__main__':
    actual_state_profile_main()
