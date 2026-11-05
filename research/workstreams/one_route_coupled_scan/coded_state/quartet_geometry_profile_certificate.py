"""Exact geometry-only whole-profile class, with the free terminal retained.

Inject quartet_profile_certificate.py definitions with its main block
removed. Only its generic exact routing, binomial cache and integer
ceil-log2 helpers are used: NO eigenvector or run-moment witness is needed.

For any fixed sequence of quartet partitions and arbitrary masks
S_0=empty,S_1,...,S_T on n coordinates, put
  L=sum|S_i|, G=sum|Gamma_i(S_(i-1))|, J=n(T-1)/4.
At closure marker t=1 their EXACT bivariate counting polynomial is
  (1+s)^n [1+z^4((1+s)^4-1)]^J.
Each of the T-1 nonterminal masks is counted in its next partition;
S_T contributes (1+s)^n and is not expanded. Consequently, for G=4k,
any positive integers a,b give the coefficient bound
  count(L,G) <= C(J,k) ((a+b)^4-b^4)^k (a+b)^n
                     / (a^L b^(G+n-L)).

This counts arbitrary masks and therefore applies separately to actual
state-support or canonical-availability profiles. For a fixed profile
(h,R,L,f,F), set G=F-f and multiply by the generic forced-birth routing
ratio, capped at one. Actual f and the additional canonical constraints
are relaxed in the history count. No field-probability factor is needed;
it is bounded by one. No availability-only screens are used here.

Every ratio and ceiling below is verified with integers. Taking the
minimum over a fixed finite s=a/b family bounds the SAME profile and
does not incur a family union penalty. There are strictly fewer than
2^96 nonzero profiles in the physical ranges (in particular f<=L).
Thus ALL profiles whose integer exponent is <=-172 form a class with
total failure probability <2^-76. This is a class theorem, NOT coverage
of every remaining profile or a full-code distance certificate. The
actual-state and availability applications are separate class bounds;
combining them requires the project's existing probability accounting.

The --checks mode uses independent small labeled-mask censuses and exact
coefficient/Chernoff comparisons. No local imports, subprocesses, data
dumps or floating mathematical comparisons. Run all jobs serially.
"""

import argparse
from collections import Counter
from fractions import Fraction
from math import comb, gcd
from operator import index
import os
import shutil
import threading
import time


GPC_WIDTH = 512
GPC_UPDATES = 4096
GPC_N = GPC_WIDTH * GPC_UPDATES
GPC_K = 1 << 20
GPC_ROWS = 8192
GPC_REGIONS = 256
GPC_UNION_BITS = 96
GPC_THRESHOLD = -172
GPC_MARKERS = ((1, 2), (1, 1), (3, 2), (2, 1), (4, 1), (8, 1))


def _gpc_require_helpers():
    for name in ('profile_certificate_route', 'profile_certificate_ceil_log2_ratio',
                 '_profile_certificate_choose'):
        if name not in globals():
            raise RuntimeError('Inject quartet_profile_certificate.py definitions first; missing ' + name)


def geometry_profile_constants():
    assert GPC_N == 1 << 21 and GPC_WIDTH % 4 == 0
    assert GPC_ROWS * GPC_K * GPC_N ** 3 == 1 << GPC_UNION_BITS
    # The ordered pair (f,L) already makes the physical profile count strict.
    assert GPC_N * (GPC_N + 1) // 2 < GPC_N ** 2
    assert GPC_THRESHOLD + GPC_UNION_BITS == -76
    assert all(a > 0 and b > 0 and gcd(a, b) == 1 for a, b in GPC_MARKERS)
    return dict(width=GPC_WIDTH, updates=GPC_UPDATES, length=GPC_N,
                nonterminal_quartets=(GPC_N - GPC_WIDTH) // 4,
                free_terminal=True, profile_union_bits=GPC_UNION_BITS,
                qualifying_threshold=GPC_THRESHOLD, class_failure_exponent=-76,
                exact_integer_arithmetic=True, eigenvector_witness_required=False)


def _gpc_parameters(L, G, a, b, width, updates):
    L, G, a, b, width, updates = map(index, (L, G, a, b, width, updates))
    if width < 4 or width % 4 or updates < 1:
        raise ValueError('Require positive updates and width divisible by four')
    if a <= 0 or b <= 0:
        raise ValueError('Marker numerator and denominator must be positive integers')
    divisor = gcd(a, b)
    a, b = a // divisor, b // divisor
    J = width * (updates - 1) // 4
    if G < 0 or G % 4 or not 0 <= G // 4 <= J:
        raise ValueError('G must be four times a feasible nonterminal quartet count')
    k = G // 4
    if not k <= L <= G + width:
        raise ValueError('L lies outside the exact coefficient support')
    return L, G, a, b, width, updates, J, k


def geometry_profile_count_ratio(L, G, a, b, width=GPC_WIDTH, updates=GPC_UPDATES, cache=None):
    """Positive integer numerator/denominator of one coefficient upper bound.

    The pair is left unreduced; no expensive giant-integer gcd is needed.
    Keep it in memory only, and use geometry_profile_count_bound for output.
    """
    _gpc_require_helpers()
    L, G, a, b, width, updates, J, k = _gpc_parameters(L, G, a, b, width, updates)
    if cache is None:
        cache = {}
    choose = _profile_certificate_choose(J, k, cache)
    hit_polynomial = (a + b) ** 4 - b ** 4
    numerator = choose * pow(hit_polynomial, k) * pow(a + b, width)
    denominator = pow(a, L) * pow(b, G + width - L)
    return numerator, denominator


def geometry_profile_count_bound(L, G, a, b, width=GPC_WIDTH, updates=GPC_UPDATES, cache=None):
    """Small metadata only; the bound uses exact powers and a minimal ceil-log."""
    L, G, a, b, width, updates, J, k = _gpc_parameters(L, G, a, b, width, updates)
    numerator, denominator = geometry_profile_count_ratio(L, G, a, b, width, updates, cache)
    bits = profile_certificate_ceil_log2_ratio(numerator, denominator)
    return dict(marker_numerator=a, marker_denominator=b, geometry_count_bits=bits,
                nonterminal_quartets=J, hit_quartets=k,
                numerator_bit_length=numerator.bit_length(),
                denominator_bit_length=denominator.bit_length())


def geometry_profile_screen(h, R, L, f, F):
    """Common necessary support conditions; no availability-specific cuts."""
    h, R, L, f, F = map(index, (h, R, L, f, F))
    C_H, C_outside = F - R + 128 * h, R - 128 * h
    G = F - f
    failures = []
    if not (1 <= h <= GPC_ROWS and 128 * h <= R <= GPC_K):
        failures.append('rank_and_row_ranges')
    if not (1 <= f <= L <= F <= GPC_N):
        failures.append('nonzero_profile_ranges')
    if not (129 * h <= C_H <= 256 * h and 0 <= C_outside <= 128 * (GPC_ROWS - h)):
        failures.append('covered_row_capacity')
    if not (f <= C_H and 0 <= 256 * h - C_H <= GPC_N - F):
        failures.append('forced_birth_routing_support')
    if not (0 <= G <= GPC_N - GPC_WIDTH and G % 4 == 0 and G // 4 <= L <= G + GPC_WIDTH):
        failures.append('quartet_and_free_terminal_support')
    return dict(h=h, R=R, L=L, f=f, F=F, G=G, C_H=C_H, C_outside=C_outside,
                empty=bool(failures), failed_necessary_conditions=failures,
                availability_specific_screens_used=False)


def geometry_profile_point(h, R, L, f, F, markers=GPC_MARKERS, cache=None):
    """Exact class criterion, simultaneously summing all death counts."""
    _gpc_require_helpers()
    constants = geometry_profile_constants()
    screen = geometry_profile_screen(h, R, L, f, F)
    if screen['empty']:
        raise ValueError('Profile fails necessary conditions: ' + ', '.join(screen['failed_necessary_conditions']))
    if not markers:
        raise ValueError('At least one positive rational marker is required')
    if cache is None:
        cache = {}
    candidates = [geometry_profile_count_bound(screen['L'], screen['G'], a, b, cache=cache)
                  for a, b in markers]
    best = min(candidates, key=lambda candidate: candidate['geometry_count_bits'])
    routing = profile_certificate_route(screen['h'], screen['R'], screen['f'], screen['F'], cache)
    exponent = best['geometry_count_bits'] + routing['probability_capped_bits']
    return dict(kind='exact_geometry_profile_class', **screen,
                geometry_count_bits=best['geometry_count_bits'],
                chosen_marker_numerator=best['marker_numerator'],
                chosen_marker_denominator=best['marker_denominator'],
                marker_candidates=candidates,
                exact_profile_exponent_upper=exponent,
                qualifies_global_profile_class=exponent <= GPC_THRESHOLD,
                qualifying_threshold=GPC_THRESHOLD, global_profile_union_bits=GPC_UNION_BITS,
                resulting_class_failure_exponent=constants['class_failure_exponent'],
                no_marker_family_union=True, all_death_counts_summed=True,
                free_terminal=True, eigenvector_witness_required=False,
                field_probability_bounded_by_one=True, exact_integer_arithmetic=True,
                full_code_certificate=False,
                routebits=routing['routebits'],
                probability_capped_routebits=routing['probability_capped_bits'],
                raw_route_ratio_at_most_one=routing['raw_ratio_at_most_one'],
                route_numerator_bit_length=routing['numerator_bit_length'],
                route_denominator_bit_length=routing['denominator_bit_length'])


def _gpc_convolve(left, right):
    result = [0] * (len(left) + len(right) - 1)
    for i, value in enumerate(left):
        for j, other in enumerate(right):
            result[i + j] += value * other
    return result


def _gpc_exact_coefficients(width, updates):
    J = width * (updates - 1) // 4
    free_terminal = [comb(width, L) for L in range(width + 1)]
    hit_polynomial = [1]
    coefficients = {}
    for k in range(J + 1):
        for L, value in enumerate(_gpc_convolve(hit_polynomial, free_terminal)):
            if value:
                coefficients[L, 4 * k] = comb(J, k) * value
        hit_polynomial = _gpc_convolve(hit_polynomial, (0, 4, 6, 4, 1))
    return coefficients


def _gpc_mask_census(width, updates):
    """Independent labeled masks under fixed disjoint quartets."""
    local = Counter()
    for mask in range(1 << width):
        hit_size = sum(4 for start in range(0, width, 4) if mask & (15 << start))
        local[mask.bit_count(), hit_size] += 1
    histories = Counter({(0, 0): 1})
    for _ in range(updates - 1):
        following = Counter()
        for (L, G), value in histories.items():
            for (new_L, new_G), multiplicity in local.items():
                following[L + new_L, G + new_G] += value * multiplicity
        histories = following
    terminal = Counter(mask.bit_count() for mask in range(1 << width))
    result = Counter()
    for (L, G), value in histories.items():
        for final_L, multiplicity in terminal.items():
            result[L + final_L, G] += value * multiplicity
    return dict(result)


def geometry_profile_checks():
    _gpc_require_helpers()
    constants = geometry_profile_constants()
    coefficients = chernoff = ceilings = terminals = 0
    for width, updates in ((4, 1), (4, 2), (4, 3), (8, 1), (8, 2)):
        exact = _gpc_exact_coefficients(width, updates)
        truth = _gpc_mask_census(width, updates)
        assert exact == truth
        assert sum(truth.values()) == 1 << (width * updates)
        if updates == 1:
            assert truth == {(L, 0): comb(width, L) for L in range(width + 1)}
            terminals += 1
        cache = {}
        for (L, G), count in truth.items():
            coefficients += 1
            for a, b in GPC_MARKERS:
                numerator, denominator = geometry_profile_count_ratio(L, G, a, b, width, updates, cache)
                assert count * denominator <= numerator
                bits = profile_certificate_ceil_log2_ratio(numerator, denominator)
                if bits >= 0:
                    assert numerator <= denominator << bits
                else:
                    assert numerator << -bits <= denominator
                predecessor = bits - 1
                if predecessor >= 0:
                    assert numerator > denominator << predecessor
                else:
                    assert numerator << -predecessor > denominator
                chernoff += 1
                ceilings += 1
    # A rational marker is invariant under common integer rescaling.
    for a, b in GPC_MARKERS:
        reference = geometry_profile_count_ratio(5, 8, a, b, 4, 3)
        scaled = geometry_profile_count_ratio(5, 8, 3 * a, 3 * b, 4, 3)
        assert reference == scaled
    screen = geometry_profile_screen(3584, 524288, 526950, 4, 954368)
    assert not screen['empty']
    print(dict(exact_mask_coefficient_checks=coefficients,
               exact_rational_chernoff_checks=chernoff, exact_minimal_ceil_checks=ceilings,
               free_terminal_checks=terminals, marker_scale_checks=len(GPC_MARKERS),
               default_profile_support_checked=True,
               profile_union_bits=constants['profile_union_bits'],
               no_eigenvector_witness=True, full_code_certificate=False), flush=True)


def geometry_profile_main():
    parser = argparse.ArgumentParser(description=__doc__)
    mode = parser.add_mutually_exclusive_group(required=True)
    mode.add_argument('--checks', action='store_true')
    mode.add_argument('--certify', action='store_true')
    parser.add_argument('--profile', type=int, nargs=5,
                        default=(3584, 524288, 526950, 4, 954368),
                        metavar=('h', 'R', 'L', 'f', 'F'))
    parser.add_argument('--markers', nargs='+', default=('1/2', '1', '3/2', '2', '4', '8'),
                        help='Positive exact rational s markers; no family union is charged')
    args = parser.parse_args()
    import resource
    if shutil.disk_usage('/tmp').free < 20 * (1 << 30):
        raise RuntimeError('20GiB /tmp free-space guard')
    resource.setrlimit(resource.RLIMIT_AS, (2 * (1 << 30), 2 * (1 << 30)))
    timer = threading.Timer(175, lambda: os._exit(124))
    timer.daemon = True
    timer.start()
    try:
        if args.checks:
            geometry_profile_checks()
            return
        rational_markers = tuple(Fraction(value) for value in args.markers)
        if any(marker <= 0 for marker in rational_markers):
            raise ValueError('Every marker must be strictly positive')
        markers = tuple((marker.numerator, marker.denominator) for marker in rational_markers)
        started = time.monotonic()
        result = geometry_profile_point(*args.profile, markers=markers)
        print(dict(**result, seconds=round(time.monotonic() - started, 3)), flush=True)
    finally:
        timer.cancel()


if __name__ == '__main__':
    geometry_profile_main()
