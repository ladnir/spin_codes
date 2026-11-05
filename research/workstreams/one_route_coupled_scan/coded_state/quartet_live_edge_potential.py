"""Bounded five-state FLOAT control for live-selector temporal overlap.

Fix all quartet partitions. Each quartet has an iid eligible count j with
p_j=C(4,j)p^j(1-p)^(4-j), before conditioning regional slot quotas. For
binary adjacent selector bits a,b define

 g(k,j)=max_(a,b) [ab+(mu/8)(a(k-r)+b(j-r))].

The integer exponent e=4g is max(0,(mu/2)(k-r),(mu/2)(j-r),
4+(mu/2)(k+j-2r)). Half-integer r is supported; then mu is a multiple4.
Set L_kj=p_j*beta^e_kj. Any positive v and lambda>=max(Lv/v) give

 log2 E beta^sum_edges g
 <=128(T-1)log2lambda+128log2(sum_j p_j*v_j/min(v)).

This follows by Holder at each new quartet and its four incoming stubs.
It is conditional on every partition sequence and requires no support
union. The necessary actual-profile score is

 Ls-512-f+mu*((1-r/4)*(F-f)-(R-128h))
       -(mu/2)*128*(4-r).

The first selector layer is empty. The two subtractions retain the free
terminal and last-layer edge boundary. The probability bound adds one
H-set union and256 iid quota penalties, then subtracts log2(beta)*score.

The probe uses ordinary floating eigenvectors/Collatz ratios and is NOT
a directed certificate. Exact tiny checks use integer/Fraction arithmetic.
Combining cut classes for distinct(h,mu,r) requires a finite union;
p,beta,v choices for one fixed statistic require no additional union.
No geometry helper or legacy source is imported or modified.
"""

import argparse
from collections import Counter
from fractions import Fraction
from itertools import combinations, combinations_with_replacement, product
from math import comb, isfinite, lgamma, log, log2
from operator import index
import os
import shutil
import threading
import time


LE_M = 8192
LE_WIDTH = 512
LE_UPDATES = 4096
LE_GROUPS = LE_WIDTH // 4
LE_REGIONS = 256
LE_PROFILES = ((2048, 819200, 819200, 4, 950272),
               (1536, 753664, 730726, 4, 851968),
               (2560, 786432, 779879, 4, 950272))


def live_edge_exponents(r=Fraction(3, 2), mu=8):
    """Exact5x5 exponents; accepted floating r values must be half-integers."""
    r = Fraction(str(r)) if isinstance(r, float) else Fraction(r)
    mu = index(mu)
    if not 0 <= r <= 4 or (2 * r).denominator != 1 or mu <= 0 or mu % 2:
        raise ValueError('Require half-integer0<=r<=4 and positive even mu')
    if r.denominator == 2 and mu % 4:
        raise ValueError('Half-integer r requires mu divisible by4')
    result = []
    for k in range(5):
        row = []
        for j in range(5):
            value = max(Fraction(0), Fraction(mu, 2) * (k - r),
                        Fraction(mu, 2) * (j - r),
                        4 + Fraction(mu, 2) * (k + j - 2 * r))
            assert value.denominator == 1
            row.append(value.numerator)
        result.append(tuple(row))
    return tuple(result)


def live_edge_score(profile=LE_PROFILES[0], r=Fraction(3, 2), mu=8):
    """Exact actual-profile statistic; no availability-only screen is used."""
    h, rank, volume, births, footprint = map(index, profile)
    r = Fraction(str(r)) if isinstance(r, float) else Fraction(r)
    live_edge_exponents(r, mu)
    if not (1 <= h < LE_M and rank >= 128 * h
            and 0 <= births <= volume <= footprint <= LE_WIDTH * LE_UPDATES):
        raise ValueError('Invalid actual-state profile')
    G, outside = footprint - births, rank - 128 * h
    if G % 4:
        raise ValueError('Actual G=F-f must be divisible by four')
    z = (4 - r) * (G // 4) - outside
    boundary = Fraction(mu, 2) * LE_GROUPS * (4 - r)
    score = volume - LE_WIDTH - births + mu * z - boundary
    return dict(z=z, boundary=boundary, score=score, G=G, outside=outside)


def _le_falling(value, length):
    result = 1
    for i in range(length):
        result *= value - i
    return result


def _le_mean_vertices():
    """All arrangement vertices on mu>=0,0<=nu<=4mu, including the origin."""
    lines = {(Fraction(k), Fraction(0)) for k in range(5)}
    lines.update((Fraction(j + k, 2), Fraction(4))
                 for j in range(5) for k in range(j, 5))
    lines.update((Fraction(j), Fraction(8)) for j in range(5))
    vertices = {(Fraction(0), Fraction(0))}
    for (a, b), (c, d) in combinations(sorted(lines), 2):
        if a == c:
            continue
        mu = (d - b) / (a - c)
        nu = a * mu + b
        if mu >= 0 and 0 <= nu <= 4 * mu:
            vertices.add((mu, nu))
    return tuple(sorted(vertices))


def _le_conditional_mean_data(profile):
    """Exact finite vertices and recession slopes for the conditional mean gap."""
    profile = tuple(map(index, profile))
    physical = live_edge_score(profile, 0, 4)
    h, _, volume, births, _ = profile
    G, outside = physical['G'], physical['outside']
    marginal = tuple(Fraction(comb(4, j) * _le_falling(h, j)
                               * _le_falling(LE_M - h, 4 - j), _le_falling(LE_M, 4))
                     for j in range(5))
    same_edges = LE_REGIONS * (LE_UPDATES // LE_REGIONS - 1) * LE_WIDTH
    cross_edges = (LE_REGIONS - 1) * LE_WIDTH
    weights = []
    for j in range(5):
        for k in range(j, 5):
            same = Fraction(comb(4, j) * comb(4, k) * _le_falling(h, j + k)
                            * _le_falling(LE_M - h, 8 - j - k), _le_falling(LE_M, 8))
            weight = same_edges * same + cross_edges * marginal[j] * marginal[k]
            weights.append((j, k, weight * (1 if j == k else 2)))
    assert sum(weight for _, _, weight in weights) == LE_WIDTH * (LE_UPDATES - 1)
    base = volume - LE_WIDTH - births
    coef_mu = G - outside - LE_WIDTH // 2
    coef_nu = LE_WIDTH // 8 - G // 4
    values = []
    for mu, nu in _le_mean_vertices():
        mean = sum(weight * max(Fraction(0), (k * mu - nu) / 8,
                                 1 + (j + k) * mu / 8 - nu / 4)
                   for j, k, weight in weights)
        values.append((mu, nu, mean - base - coef_mu * mu - coef_nu * nu))
    # Recession g=( (j-r)_+ +(k-r)_+ )/8. Only integer r are breakpoints.
    rays = tuple((r, LE_GROUPS * (LE_UPDATES - 1)
                   * sum(marginal[j] * max(0, j - r) for j in range(5))
                   - coef_mu - coef_nu * r) for r in range(5))
    return tuple(values), rays


def live_edge_mean_obstruction(profile=LE_PROFILES[1]):
    """Exact Jensen obstruction for ALL mu>=0,0<=r<=4, beta>1.

    This diagnoses the edge relaxation, not the code. If the exact
    conditional mean of its setup statistic always exceeds the required
    score, even its exact Chernoff bound is at least1. A negative recession
    is reported explicitly and must not be replaced by a finite minimum.
    """
    values, rays = _le_conditional_mean_data(profile)
    minimum = min(values, key=lambda item: item[2])
    minimum_ray = min(rays, key=lambda item: item[1])
    nonnegative_rays = minimum_ray[1] >= 0
    obstruction = nonnegative_rays and minimum[2] >= 0
    return dict(profile=tuple(profile), vertex_count=len(values),
                minimizing_mu=str(minimum[0]), minimizing_nu=str(minimum[1]),
                minimizing_r=None if not minimum[0] else str(minimum[1] / minimum[0]),
                minimum_vertex_gap_exact=str(minimum[2]), minimum_vertex_gap=float(minimum[2]),
                recession_gaps_exact=[(r, str(gap)) for r, gap in rays],
                minimum_recession_gap=float(minimum_ray[1]),
                all_recession_slopes_nonnegative=nonnegative_rays,
                global_gap_unbounded_below=not nonnegative_rays,
                all_continuous_mu_r_obstructed=obstruction,
                domain='mu>=0,0<=r<=4,beta>1',
                exact_integer_rational_diagnostic=True,
                construction_counterexample=False, full_code_certificate=False)


def live_edge_float_bound(profile=LE_PROFILES[0], r=Fraction(3, 2), mu=8,
                          log_beta=1/8, p=None):
    """One small FLOAT evaluation; no marker search or numerical certificate."""
    import numpy as np
    profile = tuple(map(index, profile))
    h = profile[0]
    p = (7 / 8) * (h / LE_M) if p is None else float(p)
    log_beta = float(log_beta)
    if not (isfinite(p) and 0 < p < 1 and isfinite(log_beta) and log_beta > 0):
        raise ValueError('Require0<p<1 and finite positive log2(beta)')
    exponents = np.asarray(live_edge_exponents(r, mu), dtype=np.float64)
    score = live_edge_score(profile, r, mu)
    probabilities = np.asarray([comb(4, j) * p ** j * (1 - p) ** (4 - j)
                                for j in range(5)], dtype=np.float64)
    with np.errstate(over='raise', under='raise', invalid='raise'):
        matrix = np.exp2(log_beta * exponents) * probabilities[np.newaxis, :]
    if not np.all(np.isfinite(matrix)) or not np.all(matrix > 0):
        raise ArithmeticError('A positive finite five-state matrix is required')
    eigenvalues, eigenvectors = np.linalg.eig(matrix)
    principal = int(np.argmax(eigenvalues.real))
    raw_vector = eigenvectors[:, principal]
    if np.max(np.abs(raw_vector.imag)) > 1e-10:
        raise ArithmeticError('The computed Perron vector is not real')
    vector = raw_vector.real.copy()
    if vector.sum() < 0:
        vector *= -1
    if not np.all(np.isfinite(vector)) or not np.all(vector > 0):
        raise ArithmeticError('The computed Perron vector is not strictly positive')
    vector /= vector.max()
    ratios = (matrix @ vector) / vector
    if not np.all(np.isfinite(ratios)) or not np.all(ratios > 0):
        raise ArithmeticError('Invalid Collatz ratios')
    lam = float(np.nextafter(ratios.max(), np.inf))
    endpoint = LE_GROUPS * log2(float(probabilities @ vector) / float(vector.min()))
    bulk = LE_GROUPS * (LE_UPDATES - 1) * log2(lam)
    log_choose = (lgamma(LE_M + 1) - lgamma(h + 1) - lgamma(LE_M - h + 1)) / log(2)
    log_quota = log_choose + h * log2(p) + (LE_M - h) * log2(1 - p)
    field_free_prefactor = log_choose - LE_REGIONS * log_quota
    raw_bound = bulk + endpoint + field_free_prefactor - log_beta * float(score['score'])
    if not all(isfinite(value) for value in (bulk, endpoint, field_free_prefactor, raw_bound)):
        raise ArithmeticError('Nonfinite FLOAT result')
    return dict(profile=profile, statistics='actual_state_support',
                r=str(Fraction(str(r))), mu=mu, log2_beta=log_beta, p=p,
                selector_statistic=str(score['z']), required_score=str(score['score']),
                last_layer_boundary=str(score['boundary']),
                collatz_lambda=lam, minimum_ratio=float(ratios.min()),
                maximum_ratio=float(ratios.max()),
                relative_ratio_spread=float((ratios.max() - ratios.min()) / ratios.max()),
                positive_vector=vector.tolist(), log2_bulk=bulk,
                log2_endpoint=endpoint, log2_H_union=log_choose,
                log2_single_region_quota=log_quota,
                log2_H_union_and_quota=field_free_prefactor,
                raw_log2_probability_bound=raw_bound,
                capped_log2_probability_bound=min(0.0, raw_bound),
                conditional_on_every_partition_sequence=True,
                setup_statistic_key=(h, str(Fraction(str(r))), mu),
                support_union_used=False, floating_diagnostic_only=True,
                full_code_certificate=False)


def _le_probabilities(p):
    return tuple(comb(4, j) * p ** j * (1 - p) ** (4 - j) for j in range(5))


def _le_edge_counts(old_partition, new_partition):
    return tuple(tuple((old & new).bit_count() for new in new_partition)
                 for old in old_partition)


def _le_exact_edge_mgf(edges, exponents, root_beta, p):
    """Count-state enumeration for beta=root_beta^4, so all powers are integers."""
    groups = len(edges)
    probabilities = _le_probabilities(p)
    total = Fraction(0)
    for labels in product(range(5), repeat=2 * groups):
        old, new = labels[:groups], labels[groups:]
        probability = Fraction(1)
        for label in labels:
            probability *= probabilities[label]
        exponent = sum(edges[i][j] * exponents[old[i]][new[j]]
                       for i in range(groups) for j in range(groups))
        total += probability * root_beta ** exponent
    return total


def _le_physical_edge_mgf(width, old_partition, new_partition, exponents, root_beta, p):
    """Independent eligible-slot enumeration; integer histogram before weighting."""
    histogram = Counter()
    edges = _le_edge_counts(old_partition, new_partition)
    populations_old = [tuple((mask & block).bit_count() for block in old_partition)
                       for mask in range(1 << width)]
    populations_new = [tuple((mask & block).bit_count() for block in new_partition)
                       for mask in range(1 << width)]
    for old_mask, old in enumerate(populations_old):
        for new_mask, new in enumerate(populations_new):
            exponent = sum(edges[i][j] * exponents[old[i]][new[j]]
                           for i in range(len(edges)) for j in range(len(edges)))
            histogram[(old_mask.bit_count() + new_mask.bit_count(), exponent)] += 1
    return sum(count * p ** eligible * (1 - p) ** (2 * width - eligible)
               * root_beta ** exponent for (eligible, exponent), count in histogram.items())


def live_edge_checks():
    local_checks = holder_checks = edge_checks = boundary_checks = 0
    settings = ((Fraction(1), 4), (Fraction(3, 2), 4), (Fraction(2), 8))
    p, root_beta = Fraction(1, 3), Fraction(3, 2)
    probabilities = _le_probabilities(p)
    vector = tuple(Fraction(j + 1) for j in range(5))
    for r, mu in settings:
        exponents = live_edge_exponents(r, mu)
        for k, j in product(range(5), repeat=2):
            exact = max(Fraction(a * b) + Fraction(mu, 8) * (a * (k - r) + b * (j - r))
                        for a, b in product((0, 1), repeat=2))
            assert 4 * exact == exponents[k][j]
            local_checks += 1
        weighted_rows = [sum(probabilities[j] * vector[j] * root_beta ** (4 * exponents[k][j])
                             for j in range(5)) for k in range(5)]
        for incoming in combinations_with_replacement(range(5), 4):
            lhs = sum(probabilities[j] * vector[j]
                      * root_beta ** sum(exponents[k][j] for k in incoming) for j in range(5))
            rhs = Fraction(1)
            for k in incoming:
                rhs *= weighted_rows[k]
            assert lhs ** 4 <= rhs
            holder_checks += 1
        lam = max(weighted_rows[k] / vector[k] for k in range(5))
        initial = sum(probabilities[j] * vector[j] for j in range(5)) / min(vector)
        for width in (4, 8):
            old_partition = (15,) if width == 4 else (15, 240)
            partitions = [(15,)] if width == 4 else [
                (first, 255 ^ first) for rest in combinations(range(1, 8), 3)
                for first in (1 + sum(1 << i for i in rest),)]
            for new_partition in partitions:
                exact = _le_exact_edge_mgf(_le_edge_counts(old_partition, new_partition),
                                            exponents, root_beta, p)
                assert exact <= (lam * initial) ** (width // 4)
                edge_checks += 1
            physical_partition = (15,) if width == 4 else (85, 170)
            assert _le_physical_edge_mgf(width, old_partition, physical_partition,
                                         exponents, root_beta, p) == _le_exact_edge_mgf(
                                             _le_edge_counts(old_partition, physical_partition),
                                             exponents, root_beta, p)
            edge_checks += 1
        # Three one-quartet layers: first selector fixed empty, last is free.
        for labels in product(range(5), repeat=3):
            potential = sum(exponents[labels[t]][labels[t + 1]] for t in range(2))
            for second, final in product((0, 1), repeat=2):
                selectors = (0, second, final)
                overlap = 4 * second * final
                z = sum(selectors[t] * (labels[t] - r) for t in range(3))
                boundary = Fraction(mu, 2) * (4 - r)
                assert overlap + mu * z - boundary <= potential
                boundary_checks += 1
    # Omitting the final edge correction is genuinely invalid.
    exponents = live_edge_exponents(Fraction(3, 2), 8)
    potential = exponents[0][0] + exponents[0][4]
    uncorrected = 8 * (4 - Fraction(3, 2))
    assert uncorrected > potential
    assert uncorrected - 4 * (4 - Fraction(3, 2)) == potential
    # Three full actual states: L-f exceeds physical adjacent overlap by width.
    assert 12 - 4 > 4 and 12 - 4 - 4 == 4
    # A constant positive row probability has exact Collatz value1 at beta1.
    assert sum(probabilities) == 1
    for r, mu in settings:
        score = live_edge_score(LE_PROFILES[0], r, mu)
        assert score['score'].denominator == 1
    assert len(_le_mean_vertices()) == 21
    first_values, first_rays = _le_conditional_mean_data(LE_PROFILES[0])
    second_values, second_rays = _le_conditional_mean_data(LE_PROFILES[1])
    _, third_rays = _le_conditional_mean_data(LE_PROFILES[2])
    first_min = min(first_values, key=lambda item: item[2])
    second_min = min(second_values, key=lambda item: item[2])
    assert first_min[:2] == second_min[:2] == (Fraction(4), Fraction(8))
    assert -947 < first_min[2] < -946 and all(gap >= 0 for _, gap in first_rays)
    assert 36014 < second_min[2] < 36015
    assert 15628 < min(gap for _, gap in second_rays) < 15629
    assert any(gap < 0 for _, gap in third_rays)
    try:
        live_edge_exponents(Fraction(3, 2), 6)
    except ValueError:
        pass
    else:
        raise AssertionError('Half-integer exponent validation failed')
    print(dict(exact_local_maxima=local_checks, exact_Holder_multisets=holder_checks,
               exact_width4_width8_edge_checks=edge_checks,
               exact_boundary_checks=boundary_checks,
               missing_boundary_negative_controls=True,
               half_integer_r_checked=True, exact_mean_vertices=21,
               exact_mean_profile_controls=3, full_code_certificate=False), flush=True)


def live_edge_main():
    parser = argparse.ArgumentParser(description=__doc__)
    mode = parser.add_mutually_exclusive_group(required=True)
    mode.add_argument('--checks', action='store_true')
    mode.add_argument('--probe', action='store_true')
    mode.add_argument('--mean-obstruction', action='store_true')
    parser.add_argument('--profile', type=int, nargs=5, default=LE_PROFILES[0],
                        metavar=('h', 'R', 'Ls', 'f', 'F'))
    parser.add_argument('--r', type=Fraction, default=Fraction(3, 2))
    parser.add_argument('--mu', type=int, default=8)
    parser.add_argument('--log-beta', type=float, default=1/8)
    parser.add_argument('--p', type=float)
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
            live_edge_checks()
        elif args.mean_obstruction:
            print(live_edge_mean_obstruction(args.profile), flush=True)
        else:
            print(live_edge_float_bound(args.profile, args.r, args.mu, args.log_beta, args.p), flush=True)
        print(dict(seconds=round(time.monotonic() - started, 3)), flush=True)
    finally:
        timer.cancel()


if __name__ == '__main__':
    live_edge_main()
