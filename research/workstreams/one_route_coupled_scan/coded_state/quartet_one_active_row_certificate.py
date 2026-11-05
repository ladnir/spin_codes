"""Exact whole-class certificate for one active outer RS row.

This counts true nonzero input rows, not footprint-eligible rows.  Inject
quartet_selector_capacity.py and quartet_adaptive_zero_budget.py for the
already proved directed h=1 joint-zero threshold.  The fixed-message
zero-budget proof applies to one prescribed input row regardless of how
many rows the resulting footprint makes eligible.

Thirteen exact partition-growth inequalities and128 exact budget
comparisons prove that any bad one-row message lies in one of385
families, each with probability below2^-92.  The conservative whole-class
budget is2^-82.  No local state DP, matrix witness, numerical optimization,
state-support union or construction change is used.  Parent-only jobs.
"""

import argparse
from math import comb
from operator import index
import os
import shutil
import threading
import time


OR_N = 1 << 21
OR_T = 4096
OR_WIDTH = 512
OR_W = 943718
OR_GROWTH_BITS = 20
OR_FAMILY_BITS = 92


def _or_ceildiv(a, b):
    return -((-a) // b)


def one_row_level(u):
    u = index(u)
    if not 0 <= u <= OR_WIDTH:
        raise ValueError('State population must lie in0..512')
    if u == 0:
        return 0
    if u <= 256:
        return u.bit_length()
    if u == 512:
        return 18
    holes = 512 - u
    return 17 - (holes - 1).bit_length()


def _or_floor_log_ratio(numerator, denominator):
    """Floor log2(numerator/denominator), using integers only."""
    exponent = numerator.bit_length() - denominator.bit_length()
    if exponent >= 0:
        if denominator << exponent > numerator:
            exponent -= 1
    elif numerator << (-exponent) < denominator:
        exponent -= 1
    return exponent


def one_row_growth_certificate():
    """Uniform probability<=2^-20 of no level gain at a partial state."""
    comparisons = []
    for u in (4, 8, 16, 32, 64, 128):
        k = u // 2 - 1
        numerator = comb(128, k) * comb(4 * k, u)
        denominator = comb(512, u)
        if (numerator << OR_GROWTH_BITS) > denominator:
            raise ArithmeticError('Low-population growth comparison failed')
        comparisons.append(dict(kind='support', endpoint=u, quartets=k,
                                inverse_probability_floor_bits=_or_floor_log_ratio(denominator, numerator)))
    for holes in (4, 8, 16, 32, 64, 128, 256):
        k = holes // 8 + 1
        numerator = comb(128, k) * comb(512 - 4 * k, holes - 4 * k)
        denominator = comb(512, holes)
        if (numerator << OR_GROWTH_BITS) > denominator:
            raise ArithmeticError('High-population growth comparison failed')
        comparisons.append(dict(kind='holes', endpoint=holes, quartets=k,
                                inverse_probability_floor_bits=_or_floor_log_ratio(denominator, numerator)))
    return comparisons


def one_row_budget_row(d):
    """Three fixed-message tails and their deterministic contradiction."""
    d = index(d)
    if not 1 <= d <= 128:
        raise ValueError('Actual row weight is128+d, with1<=d<=128')
    cost = 32 * d + 362
    state = _or_ceildiv(cost, 14)
    extinctions = _or_ceildiv(cost, 116)
    partitions = max(48, _or_ceildiv(cost, 12))
    assert state >= 24 and 3 * OR_N <= state * (1 << 18)
    assert partitions >= 48 and 3 * OR_T <= partitions * (1 << 8)
    base = 270 + 32 * d  # gamma<2,8192 row choices,2^256 supports,q^d messages
    exponents = (base - 14 * state, base - 116 * extinctions,
                 base - 12 * partitions)
    assert max(exponents) <= -OR_FAMILY_BITS
    maximum = state - 1 + partitions - 1 + 33 * (extinctions - 1)
    required = 170 + 16 * d
    if maximum >= required:
        raise ArithmeticError('The one-row occupation contradiction failed')
    return dict(d=d, row_weight=128 + d, state_zero_threshold=state,
                extinction_threshold=extinctions, bad_partition_threshold=partitions,
                three_tail_exponents=exponents,
                complement_score_upper=maximum, required_score_lower=required)


def one_row_certificate():
    if 'adaptive_zero_threshold' not in globals():
        raise RuntimeError('Inject the adaptive-zero and selector-log helpers')
    joint = adaptive_zero_threshold(1, 'joint')
    if joint['threshold'] > 254 or joint['exact_log2_tail_upper'] > -92:
        raise ArithmeticError('The one-input-row joint-zero budget was not certified')
    assert (OR_W + 253) // 512 == 1843
    assert 1843 + 18 + 1 == 1862
    assert 16 * 128 - 16 - 1862 == 170
    growth = one_row_growth_certificate()
    budgets = [one_row_budget_row(d) for d in range(1, 129)]
    assert 1 + 3 * len(budgets) == 385 < 1 << 10
    # Each family is strictly below2^-92, including the accepted density.
    assert 10 - OR_FAMILY_BITS == -82
    return dict(true_active_input_rows=1, footprint_eligible_rows_unrestricted=True,
                exact_growth_comparisons=len(growth), exact_weight_budget_comparisons=len(budgets),
                growth_probability_exponent=-OR_GROWTH_BITS,
                joint_bad_threshold=joint['threshold'], full_incoming_updates_upper=1843,
                family_count=385, per_family_exponent=-92,
                conservative_global_failure_exponent=-82,
                minimum_integer_contradiction_margin=min(row['required_score_lower'] - row['complement_score_upper']
                                                         for row in budgets),
                growth_comparisons=growth, exact_integer_certificate=True,
                no_state_support_union=True, all_input_restarts_retained=True,
                full_code_certificate=False)


def one_row_checks():
    level_checks = coverage_checks = occupation_checks = 0
    levels = [one_row_level(u) for u in range(513)]
    assert levels[0] == 0 and levels[-1] == 18
    for u in range(512):
        assert 0 <= levels[u + 1] - levels[u] <= 1
        level_checks += 1
    # Check that the thirteen endpoint tails contain every bad level event.
    for u in range(1, 512):
        for covered in range(4 * ((u + 3) // 4), 513, 4):
            bad = levels[covered] == levels[u]
            if u < 256:
                endpoint = 1 << (u.bit_length() - 1)
                if endpoint <= 2:
                    assert not bad
                else:
                    assert bad == (covered <= 2 * endpoint - 4)
            else:
                holes = 512 - u
                endpoint = 1 << (holes - 1).bit_length()
                if endpoint <= 2:
                    assert not bad
                else:
                    k = endpoint // 8 + 1
                    assert bad == (512 - covered >= 4 * k)
            coverage_checks += 1
    # Independent binary-run checks of the occupied-region inequality.
    for updates, regional_updates in ((4, 2), (6, 2), (8, 2), (8, 4), (10, 2)):
        for mask in range(1 << updates):
            bits = [(mask >> t) & 1 for t in range(updates)]
            active = sum(bits)
            deaths = sum(bits[t - 1] and not bits[t] for t in range(1, updates))
            covered = [bool(bits[t] or (bits[t - 1] if t else 0)) for t in range(updates)]
            assert sum(covered) == active + deaths
            occupied = sum(any(covered[start:start + regional_updates])
                           for start in range(0, updates, regional_updates))
            assert occupied <= (active + deaths) // regional_updates + 2 * deaths + 1
            occupation_checks += 1
    certificate = one_row_certificate()
    print(dict(exact_level_lipschitz_checks=level_checks,
               exact_growth_tail_coverage_checks=coverage_checks,
               exact_binary_region_occupation_checks=occupation_checks,
               **certificate), flush=True)


def one_row_main():
    parser = argparse.ArgumentParser(description=__doc__)
    mode = parser.add_mutually_exclusive_group(required=True)
    mode.add_argument('--checks', action='store_true')
    mode.add_argument('--certify', action='store_true')
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
            one_row_checks()
        else:
            print(one_row_certificate(), flush=True)
        print(dict(seconds=round(time.monotonic() - started, 3)), flush=True)
    finally:
        timer.cancel()


if __name__ == '__main__':
    one_row_main()
