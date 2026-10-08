"""Exact regional coefficients for the randomized-outer quiet bound.

Only a bounded list of occupancies is investigated. Float arithmetic selects
rational witnesses; Arb evaluates the positive polynomial coefficients. No
native code, timings, or whole-code distance claim is involved.
"""

from fractions import Fraction as F
from functools import lru_cache
from math import comb, exp2, log2
import argparse
import json

from hd_randomized_outer_joint_quiet import (
    directed_interval as directed_marker_interval,
    joint_moment_enclosures, joint_moment_logs,
)
from hd_randomized_outer_proof import Parameters


@lru_cache(None)
def moment_logs(mark, field=65536):
    return joint_moment_logs(F(mark), field)


def weighted_logsum_mean(log_coefficients, log_odds):
    terms = tuple(value+j*log_odds for j, value in enumerate(log_coefficients))
    maximum = max(terms)
    weights = tuple(exp2(value-maximum) for value in terms)
    total = sum(weights)
    return maximum+log2(total), sum(j*w for j, w in enumerate(weights))/total


def optimized_marker(active_rows, mark, alpha, parameters=Parameters()):
    """Float saddle point; optimize the marker without a coarse p-grid.

    The free final factor is included in the stationary equation. The result
    is only a witness proposal and does not certify a bound or an optimum.
    """
    p = parameters
    if not 0 < active_rows < p.rows or not 0 < alpha <= 1:
        raise ValueError('Need interior occupancy and 0<alpha<=1')
    quiet, free = moment_logs(F(mark), p.field)
    log_choose = tuple(log2(comb(p.width, j)) for j in range(p.width+1))
    a = float(alpha)
    ordinary = tuple(c+a*k for c, k in zip(log_choose, quiet))
    terminal = tuple(c+a*k for c, k in zip(log_choose, free))
    lower, upper = -128.0, 128.0
    for unused in range(80):
        middle = (lower+upper)/2
        _, mean = weighted_logsum_mean(ordinary, middle)
        _, final_mean = weighted_logsum_mean(terminal, middle)
        if (p.packets-1)*mean+final_mean < p.outer_n*active_rows:
            lower = middle
        else:
            upper = middle
    log_odds = (lower+upper)/2
    log_ordinary, _ = weighted_logsum_mean(ordinary, log_odds)
    log_terminal, _ = weighted_logsum_mean(terminal, log_odds)
    choose = log2(comb(p.rows, active_rows))
    log_beta = p.outer_n*log2(p.field)-(p.outer_n-p.outer_k)*log2(p.field-1)
    bound = (a*(active_rows*log_beta-log2(p.field-1)-p.limit*log2(float(mark)))
             +(p.packets-1)*log_ordinary+log_terminal
             -p.outer_n*active_rows*log_odds-(p.outer_n-1)*choose)
    marker = F(1/(1+exp2(-log_odds))).limit_denominator(1 << 40)
    if not 0 < marker < 1:
        raise ArithmeticError('Saddle marker rounded to a boundary')
    return dict(active_rows=active_rows, mark=str(F(mark)), alpha=str(F(alpha)),
                marker=str(marker), log_odds=log_odds, marker_screen_log2=bound)


def screen_selected(occupancies=(954, 1024, 1536, 2048, 4096, 5120, 5939),
                    parameters=Parameters()):
    """Bounded mark/exponent grid with an optimized marker for every entry."""
    p = parameters
    marks = tuple(F(1, 1 << shift) for shift in (0, 1, 2, 4, 6, 8, 10, 12, 14, 16))
    marks += (F(2, 3*(p.field-1)), F(133, 128*(p.field-1)))
    alphas = tuple(F(a, 32) for a in (1, 2, 3, 4, 6, 8, 12, 16, 24, 32))
    best = {r: None for r in occupancies}
    for mark in marks:
        for alpha in alphas:
            for active in occupancies:
                proposal = optimized_marker(active, mark, alpha, p)
                if best[active] is None or proposal['marker_screen_log2'] < best[active]['marker_screen_log2']:
                    best[active] = proposal
    return tuple(best[r] for r in occupancies)


def _arb_fraction(value):
    from flint import arb
    value = F(value)
    return arb(value.numerator)/value.denominator


def regional_coefficients(quiet, free, alpha, marker, packets, active):
    """Return directed fixed-count moments and normalized coefficient masses.

    quiet/free contain the local j-slot moments for j=0,...,width. The last
    packet is free. The marker only scales coefficients; exact returned
    moments are independent of it. All polynomial coefficients are positive.
    """
    from flint import arb, arb_poly
    if len(quiet) != len(free) or len(quiet) < 2 or packets < 1:
        raise ValueError('Invalid local polynomial or packet count')
    width = len(quiet)-1
    rows = width*packets
    if not 0 <= active <= rows or not 0 < alpha <= 1 or not 0 < marker < 1:
        raise ValueError('Invalid occupancy, exponent, or scaling marker')
    a, m = _arb_fraction(alpha), _arb_fraction(marker)
    convert = lambda value: _arb_fraction(value) if isinstance(value, F) else arb(value)
    coefficients = tuple(comb(width, j)*m**j*(1-m)**(width-j)*convert(k)**a
                         for j, k in enumerate(quiet))
    last_coefficients = tuple(comb(width, j)*m**j*(1-m)**(width-j)*convert(k)**a
                              for j, k in enumerate(free))
    ordinary = sum(coefficients, arb(0))
    terminal = sum(last_coefficients, arb(0))
    polynomial = arb_poly([value/ordinary for value in coefficients])
    last_polynomial = arb_poly([value/terminal for value in last_coefficients])
    prefix = polynomial**(packets-1)
    quiet_probability = (prefix*polynomial)[active]
    last_probability = (prefix*last_polynomial)[active]
    if not quiet_probability > 0 or not last_probability > 0:
        raise ArithmeticError('Coefficient enclosure lost positivity; increase precision or improve marker')
    binomial_mass = arb(comb(rows, active))*m**active*(1-m)**(rows-active)
    return dict(quiet=quiet_probability*ordinary**packets/binomial_mass,
                last=last_probability*ordinary**(packets-1)*terminal/binomial_mass,
                quiet_probability=quiet_probability, last_probability=last_probability,
                ordinary=ordinary, terminal=terminal, binomial_mass=binomial_mass)


def directed_endpoint(active_rows, mark, alpha, marker=None, parameters=Parameters()):
    """Exact regional composition for one occupancy, with Arb enclosures."""
    from flint import arb
    p = parameters
    mark, alpha = F(mark), F(alpha)
    if p.width != 64 or p.rows % p.width:
        raise ValueError('Current local moments require width 64 and whole regional packets')
    if not 0 < active_rows < p.rows:
        raise ValueError('Need an interior active-row count')
    proposal = optimized_marker(active_rows, mark, alpha, p)
    marker = F(proposal['marker']) if marker is None else F(marker)
    quiet, free = joint_moment_enclosures(mark, p.field)
    region = regional_coefficients(quiet, free, alpha, marker, p.rows//p.width, active_rows)
    logarb = lambda value: value.log()/arb(2).log()
    choose = logarb(arb(comb(p.rows, active_rows)))
    log_beta = p.outer_n*logarb(arb(p.field))-(p.outer_n-p.outer_k)*logarb(arb(p.field-1))
    prefix = (choose+_arb_fraction(alpha)*(active_rows*log_beta-logarb(arb(p.field-1))
                                          -p.limit*logarb(_arb_fraction(mark))))
    exact = prefix+(p.outer_n-1)*logarb(region['quiet'])+logarb(region['last'])
    marker_bound = (prefix+(p.packets-1)*logarb(region['ordinary'])
                    +logarb(region['terminal'])-p.outer_n*logarb(region['binomial_mass']))
    improvement = marker_bound-exact
    return dict(active_rows=active_rows, mark=str(mark), alpha=str(alpha), marker=str(marker),
                exact_coefficient_log2_upper=str(exact), marker_log2_upper=str(marker_bound),
                exact_coefficient_gain_bits=str(improvement),
                quiet_coefficient_mass=str(region['quiet_probability']),
                last_coefficient_mass=str(region['last_probability']),
                below_2_to_minus_40=bool(exact < -40), directed_acceptance=True,
                covers_only_this_occupancy=True, whole_code_certificate=False)


def rank_refinement_floor(receipt, parameters=Parameters()):
    """Limitation at this witness, not a construction or global impossibility.

    The current residual rank rho is exact at rho=0,1, and actual rank is
    at most four. Thus actual_rank-rho<=2 termwise in the positive expansion.
    At z=1 the existing syndrome-rank moment is already exact.
    """
    from flint import arb
    p = parameters
    alpha = _arb_fraction(F(receipt['alpha']))
    gain = (arb(0) if F(receipt['mark']) == 1 else
            2*alpha*(p.packets-1)*arb(p.field).log()/arb(2).log())
    floor = arb(receipt['exact_coefficient_log2_upper'])-gain
    return dict(active_rows=receipt['active_rows'], mark=receipt['mark'],
                alpha=receipt['alpha'], maximum_rank_only_gain_bits=str(gain),
                best_rank_only_log2_lower=str(floor),
                same_witness_cannot_certify=bool(floor > -40),
                fixed_witness_only=True, construction_obstruction=False)


def directed_sparse_extension(first=954, last=962, alpha=F(3, 32)):
    """Cover a short whole interval from one exact-coefficient evaluation.

    At z=1 each local factor q^(-alpha min(j,4)) is nonincreasing in j;
    the final free factor is one. Nested uniform subsets couple successive
    regional occupancies, so their entire routing moment is nonincreasing.
    Only the increasing explicit outer prefactor is charged at the last r.
    """
    from flint import arb
    p = Parameters()
    if not 1 <= first <= last <= p.rows//2:
        raise ValueError('Need an increasing-binomial occupancy interval')
    endpoint = directed_endpoint(first, F(1), alpha)
    logarb = lambda value: value.log()/arb(2).log()
    log_beta = p.outer_n*logarb(arb(p.field))-(p.outer_n-p.outer_k)*logarb(arb(p.field-1))
    maximum = (arb(endpoint['exact_coefficient_log2_upper'])
               +logarb(arb(comb(p.rows, last))/comb(p.rows, first))
               +_arb_fraction(alpha)*(last-first)*log_beta)
    union = maximum+logarb(arb(last-first+1))
    return dict(class_kind='all_nonzero_quiet_words', active_rows_first=first,
                active_rows_last=last, mark='1', alpha=str(alpha),
                first_endpoint_log2_upper=endpoint['exact_coefficient_log2_upper'],
                maximum_prefactor_envelope_log2_upper=str(maximum),
                union_log2_upper=str(union), below_2_to_minus_40=bool(union < -40),
                interior_justification='nested-subset monotonicity at z=1',
                exact_regional_occupancies_evaluated=1, whole_code_certificate=False)


def directed_dense_extension():
    """Four exact points plus a convex marker interval cover 5865..8192."""
    from flint import arb
    mark, alpha = F(1, 16384), F(3, 8)
    points = tuple(directed_endpoint(r, mark, alpha) for r in range(5865, 5869))
    interval = directed_marker_interval(5869, 8192, mark,
                                        F(12860472643, 12983074412), alpha)
    probability = arb(2)**arb(interval['union_log2_upper'])
    for point in points:
        probability += arb(2)**arb(point['exact_coefficient_log2_upper'])
    union = probability.log()/arb(2).log()
    return dict(class_kind='quiet_low_weight', active_rows_first=5865,
                active_rows_last=8192, union_log2_upper=str(union),
                below_2_to_minus_40=bool(union < -40), exact_points=points,
                marker_interval=interval, whole_code_certificate=False)


def main():
    from flint import ctx
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--rows', nargs='+', type=int,
                        default=[954, 1024, 1536, 2048, 4096, 5120, 5939])
    parser.add_argument('--precision', type=int, default=192)
    parser.add_argument('--extensions-only', action='store_true')
    args = parser.parse_args()
    ctx.prec = args.precision
    if not args.extensions_only:
        for proposal in screen_selected(tuple(args.rows)):
            receipt = directed_endpoint(proposal['active_rows'], proposal['mark'],
                                        proposal['alpha'], proposal['marker'])
            print(json.dumps(receipt), flush=True)
            print(json.dumps(rank_refinement_floor(receipt)), flush=True)
    print(json.dumps(directed_sparse_extension()), flush=True)
    print(json.dumps(directed_dense_extension()), flush=True)


if __name__ == '__main__':
    main()
