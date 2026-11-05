"""FLOAT screen for a sound product comparison of XOR2 outer columns.

The exact column-support kernel is dominated by an iid coordinate law:
with probability p a coordinate is uniform GF256, otherwise it is zero.
The maximum density ratio is paid once per column. This removes column-mask
correlations by domination, not an independence assumption. These floating
evaluations are not directed receipts or a whole-code certificate.
"""
from functools import lru_cache
from math import comb, log2
from time import perf_counter

import numpy as np
from scipy.optimize import minimize, minimize_scalar

from hd_bucket_feedback import profiles
from hd_common_refresh_screen import local_return_cap
from hd_grouped_quiet import quiet_screen
from hd_grouped_subfield_screen import bound, density_bits
from hd_xor_column_probe import (column_support_counts, mds_prescribed_support,
                                 uniform_column_support_counts)
from hd_randomized_outer_state_profile import coefficient_one, multiply, power, total

Q = 256
GROUP = 8
COLUMNS = 32
BASELINE_BITS = 256*log2(Q)-128*log2(Q-1)


@lru_cache(None)
def comparison_masks(active):
    counts = uniform_column_support_counts(Q, GROUP, (1 << active)-1)
    return tuple((mask.bit_count(), log2(count)-active*log2(Q))
                 for mask, count in enumerate(counts) if count)


@lru_cache(None)
def marked_comparison_masks(active):
    by_size = [[0]*(1 << GROUP) for _ in range(active+1)]
    for source in range(1 << active):
        counts = column_support_counts(Q, GROUP, source)
        by_size[source.bit_count()] = [a+b for a, b in zip(by_size[source.bit_count()], counts)]
    return tuple((target.bit_count(), tuple((size, log2(by_size[size][target])-active*log2(Q))
                                            for size in range(active+1) if by_size[size][target]))
                 for target in range(1 << GROUP) if any(row[target] for row in by_size))


def column_density_bits(active, probability, input_log_mark=0.):
    if not 1 <= active <= GROUP or not 0 < probability <= 1:
        raise ValueError('Invalid active count or probability')
    zero = 1-probability+probability/Q
    nonzero = probability/Q
    entries = (comparison_masks(active) if input_log_mark == 0 else
               tuple((weight, total(value+size*input_log_mark for size, value in terms))
                     for weight, terms in marked_comparison_masks(active)))
    return max(logmass-weight*log2(Q-1)-weight*log2(nonzero)
               -(GROUP-weight)*log2(zero)
               for weight, logmass in entries)


def group_density_bits(active, probability, include_row_choices=True, input_log_mark=0.):
    choices = log2(comb(GROUP, active)) if include_row_choices else 0.
    return (32*active*log2(Q)-16*active*log2(Q-1)
            +COLUMNS*column_density_bits(active, probability, input_log_mark)+choices
            -17*active*input_log_mark)


def all_row_sets_density_bits(probability):
    values = [group_density_bits(active, probability) for active in range(1, 9)]
    peak = max(values)
    return peak+log2(sum(2**(value-peak) for value in values))


@lru_cache(None)
def profile_arrays():
    items = profiles(16, 64, 4)
    return (np.array([row for row, count in items]),
            np.array([count for row, count in items], dtype=float))


def thinned_matrices(z, probability):
    """Profile-wise average of the ordinary local operator before alpha."""
    q = float(Q)
    state, width, h, length, blocks = 16, 64, 4, 16, 4
    occ, counts = profile_arrays()
    degrees = occ.sum(axis=1)
    quiet = np.array([sum(comb(j, c)*z**(j-c)*(1-z)**c*q**(-c-min(h, j-c))
                          for c in range(j+1)) for j in range(length+1)])
    free = (z+(1-z)/q)**np.arange(length+1)
    base = q**(-h)
    excess = np.array([sum(comb(j, c)*z**(j-c)*(1-z)**c*(q**(-c)-base)
                           for c in range(min(h, j+1))) for j in range(length+1)])
    local_return = np.minimum(np.maximum(0, free-quiet),
        (1-q**(-np.minimum(h, np.arange(length+1)).astype(float)))
        *z**np.maximum(0, length-h+1-np.arange(length+1)))
    local_return = np.minimum(local_return,
        [local_return_cap(j, z, h, length, q) for j in range(length+1)])
    thinning = np.zeros((length+1, length+1))
    for j in range(length+1):
        for t in range(j+1):
            thinning[j, t] = comb(j, t)*probability**t*(1-probability)**(j-t)
    k_local = thinning @ quiet
    t_local = thinning @ free
    r_local = thinning @ local_return
    v_excess = thinning @ (free*excess[::-1])
    nonempty_thinning = thinning.copy()
    nonempty_thinning[:, 0] = 0
    some_base = base*(nonempty_thinning @ free)
    some_excess = nonempty_thinning @ (free*excess[::-1])
    k = np.prod(k_local[occ], axis=1)
    t = np.prod(t_local[occ], axis=1)
    running = np.ones(len(occ))
    ret = np.zeros(len(occ))
    live_running = np.ones(len(occ))
    live_nonzero = np.zeros(len(occ))
    joint = np.zeros((len(occ), 2, 2))
    joint[:, 0, 0] = 1
    for block in range(blocks):
        kb, rb = k_local[occ[:, block]], r_local[occ[:, block]]
        ret = ret*kb+running*rb
        running *= kb+rb
        vb, ve = base*t_local[occ[:, block]], v_excess[occ[:, block]]
        live_nonzero = live_nonzero*vb+live_running*ve
        live_running *= vb+ve
        empty_block = (1-probability)**occ[:, block]
        block_values = ((empty_block*base, empty_block*excess[-1]),
                        (some_base[occ[:, block]], some_excess[occ[:, block]]))
        updated = np.zeros_like(joint)
        for old_input in range(2):
            for old_state in range(2):
                for new_input in range(2):
                    for new_state in range(2):
                        updated[:, old_input | new_input, old_state | new_state] += (
                            joint[:, old_input, old_state]*block_values[new_input][new_state])
        joint = updated
    d = q**state-1
    live = live_nonzero/(1-q**(-state))
    ret = np.minimum(ret, d*live)
    density = 2**density_bits(Q, state)
    # Separate empty and nonempty thinning with a positive recurrence.
    empty_live, empty_running = 0., 1.
    for _ in range(blocks):
        empty_live = empty_live*base+empty_running*excess[-1]
        empty_running *= base+excess[-1]
    empty_live /= 1-q**(-state)
    inverse = min(density*empty_live, quiet_screen(z, q=Q, state=state, width=width, h=h))
    actual_live = density*joint[:, 1, 1]/(1-q**(-state)) + (1-probability)**degrees*inverse
    values = np.stack((k, density*np.maximum(0, t-k), ret/d,
                       np.maximum(0, actual_live)), axis=1).reshape((-1, 2, 2))
    return degrees, counts, values


def thinned_coefficients(z, alpha, probability):
    degrees, counts, values = thinned_matrices(z, probability)
    result = np.zeros((65, 2, 2))
    for row in range(2):
        for column in range(2):
            result[:, row, column] = np.bincount(degrees,
                weights=counts*values[:, row, column]**alpha, minlength=65)
    with np.errstate(divide='ignore'):
        return np.log2(result)


def evaluate(rows, active, z, alpha, probability, regional=False, input_log_mark=0.):
    coeff = thinned_coefficients(z, alpha, probability)
    result, odds = bound(rows, z, alpha, coeff, q=Q, regional=regional)
    density = (all_row_sets_density_bits(probability) if active == 0 else
               group_density_bits(active, probability, input_log_mark=input_log_mark))
    return result+alpha*rows*(density-BASELINE_BITS), odds


def proposal(rows, active, iterations=60, regional=False, optimize_probability=True, input_mark=False):
    start = perf_counter()
    best = None
    density = all_row_sets_density_bits if active == 0 else lambda p: group_density_bits(active, p)
    initial = minimize_scalar(density, bounds=(.05, .999999), method='bounded')
    p0 = float(initial.x)
    def objective(theta):
        nonlocal best
        ell, alpha = theta[:2]
        probability = 1/(1+2**(-theta[2])) if len(theta) > 2 else p0
        input_log_mark = theta[3] if len(theta) > 3 else 0.
        if not .0001 <= ell <= 16 or not .005 <= alpha <= 1 or probability < .01 or input_log_mark < 0:
            return 1e12
        z = 2**(-ell)
        try:
            value, odds = evaluate(rows, active, z, alpha, probability, regional, input_log_mark)
        except (ArithmeticError, ValueError, FloatingPointError):
            return 1e12
        if not np.isfinite(value):
            return 1e12
        if best is None or value < best['log2_upper_screen']:
            best = dict(active_groups=rows, active_rows_per_group=active,
                        all_row_set_mixtures=active == 0, mark=z, alpha=alpha,
                        thinning_probability=probability, input_log_mark=input_log_mark,
                        group_density_bits=(density(probability) if not input_mark else
                            group_density_bits(active, probability, input_log_mark=input_log_mark)),
                        marker_log_odds=odds, log2_upper_screen=value, regional=regional,
                        directed=False, whole_code_certificate=False)
        return value
    for alpha in ((.0625, .25, 1.) if rows <= 16 else (.0625,)):
        seed = [max(.002, rows/512), alpha]
        if optimize_probability:
            seed.append(log2(p0/(1-p0)))
        if input_mark:
            if not optimize_probability or active == 0:
                raise ValueError('Input mark currently requires one active-row class and variable probability')
            seed[2] = log2(.98/.02)
            seed.append(20.)
        minimize(objective, seed, method='Nelder-Mead', options={'maxiter': iterations})
    if best is None:
        raise ArithmeticError('No usable FLOAT proposal')
    best['elapsed_seconds'] = perf_counter()-start
    return best


@lru_cache(None)
def masks_for_active(active):
    return tuple((mask, tuple((target, log2(count)-active*log2(Q)
                              +target.bit_count()*log2(Q/(Q-1)))
                             for target, count in enumerate(uniform_column_support_counts(Q, GROUP, mask))
                             if count))
                 for mask in range(1, 1 << GROUP) if mask.bit_count() == active)


def one_group_region_words(z, nonzero=True):
    coeff = thinned_coefficients(z, 1., 1.)
    packet_zero = tuple(tuple(row) for row in coeff[0])
    packet_one = tuple(tuple(row-6.) for row in coeff[1])
    if nonzero:
        # A forced nonzero one-coordinate input has nonzero feedback.
        # Keep simple conditional caps for the two flat-state entries.
        kappa = log2(Q/(Q-1))
        packet_one = ((-float('inf'), density_bits(Q, 16)+log2(z)),
                      (packet_one[1][0]+kappa, packet_one[1][1]+kappa))
    region_zero = power(packet_zero, 128)
    region_one = tuple(tuple(value-7. for value in row)
                       for row in coefficient_one(packet_zero, packet_one, 128))
    identity = ((0., -float('inf')), (-float('inf'), 0.))
    words = [identity]
    for bit in range(GROUP):
        words = [multiply(words[mask & ((1 << bit)-1)],
                          region_one if mask >> bit & 1 else region_zero)
                 for mask in range(1 << (bit+1))]
    return words


def one_group_correlated_bound(active, z, alpha):
    """One active group; exact column masks, Jensen over routes per 8-region block.

    The mask mixture stays before alpha. Route averaging per block precedes
    alpha by Jensen; this may loosen the bound but preserves soundness.
    The source-row-set union is outside alpha in this diagnostic.
    """
    words = one_group_region_words(z)
    cases = []
    prefix = 13.+alpha*(32*active*log2(Q)-16*active*log2(Q-1)
                       -log2(Q-1)-(2**21*7//20)*log2(z))
    for mask, support in masks_for_active(active):
        block = tuple(tuple(total(logmass-target.bit_count()*log2(Q/(Q-1))
                                  +words[target][row][column]
                                  for target, logmass in support)
                            for column in range(2)) for row in range(2))
        fractional = tuple(tuple(alpha*value for value in row) for row in block)
        cases.append(prefix+total(power(fractional, 32)[0]))
    return total(cases)


def one_group_proposal(active, iterations=60):
    start = perf_counter()
    best = None
    def objective(theta):
        nonlocal best
        ell, alpha = theta
        if not .00001 <= ell <= .1 or not .005 <= alpha <= 1:
            return 1e12
        z = 2**(-ell)
        value = one_group_correlated_bound(active, z, alpha)
        if best is None or value < best['log2_upper_screen']:
            best = dict(active_groups=1, active_rows_per_group=active,
                        mark=z, alpha=alpha, log2_upper_screen=value,
                        exact_column_mask_mixture=True, route_block_jensen=True,
                        directed=False, whole_code_certificate=False)
        return value
    for alpha in (.25, 1.):
        minimize(objective, (.001, alpha), method='Nelder-Mead', options={'maxiter': iterations})
    best['elapsed_seconds'] = perf_counter()-start
    return best


@lru_cache(None)
def exact_mask_kernels(active):
    if active not in (1, 2):
        raise ValueError('Univariate/bivariate screen only')
    result = []
    for mask in range(1, 1 << GROUP):
        if mask.bit_count() != active:
            continue
        indices = [i for i in range(GROUP) if mask >> i & 1]
        kernels = []
        for small in range(1 << active):
            support = sum(1 << indices[i] for i in range(active) if small >> i & 1)
            counts = column_support_counts(Q, GROUP, support)
            entries = tuple((target, log2(count)-small.bit_count()*log2(Q-1)
                             +target.bit_count()*log2(Q/(Q-1)))
                            for target, count in enumerate(counts) if count)
            kernels.append((tuple(int(small >> i & 1) for i in range(active)), entries))
        result.append((mask, tuple(kernels)))
    return tuple(result)


def matrix_polynomial_power(local, variables, exponent):
    result = np.full((1,)*variables+(2, 2), -np.inf)
    result[(0,)*variables+(0, 0)] = 0.
    result[(0,)*variables+(1, 1)] = 0.
    for degree in range(exponent):
        following = np.full((degree+2,)*variables+(2, 2), -np.inf)
        for powers, matrix in local:
            target = tuple(slice(power, power+degree+1) for power in powers)
            for row in range(2):
                for column in range(2):
                    product = np.logaddexp2(result[..., row, 0]+matrix[0][column],
                                            result[..., row, 1]+matrix[1][column])
                    old = following[target+(row, column)]
                    following[target+(row, column)] = np.logaddexp2(old, product)
        result = following
    return result


def one_group_exact_bound(active, z):
    """First moment with exact MDS row-support counts; no fractional power."""
    words = one_group_region_words(z)
    row_counts = np.array([log2(mds_prescribed_support(Q, 32, 16, weight))
                           for weight in range(17, 33)])
    cases = []
    for mask, kernels in exact_mask_kernels(active):
        local = []
        for powers, support in kernels:
            matrix = tuple(tuple(total(logmass-target.bit_count()*log2(Q/(Q-1))
                                        +words[target][row][column]
                                        for target, logmass in support)
                                  for column in range(2)) for row in range(2))
            local.append((powers, matrix))
        polynomial = matrix_polynomial_power(local, active, 32)
        selected = polynomial[(slice(17, 33),)*active+(0, slice(None))]
        marked = np.logaddexp2(selected[..., 0], selected[..., 1])
        for axis in range(active):
            shape = [1]*active
            shape[axis] = 16
            marked = marked+row_counts.reshape(shape)
        cases.append(total(marked.ravel()))
    return 13.-log2(Q-1)-(2**21*7//20)*log2(z)+total(cases)


def one_group_exact_proposal(active, iterations=40):
    start = perf_counter()
    chosen = minimize_scalar(lambda ell: one_group_exact_bound(active, 2**(-ell)),
                             bounds=(.00001, .01), method='bounded',
                             options={'maxiter': iterations, 'xatol': 1e-9})
    return dict(active_groups=1, active_rows_per_group=active, mark=2**(-float(chosen.x)),
                log2_upper_screen=float(chosen.fun), exact_mds_support=True,
                exact_column_mask_mixture=True, alpha=1., directed=False,
                whole_code_certificate=False, elapsed_seconds=perf_counter()-start)


if __name__ == '__main__':
    import argparse
    import json
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--rows', type=int, nargs='+', default=(1, 128, 2048))
    parser.add_argument('--active', type=int, nargs='+', default=(1, 2, 7))
    parser.add_argument('--iterations', type=int, default=60)
    parser.add_argument('--regional', action='store_true')
    parser.add_argument('--fixed-probability', action='store_true')
    parser.add_argument('--one-group-correlated', action='store_true')
    parser.add_argument('--one-group-exact', action='store_true')
    parser.add_argument('--input-mark', action='store_true')
    args = parser.parse_args()
    for active in args.active:
        if args.one_group_exact:
            print(json.dumps(one_group_exact_proposal(active, args.iterations)), flush=True)
            continue
        if args.one_group_correlated:
            print(json.dumps(one_group_proposal(active, args.iterations)), flush=True)
            continue
        for rows in args.rows:
            print(json.dumps(proposal(rows, active, args.iterations, args.regional,
                                      not args.fixed_probability, args.input_mark)), flush=True)
