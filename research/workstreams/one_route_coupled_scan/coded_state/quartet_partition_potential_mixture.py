"""Bounded FLOAT two-term potential control with fixed selector markers.

For common(beta,d,z,w), let V_j(u,a)=A_j^u B_j^a and
V=V_1+2^offset V_2. Each cached component action is bounded by its
coefficient-Chernoff row comparison, tightened by the valid uniform
pointwise comparison. Linearity bounds K V by the sum of these actions.
The maximum of (K V)/V is taken over every reachable population row.
The alternative V=min(V_1,2^offset V_2) uses the positive-operator
inequality K V<=min(K V_1,2^offset K V_2). This is a pointwise minimum
of two physical-state functions, not an entrywise minimum of
overlapping-basis operators.

Xi <= V(0,0)*lambda^(4095)
      *max_(u,a) beta^u(1+w)^(512-u)/V(u,a).

The grid keeps code markers fixed and varies potential pairs symmetrically
about the supplied centre, along u,a,and a-u/4. Three separations and three
crossovers give27 pairs plus the identical-component control, in each
of the sum and minimum modes. Component
arrays are cached with a fixed memory cap. All results are FLOAT controls,
not directed certificates or claims of optimality.
"""

import argparse
from fractions import Fraction
from itertools import combinations
from math import log2
from operator import index
import os
import shutil
import threading
import time

import numpy as np

try:
    partition_averaged_rows
except NameError:
    from quartet_partition_averaged_potential import (
        partition_averaged_rows, _pa_choose, PA_PROFILE, PA_LOGS,
        PA_WIDTH, PA_UPDATES, PA_M, PA_W)


PM_DIRECTIONS = (('u', 1.0, 0.0), ('a', 0.0, 1.0), ('a_minus_u_over4', -0.25, 1.0))
PM_DELTAS = (1/16, 1/8, 1/4)
PM_CROSSOVERS = (128, 256, 384)


class PartitionPotentialMixture:
    """Cache at most64 component comparisons in one common physical row order."""

    def __init__(self, log_beta, log_d, log_z, log_w, width=512, max_iterations=40):
        self.common = tuple(map(float, (log_beta, log_d, log_z, log_w)))
        if not np.all(np.isfinite(self.common)):
            raise ValueError('Require finite common marker logs')
        self.width, self.iterations = index(width), index(max_iterations)
        self.rows = None
        self.cache = {}

    def component(self, log_A, log_B):
        key = (float(log_A), float(log_B))
        if key in self.cache:
            return self.cache[key]
        if len(self.cache) >= 64:
            raise ValueError('Bounded cache permits at most64 components')
        result = partition_averaged_rows(*self.common, *key, width=self.width,
                                          physical_rows=True, max_iterations=self.iterations)
        if self.rows is None:
            self.rows = result['rows'].copy()
        elif not np.array_equal(self.rows, result['rows']):
            raise ArithmeticError('Mixture components must use identical physical row ordering')
        log_V = self.rows[:, 0]*key[0]+self.rows[:, 1]*key[1]
        ratio = np.minimum(result['row_log2_upper'], result['pointwise_log2_growth_upper'])
        value = dict(log_V=log_V, log_action_upper=log_V+ratio,
                     raw_chernoff_growth=result['raw_chernoff_log2_growth_upper'],
                     pointwise_growth=result['pointwise_log2_growth_upper'],
                     optimizer=result['optimizer'])
        self.cache[key] = value
        return value

    def pair(self, first, second, offset=0.0, updates=4096, mode='sum'):
        """Full-row mixture growth, exact free-endpoint formula, and initial factor."""
        offset, updates = float(offset), index(updates)
        if not np.isfinite(offset) or abs(offset) > 4096 or updates < 1 or mode not in ('sum', 'minimum'):
            raise ValueError('Require finite bounded offset and at least one update')
        left, right = self.component(*first), self.component(*second)
        combine = np.logaddexp2 if mode == 'sum' else np.minimum
        log_V = combine(left['log_V'], offset+right['log_V'])
        log_action = combine(left['log_action_upper'], offset+right['log_action_upper'])
        ratios = log_action-log_V
        u = self.rows[:, 0]
        log_beta, _, _, log_w = self.common
        terminals = u*log_beta+(self.width-u)*np.logaddexp2(0.0, log_w)-log_V
        growth_index, endpoint_index = int(np.argmax(ratios)), int(np.argmax(terminals))
        growth, endpoint = float(ratios[growth_index]), float(terminals[endpoint_index])
        initial = float(combine(0.0, offset))
        moment = initial+(updates-1)*growth+endpoint
        return dict(mode=mode, first_potential=list(map(float, first)), second_potential=list(map(float, second)),
                    log2_second_coefficient=offset, physical_row_count=len(self.rows),
                    log2_growth_upper=growth, maximum_growth_row=self.rows[growth_index].tolist(),
                    log2_terminal_ratio=endpoint, maximum_terminal_row=self.rows[endpoint_index].tolist(),
                    log2_initial_potential=initial, log2_moment_upper=moment,
                    first_component_growth=min(left['raw_chernoff_growth'], left['pointwise_growth']),
                    second_component_growth=min(right['raw_chernoff_growth'], right['pointwise_growth']),
                    maximum_inner_mean_residual=max(left['optimizer']['maximum_mean_residual'],
                                                     right['optimizer']['maximum_mean_residual']),
                    common_physical_basis=True, positive_linearity=(mode == 'sum'),
                    physical_function_minimum_comparison=(mode == 'minimum'),
                    overlapping_basis_entrywise_minimum=False, free_endpoint=True,
                    floating_diagnostic_only=True, full_code_certificate=False)


def partition_potential_mixture_probe(profile=PA_PROFILE, death_upper=None, logs=PA_LOGS,
                                      max_iterations=40):
    """Fixed27-pair grid plus control in each of two modes; no marker search."""
    h, rank, volume, births, footprint = map(index, profile)
    logs = np.asarray(logs, dtype=np.float64)
    if logs.shape != (6,) or not np.all(np.isfinite(logs)) or np.any(logs[4:] < 0):
        raise ValueError('Require six finite logs(z,w,A,B,b,k), b/k logs nonnegative')
    N = PA_WIDTH*PA_UPDATES
    if not (1 <= h <= PA_M and 128*h <= rank <= 1 << 20 and 1 <= births <= volume <= footprint <= N):
        raise ValueError('Invalid actual-state profile')
    G, m, CH = footprint-births, volume-births, footprint-rank+128*h
    if G % 4 or G > N-PA_WIDTH:
        raise ValueError('Invalid live-selector volume')
    lower, upper = max(0, births-PA_WIDTH), min(footprint-volume, PA_W//2)
    D0 = min(PA_W-volume if death_upper is None else index(death_upper), upper)
    if D0 < lower or not (129*h <= CH <= 256*h and 0 <= CH-births <= G and 0 <= 256*h-CH <= N-footprint):
        return dict(profile=tuple(profile), status='empty', full_code_certificate=False)
    log_z, log_w, log_A, log_B, log_b, log_k = logs
    route = min(0.0, _pa_choose(G, CH-births)+_pa_choose(N-footprint, 256*h-CH)-255*_pa_choose(PA_M, h))
    correction = route-G*log_z-births*log_w-m*log_b-(m-D0)*log_k
    comparison = PartitionPotentialMixture(log_b+log_k, -log_k, log_z, log_w, max_iterations=max_iterations)
    center = (log_A, log_B)
    candidates = []
    for mode in ('sum', 'minimum'):
        control = comparison.pair(center, center, 0.0, mode=mode)
        control.update(kind='identical_component_control', raw_log2_probability_bound=float(control['log2_moment_upper']+correction))
        candidates.append(control)
    controls = list(candidates)
    for name, direction_A, direction_B in PM_DIRECTIONS:
        for delta in PM_DELTAS:
            first = (log_A+delta*direction_A, log_B+delta*direction_B)
            second = (log_A-delta*direction_A, log_B-delta*direction_B)
            for crossover in PM_CROSSOVERS:
                for mode in ('sum', 'minimum'):
                    result = comparison.pair(first, second, 2*delta*crossover, mode=mode)
                    result.update(kind='symmetric_pair', direction=name, delta=delta, crossover=crossover,
                                  raw_log2_probability_bound=float(result['log2_moment_upper']+correction))
                    candidates.append(result)
    if len(candidates) > 200:
        raise ArithmeticError('Fixed candidate budget exceeded')
    ordered = sorted(candidates, key=lambda candidate: candidate['raw_log2_probability_bound'])
    return dict(profile=tuple(map(int, profile)), death_interval=(lower, D0), logs=list(map(float, logs)),
                status='bounded_mixture_control', code_markers_fixed=True, candidate_count=len(candidates),
                component_count=len(comparison.cache), cached_array_bytes=sum(
                    entry['log_V'].nbytes+entry['log_action_upper'].nbytes for entry in comparison.cache.values()),
                route_log2=route, external_log2_correction=float(correction),
                controls=controls, best=ordered[0], five_best=ordered[:5],
                births_integrated=True, H_union_already_in_route=True,
                optimizer_optimality_proved=False, scalar_family_impossibility_proved=False,
                floating_diagnostic_only=True, full_code_certificate=False)


def _pm_subsets(mask):
    subset = mask
    while True:
        yield subset
        if subset == 0:
            break
        subset = (subset-1)&mask


def _pm_transitions(old, previous, partition, beta, d, z, w):
    for selector in range(1 << len(partition)):
        blocks = [block for j, block in enumerate(partition) if selector >> j & 1]
        new = sum(blocks)
        for births in _pm_subsets(new&~old):
            if any(not (block&old) and not (block&births) for block in blocks):
                continue
            allowed = (old&new)|births
            weight = (beta**(old&new).bit_count()*d**(previous&~allowed).bit_count()
                      *z**new.bit_count()*w**births.bit_count())
            yield new, allowed, weight


def partition_potential_mixture_checks():
    linearity = minimum_checks = endpoints = float_rows = normalization = moments = 0
    beta, d, z, w = Fraction(2), Fraction(1, 2), Fraction(1, 2), Fraction(1, 8)
    first, second, coefficient = (Fraction(2), Fraction(1, 2)), (Fraction(3, 2), Fraction(2, 3)), Fraction(3, 2)
    for width in (4, 8):
        comparison = PartitionPotentialMixture(*(log2(float(x)) for x in (beta, d, z, w)), width=width, max_iterations=64)
        numeric = {mode: comparison.pair(tuple(log2(float(x)) for x in first),
                                         tuple(log2(float(x)) for x in second), log2(float(coefficient)),
                                         updates=3, mode=mode) for mode in ('sum', 'minimum')}
        partitions = [(15,)] if width == 4 else [
            (block, 255^block) for rest in combinations(range(1, 8), 3)
            for block in (1+sum(1 << j for j in rest),)]
        exact_ratios = {mode: [] for mode in numeric}
        exact_endpoints = {mode: [] for mode in numeric}
        for u, a in comparison.rows:
            u, a = int(u), int(a)
            old, previous = (1 << u)-1, (1 << a)-1
            left = right = joint = joint_min = Fraction(0)
            for partition in partitions:
                for new, allowed, weight in _pm_transitions(old, previous, partition, beta, d, z, w):
                    V1 = first[0]**new.bit_count()*first[1]**allowed.bit_count()
                    V2 = second[0]**new.bit_count()*second[1]**allowed.bit_count()
                    left += weight*V1
                    right += weight*V2
                    joint += weight*(V1+coefficient*V2)
                    joint_min += weight*min(V1, coefficient*V2)
            assert joint == left+coefficient*right
            assert joint_min <= min(left, coefficient*right)
            linearity += 1
            minimum_checks += 1
            V1, V2 = first[0]**u*first[1]**a, second[0]**u*second[1]**a
            for mode in numeric:
                denominator = V1+coefficient*V2 if mode == 'sum' else min(V1, coefficient*V2)
                ratio = (joint if mode == 'sum' else joint_min)/(len(partitions)*denominator)
                endpoint = beta**u*(1+w)**(width-u)/denominator
                assert ratio <= max(left/(len(partitions)*V1), right/(len(partitions)*V2))
                exact_ratios[mode].append(ratio)
                exact_endpoints[mode].append(endpoint)
                assert log2(float(ratio)) <= numeric[mode]['log2_growth_upper']+2e-10
                assert log2(float(endpoint)) <= numeric[mode]['log2_terminal_ratio']+2e-10
                endpoints += 1
                float_rows += 1
        centre = tuple(log2(float(x)) for x in first)
        for mode in numeric:
            same0 = comparison.pair(centre, centre, 0.0, updates=3, mode=mode)
            same9 = comparison.pair(centre, centre, 9.0, updates=3, mode=mode)
            assert abs(same0['log2_moment_upper']-same9['log2_moment_upper']) < 1e-9
            normalization += 1
        if width == 4:
            states = {(0, 0): Fraction(1)}
            for updates in range(1, 5):
                exact = sum(weight*beta**old.bit_count()*(1+w)**(width-old.bit_count())
                            for (old, previous), weight in states.items())
                for mode in numeric:
                    lam, terminal = max(exact_ratios[mode]), max(exact_endpoints[mode])
                    initial = 1+coefficient if mode == 'sum' else min(Fraction(1), coefficient)
                    assert exact <= initial*lam**(updates-1)*terminal
                    moments += 1
                following = {}
                for (old, previous), weight in states.items():
                    for new, allowed, step in _pm_transitions(old, previous, (15,), beta, d, z, w):
                        key = (new, allowed)
                        following[key] = following.get(key, Fraction(0))+weight*step
                states = following
    print(dict(exact_mixture_linearity_rows=linearity, exact_minimum_comparison_rows=minimum_checks,
               exact_free_endpoint_rows=endpoints,
               float_mixture_row_controls=float_rows, identical_component_normalizations=normalization,
               exact_temporal_moments=moments, full_code_certificate=False), flush=True)


def partition_potential_mixture_main():
    parser = argparse.ArgumentParser(description=__doc__)
    modes = parser.add_mutually_exclusive_group(required=True)
    modes.add_argument('--checks', action='store_true')
    modes.add_argument('--probe', action='store_true')
    parser.add_argument('--profile', type=int, nargs=5, default=PA_PROFILE)
    parser.add_argument('--death-upper', type=int)
    parser.add_argument('--logs', type=float, nargs=6, default=PA_LOGS)
    parser.add_argument('--iterations', type=int, default=40)
    args = parser.parse_args()
    import resource
    if shutil.disk_usage('/tmp').free < 20*(1 << 30):
        raise RuntimeError('20GiB /tmp free-space guard')
    resource.setrlimit(resource.RLIMIT_AS, (2*(1 << 30), 2*(1 << 30)))
    timer = threading.Timer(175, lambda: os._exit(124))
    timer.daemon = True
    timer.start()
    started = time.monotonic()
    try:
        if args.checks:
            partition_potential_mixture_checks()
        else:
            print(partition_potential_mixture_probe(args.profile, args.death_upper, args.logs,
                                                   max_iterations=args.iterations), flush=True)
        print(dict(seconds=round(time.monotonic()-started, 3)), flush=True)
    finally:
        timer.cancel()


if __name__ == '__main__':
    partition_potential_mixture_main()
