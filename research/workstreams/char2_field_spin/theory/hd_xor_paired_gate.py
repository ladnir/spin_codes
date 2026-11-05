"""Bounded two-constituent XOR2 gate for alternating independently swapped pairs.

Two fixed message indices are counted once: F, the largest active-active
pair count in either cyclic matching, and mu, the largest constituent-pair
multiplicity. They are not resampled between RS columns. Each rectangle
uses caps valid for its upper F and mu, and its exact number of assignments.
This is a homogeneous class at one group count, never a whole-code proof.
"""
from fractions import Fraction
from hashlib import sha256
import json
from math import comb, log2
from pathlib import Path
from time import perf_counter

import numpy as np

from hd_xor_joint_columns import joint_log_count, joint_count_arb
from hd_xor_column_potential import Q, REGION_LENGTH, LIMIT


CONSTRUCTION = 'xor2-gf256-h4-alternating-swapped-pairs-k20-rate-half'


def default_plan(groups):
    """A small complete partition; it does not fit any numerical witness."""
    def bands(end, cuts, first):
        result = []
        for high in sorted({min(end, value) for value in cuts} | {end}):
            if high >= first:
                result.append((first, high))
                first = high+1
        return result
    return [(f, mu) for f in bands(groups//2, (0, 1, 2), 0)
            for mu in bands(groups, (2, 4, 8, 16), 1)]


def validate_plan(plan, groups):
    """Require disjoint rectangles; return whether they cover the full index set."""
    if type(groups) is not int or not 1 <= groups <= 17:
        raise ValueError('The paired-route gate currently supports one through seventeen groups')
    covered = set()
    for rectangle in plan:
        if len(rectangle) != 2 or any(len(interval) != 2 for interval in rectangle):
            raise ValueError('Each rectangle needs F and mu intervals')
        f, mu = rectangle
        if (any(type(x) is not int for interval in rectangle for x in interval)
                or not 0 <= f[0] <= f[1] <= groups//2
                or not 1 <= mu[0] <= mu[1] <= groups):
            raise ValueError('Invalid geometry rectangle')
        points = {(i, j) for i in range(f[0], f[1]+1) for j in range(mu[0], mu[1]+1)}
        if covered & points:
            raise ValueError('Geometry rectangles must be disjoint')
        covered |= points
    return covered == {(i, j) for i in range(groups//2+1) for j in range(1, groups+1)}


def restricted_regions(table, groups, maximum_full_pairs, *, directed=False):
    """For F<=maximum, every output subset has f<=min(F,floor(m/2)).

    The stronger lower bound depending on exact F is deliberately not used.
    Entrywise maxima may use different f; they remain upper bounds for every
    fixed support under both ordered matchings.
    """
    if not 0 <= maximum_full_pairs <= groups//2:
        raise ValueError('Invalid active-group pair cap')
    if directed:
        from hd_xor_joint_columns import arb_upper_max
        result = []
        for m in range(groups+1):
            values = list(table[m, 0])
            for f in range(1, min(maximum_full_pairs, m//2)+1):
                values = [arb_upper_max(left, right) for left, right in zip(values, table[m, f])]
            result.append(tuple(values))
        return result
    return np.array([np.maximum.reduce([table[m, f]
                     for f in range(min(maximum_full_pairs, m//2)+1)])
                     for m in range(groups+1)])


def count_rectangles(plan, groups):
    from hd_xor_geometry_strata import geometry_stratum_counts
    from hd_xor_group_geometry import matching_max_counts
    source_counts = matching_max_counts(REGION_LENGTH, groups)
    pair_counts = geometry_stratum_counts(groups)
    if sum(source_counts.values()) != comb(REGION_LENGTH, groups) or sum(pair_counts.values()) != 28**groups:
        raise ArithmeticError('Geometry counts do not cover all assignments')
    return [(sum(source_counts.get(i, 0) for i in range(f[0], f[1]+1)),
             sum(pair_counts.get(i, 0) for i in range(mu[0], mu[1]+1))) for f, mu in plan]


def gate(groups=17, z='996548217/1000000000', log_potential='-127/2', *,
         plan=None, precision=None, progress=None):
    from hd_xor_geometry_strata import MultiplicityCapsFloat, MultiplicityCapsArb
    from hd_xor_pair_route import pair_regions_screen, arb_pair_regions
    z, log_potential = Fraction(z), Fraction(log_potential)
    if not 0 < z < 1 or log_potential > 0:
        raise ValueError('Invalid exact weight mark or potential exponent')
    if precision is not None and (type(precision) is not int or precision < 128):
        raise ValueError('Directed replay requires at least128 bits')
    plan = default_plan(groups) if plan is None else plan
    complete = validate_plan(plan, groups)
    counts = count_rectangles(plan, groups)
    started = perf_counter()
    entries = []
    common = dict(construction=CONSTRUCTION, active_groups=groups,
        active_rows_per_group=2, source_dimension=2**20, output_length=2**21,
        forbidden_weight=LIMIT, mark=str(z), log2_potential=str(log_potential),
        pair_swaps=True, source_pair_geometry='alternating cyclic matchings',
        constituent_pair_geometry='fixed across all32 columns; multiplicity strata',
        complete_geometry_partition=complete, whole_code_certificate=False,
        support_scope='one homogeneous two-constituent active-group class',
        projective_normalizations=1, individual_support_floors_relaxed=True)
    if precision is None:
        table = pair_regions_screen(float(z), groups)
        prefix = -log2(Q-1)-LIMIT*log2(float(z))-float(log_potential)
        all_bounds = []
        evaluators = {}
        for rectangle, (supports, pairs) in zip(plan, counts):
            if not supports or not pairs:
                continue
            full_cap = rectangle[0][1]
            if full_cap not in evaluators:
                evaluators[full_cap] = MultiplicityCapsFloat(
                    restricted_regions(table, groups, full_cap), float(log_potential), groups)
            caps = evaluators[full_cap](rectangle[1][1])
            value = prefix+log2(supports)+log2(pairs)+joint_log_count(caps, groups)
            row = dict(ranges=rectangle, source_supports=supports, constituent_pairs=pairs,
                       failure_log2_screen=value, worst_mu=rectangle[1][1], worst_F=rectangle[0][1])
            entries.append(row)
            all_bounds.append(value)
            if progress:
                progress(row)
        bound = float(np.logaddexp2.reduce(all_bounds)) if all_bounds else None
        return dict(common, entries=entries, directed=False, failure_log2_screen=bound,
                    elapsed_seconds=perf_counter()-started)

    from flint import arb, ctx
    from hd_shifted_coverage import rational_ceiling
    with ctx.workprec(precision):
        table = arb_pair_regions(z, groups)
        potential = arb(2)**(arb(log_potential.numerator)/log_potential.denominator)
        zz = arb(z.numerator)/z.denominator
        prefix = zz**(-LIMIT)/(Q-1)/potential
        total = arb(0)
        evaluators = {}
        for rectangle, (supports, pairs) in zip(plan, counts):
            if not supports or not pairs:
                continue
            full_cap = rectangle[0][1]
            if full_cap not in evaluators:
                evaluators[full_cap] = MultiplicityCapsArb(
                    restricted_regions(table, groups, full_cap, directed=True), potential, groups)
            caps = evaluators[full_cap](rectangle[1][1])
            value = prefix*supports*pairs*joint_count_arb(caps, groups)
            if not value > 0 or not value.is_finite():
                raise ArithmeticError('Nonpositive or nonfinite geometry contribution')
            logarithm = value.log()/arb(2).log()
            row = dict(ranges=rectangle, source_supports=supports, constituent_pairs=pairs,
                       failure_log2_upper=str(rational_ceiling(logarithm)),
                       worst_mu=rectangle[1][1], worst_F=rectangle[0][1])
            entries.append(row)
            total += value
            if progress:
                progress(row)
        if not total > 0 or not total.is_finite():
            raise ArithmeticError('Empty or nonfinite geometry union')
        logarithm = total.log()/arb(2).log()
        return dict(common, entries=entries, directed=True, precision_bits=precision,
                    failure_log2_upper=str(rational_ceiling(logarithm)),
                    homogeneous_class_reaches_40_bits=bool(complete and logarithm < -40),
                    elapsed_seconds=perf_counter()-started)


if __name__ == '__main__':
    import argparse
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--groups', type=int, default=17)
    parser.add_argument('--mark', default='996548217/1000000000')
    parser.add_argument('--potential', default='-127/2')
    parser.add_argument('--plan', type=Path, help='Explicit JSON list of disjoint F/mu rectangles')
    parser.add_argument('--precision', type=int, help='Fresh directed replay, not FLOAT')
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    if args.output.exists():
        parser.error('Frozen gate results are never overwritten')
    plan = None if args.plan is None else json.loads(args.plan.read_bytes())
    result = gate(args.groups, args.mark, args.potential, plan=plan, precision=args.precision,
                  progress=lambda row: print(json.dumps(row, allow_nan=False), flush=True))
    result['source_sha256'] = sha256(Path(__file__).read_bytes()).hexdigest()
    if args.plan:
        result['plan_sha256'] = sha256(args.plan.read_bytes()).hexdigest()
    with args.output.open('x', encoding='utf-8') as target:
        json.dump(result, target, allow_nan=False, indent=2)
        target.write('\n')
    print(json.dumps(dict(path=str(args.output.resolve()), directed=result['directed'],
                          bound=result.get('failure_log2_upper', result.get('failure_log2_screen')),
                          elapsed_seconds=result['elapsed_seconds'])), flush=True)
