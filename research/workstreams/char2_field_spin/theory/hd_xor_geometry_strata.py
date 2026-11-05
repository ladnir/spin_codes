"""Fixed-geometry multiplicity strata for the XOR2 two-constituent family.

The r labelled constituent pairs are chosen once, not again for each of the
32 RS columns.  For each maximum multiplicity mu, a universal column cap is
valid for every fixed geometry in that stratum.  Its positive support
polynomial is therefore raised to power 32 before multiplying by the exact
number of geometries.  This module never averages independently over a new
geometry in each column.

Regional matrices are inputs.  The optional scalar-route screen is only a
diagnostic for independent uniform scalar regional permutations; it is not
a certificate for the alternating paired route.
"""
from fractions import Fraction
from functools import lru_cache
from itertools import combinations
from math import comb, log2
from time import perf_counter

import numpy as np

from hd_xor_joint_columns import (
    LIMIT, Q, REGION_LENGTH, arb_upper_max, categories, joint_count_arb,
    joint_log_count, resource_envelope_arb, resource_envelope_float,
)

MAX_GROUPS = 32


@lru_cache(None)
def _assignment_count(groups, choices, maximum):
    """Labelled balls in labelled bins, with each bin bounded by maximum."""
    previous = [1]+[0]*groups
    for _ in range(choices):
        previous = [sum(comb(n, j)*previous[n-j]
                        for j in range(min(maximum, n)+1))
                    for n in range(groups+1)]
    return previous[groups]


def geometry_stratum_counts(groups, choices=28):
    """Exact counts by maximum multiplicity among choices labelled pairs.

    Equivalently, the cumulative count is
    r! [x^r] (sum_{j=0}^mu x^j/j!)^choices.  The recurrence uses integers only.
    """
    if not isinstance(groups, int) or not 1 <= groups <= MAX_GROUPS or not isinstance(choices, int) or choices < 1:
        raise ValueError('One to thirty-two groups and a positive choice count required')
    counts, previous = {}, 0
    for maximum in range(1, groups+1):
        cumulative = _assignment_count(groups, choices, maximum)
        counts[maximum] = cumulative-previous
        previous = cumulative
    if sum(counts.values()) != choices**groups:
        raise ArithmeticError('Geometry strata do not partition all labelled assignments')
    return counts


def _float_tables(regions, length):
    regions = np.asarray(regions, dtype=float)
    if regions.ndim == 3 and regions.shape[-2:] == (2, 2):
        tables = (regions,)*length
    elif regions.ndim == 4 and regions.shape[0] == length and regions.shape[-2:] == (2, 2):
        tables = tuple(regions)
    else:
        raise ValueError('Expected occupancy-indexed log matrices, optionally for each labelled region')
    if np.any(np.isnan(regions)) or np.any(np.isposinf(regions)):
        raise ValueError('Log matrices must be finite or negative infinity')
    return tables


def _arb_tables(regions, length, labelled_regions):
    tables = tuple(tuple(table) for table in regions) if labelled_regions else (tuple(regions),)*length
    if len(tables) != length or any(any(len(matrix) != 4 for matrix in table) for table in tables):
        raise ValueError('Expected occupancy-indexed flat-four Arb matrices')
    if any(not value.is_finite() or not value.upper() >= 0 for table in tables for matrix in table for value in matrix):
        raise ValueError('Regional matrices must have finite nonnegative entries')
    # Only upper bounds are used. Exact endpoints preserve nonnegative
    # monotonicity even when a coefficient enclosure has a small radius.
    return tuple(tuple(tuple(value.upper() for value in matrix) for matrix in table) for table in tables)


def _validate_interval(tables, lower, upper, pair_limit, length):
    if (not isinstance(length, int) or length < 1 or not 0 <= lower <= upper
            or any(upper >= len(table) for table in tables)
            or not isinstance(pair_limit, int)):
        raise ValueError('Invalid occupancy interval, pair limit, or word length')


def pair_resource_float(regions, log_potential, lower, upper, pair_limit, length=8):
    """Ordered upper vectors indexed by sum(m_i-lower), in log2 form.

    Admissible words satisfy lower<=m_i<=upper and m_i+m_j<=pair_limit
    for all distinct positions.  A state records shifted sum and largest
    previous occupancy.  Prepending a matrix preserves labelled chronology.
    Either one common table or ``(length, occupancy, 2, 2)`` tables are allowed.
    """
    tables = _float_tables(regions, length)
    _validate_interval(tables, lower, upper, pair_limit, length)
    if not np.isfinite(log_potential):
        raise ValueError('A finite log potential is required')
    width, limit = upper-lower, pair_limit-2*lower
    if np.asarray(regions).ndim == 3 and limit >= 2*width:
        return resource_envelope_float(regions, log_potential, lower, upper, length)
    # First populate a one-letter suffix: no pair restriction applies yet.
    values = np.full((width+1, width+1, 2), -np.inf)
    for delta in range(width+1):
        matrix = tables[-1][lower+delta]
        values[delta, delta] = np.logaddexp2(matrix[:, 0], matrix[:, 1]+log_potential)
    for table in reversed(tables[:-1]):
        following = np.full((width+1, values.shape[1]+width, 2), -np.inf)
        for delta in range(width+1):
            allowed = min(width, limit-delta)
            if allowed < 0:
                continue
            matrix = table[lower+delta]
            action = np.stack([np.logaddexp2(matrix[row, 0]+values[:allowed+1, :, 0],
                                            matrix[row, 1]+values[:allowed+1, :, 1])
                               for row in range(2)], axis=2)
            stop = min(delta, allowed)+1
            following[delta, delta:delta+values.shape[1]] = np.maximum(
                following[delta, delta:delta+values.shape[1]], np.max(action[:stop], axis=0))
            if allowed > delta:
                following[delta+1:allowed+1, delta:delta+values.shape[1]] = np.maximum(
                    following[delta+1:allowed+1, delta:delta+values.shape[1]], action[delta+1:])
        values = following
    return np.max(values, axis=0)


def pair_resource_arb(regions, potential, lower, upper, pair_limit, length=8, *, labelled_regions=False):
    """Directed counterpart; caller owns the precision context.

    Every max compares exact upper endpoints, not overlapping Arb balls.
    ``labelled_regions=True`` supplies one occupancy table per position.
    """
    from flint import arb
    tables = _arb_tables(regions, length, labelled_regions)
    _validate_interval(tables, lower, upper, pair_limit, length)
    if not potential > 0 or not potential.is_finite():
        raise ValueError('A positive finite potential is required')
    width, limit = upper-lower, pair_limit-2*lower
    if not labelled_regions and limit >= 2*width:
        return resource_envelope_arb(tables[0], potential, lower, upper, length)
    zero = arb(0)
    values = {}
    for delta in range(width+1):
        a, b, c, d = tables[-1][lower+delta]
        values[delta, delta] = ((a+b*potential).upper(), (c+d*potential).upper())
    for table in reversed(tables[:-1]):
        following = {}
        for (largest, degree), (first, second) in values.items():
            for delta in range(min(width, limit-largest)+1):
                a, b, c, d = table[lower+delta]
                key = max(largest, delta), degree+delta
                left, right = following.get(key, (zero, zero))
                following[key] = (arb_upper_max(left, a*first+b*second),
                                  arb_upper_max(right, c*first+d*second))
        values = following
    result = [(zero, zero) for _ in range(length*width+1)]
    for (_, degree), (first, second) in values.items():
        left, right = result[degree]
        result[degree] = arb_upper_max(left, first), arb_upper_max(right, second)
    return result


def category_resources(single, double, cancelled, mu):
    """Occupancy box, pairwise upper sum, and shifted total for one branch.

    If q_i singletons miss vertex i and e_i cancelled edges meet i, then
    m_i=s+d-k-q_i+e_i. Thus sum m_i=7s+8d-6k. Each cancelled edge type has
    multiplicity at most mu, so e_i<=min(k,7mu) and
    e_i+e_j<=k+min(k,mu). Dropping q_i>=0 gives the upper bounds below;
    q_i<=s gives m_i>=d-k. These bounds apply to one fixed geometry.
    """
    lower = double-cancelled
    # At most seven edge labels meet a vertex, each with multiplicity <=mu.
    upper = lower+single+min(cancelled, 7*mu)
    pair_limit = 2*single+2*double-cancelled+min(cancelled, mu)
    target = 7*single+2*cancelled
    return lower, upper, min(pair_limit, 2*upper), target


class MultiplicityCapsFloat:
    """Cache resource tables across mu, for one immutable regional operator bank."""
    def __init__(self, regions, log_potential, groups):
        if not 1 <= groups <= MAX_GROUPS:
            raise ValueError('One to thirty-two groups required')
        self.regions = np.array(regions, dtype=float, copy=True)
        self.regions.flags.writeable = False
        self.log_potential, self.groups = log_potential, groups
        tables = _float_tables(self.regions, 8)
        if any(len(table) <= groups for table in tables) or not np.isfinite(log_potential):
            raise ValueError('Incomplete regional tables or invalid potential')
        self._envelope = lru_cache(None)(lambda lower, upper, limit:
            pair_resource_float(self.regions, log_potential, lower, upper, limit))

    def __call__(self, mu):
        if not isinstance(mu, int) or not 1 <= mu <= self.groups:
            raise ValueError('Invalid maximum pair multiplicity')
        costs = {}
        for single, double in categories(self.groups):
            values = np.full(2, -np.inf)
            for cancelled in range(double+1):
                lower, upper, limit, target = category_resources(single, double, cancelled, mu)
                vectors = self._envelope(lower, upper, limit)
                if target >= len(vectors):
                    continue
                weight = (log2(comb(double, cancelled))
                          +(double-cancelled)*log2(Q-2)-double*log2(Q-1))
                values = np.logaddexp2(values, vectors[target]+weight)
            costs[single, double] = float(np.max(values-[0., self.log_potential]))
        return costs


class MultiplicityCapsArb:
    """Directed evaluator; construct and use inside one precision context."""
    def __init__(self, regions, potential, groups, *, labelled_regions=False):
        from flint import ctx
        if not 1 <= groups <= MAX_GROUPS:
            raise ValueError('One to thirty-two groups required')
        tables = _arb_tables(regions, 8, labelled_regions)
        if any(len(table) <= groups for table in tables) or not potential > 0 or not potential.is_finite():
            raise ValueError('Incomplete regional tables or invalid potential')
        self.regions = tables if labelled_regions else tables[0]
        self.potential, self.groups = potential, groups
        self.precision = ctx.prec
        self._envelope = lru_cache(None)(lambda lower, upper, limit:
            pair_resource_arb(self.regions, potential, lower, upper, limit,
                              labelled_regions=labelled_regions))

    def __call__(self, mu):
        from flint import arb, ctx
        if ctx.prec != self.precision:
            raise ValueError('Use a directed evaluator only in its construction precision context')
        if not isinstance(mu, int) or not 1 <= mu <= self.groups:
            raise ValueError('Invalid maximum pair multiplicity')
        caps = {}
        for single, double in categories(self.groups):
            first, second = arb(0), arb(0)
            for cancelled in range(double+1):
                lower, upper, limit, target = category_resources(single, double, cancelled, mu)
                vectors = self._envelope(lower, upper, limit)
                if target >= len(vectors):
                    continue
                weight = arb(comb(double, cancelled)*(Q-2)**(double-cancelled))/(Q-1)**double
                left, right = vectors[target]
                first += weight*left
                second += weight*right
            # Mix cancellation branches before taking the potential quotient.
            caps[single, double] = arb_upper_max(first, second/self.potential)
        return caps


def multiplicity_caps_float(regions, log_potential, groups, mu):
    return MultiplicityCapsFloat(regions, log_potential, groups)(mu)


def multiplicity_caps_arb(regions, potential, groups, mu, *, labelled_regions=False):
    return MultiplicityCapsArb(regions, potential, groups, labelled_regions=labelled_regions)(mu)


def _common_pair_cases(single):
    return tuple((left, right, at_left) for left, right in combinations(range(8), 2)
                 for at_left in range(single+1))


def common_pair_caps_float(regions, log_potential, groups):
    """Exact column-law caps for mu=r only, with every group on one edge.

    A category's singleton rows split between that edge's two endpoints.
    All 28 fixed edges are checked, not resampled across RS columns.  The
    supplied log-matrix tables may depend on the eight labelled positions.
    With j singletons on endpoint a and s-j on b, the occupancies are
    m_a=s-j+d, m_b=j+d, and m_i=s+d-k on the other six regions. The exact
    Binomial(d,1/255) law of k is mixed separately for every fixed edge/split.
    """
    if not 1 <= groups <= MAX_GROUPS or not np.isfinite(log_potential):
        raise ValueError('Valid group count and finite potential required')
    tables = _float_tables(regions, 8)
    if any(len(table) <= groups for table in tables):
        raise ValueError('Incomplete regional tables')
    costs = {}
    for single, double in categories(groups):
        cases = np.array(_common_pair_cases(single), dtype=int)
        indices = np.arange(len(cases))
        mixed = np.full((len(cases), 2), -np.inf)
        for cancelled in range(double+1):
            occupancies = np.full((len(cases), 8), single+double-cancelled, dtype=int)
            occupancies[indices, cases[:, 0]] = single-cases[:, 2]+double
            occupancies[indices, cases[:, 1]] = cases[:, 2]+double
            vectors = np.tile([0., log_potential], (len(cases), 1))
            for position in reversed(range(8)):
                matrices = tables[position][occupancies[:, position]]
                vectors = np.logaddexp2(matrices[:, :, 0]+vectors[:, 0, None],
                                       matrices[:, :, 1]+vectors[:, 1, None])
            weight = log2(comb(double, cancelled))+(double-cancelled)*log2(Q-2)-double*log2(Q-1)
            mixed = np.logaddexp2(mixed, vectors+weight)
        costs[single, double] = float(np.max(mixed-[0., log_potential]))
    return costs


def common_pair_caps_arb(regions, potential, groups, *, labelled_regions=False):
    """Directed common-pair cap; valid only for the maximum-multiplicity r stratum."""
    from flint import arb
    if not 1 <= groups <= MAX_GROUPS or not potential > 0 or not potential.is_finite():
        raise ValueError('Valid group count and positive finite potential required')
    tables = _arb_tables(regions, 8, labelled_regions)
    if any(len(table) <= groups for table in tables):
        raise ValueError('Incomplete regional tables')
    caps = {}
    for single, double in categories(groups):
        maximum = arb(0)
        for left, right, at_left in _common_pair_cases(single):
            mixed_first, mixed_second = arb(0), arb(0)
            for cancelled in range(double+1):
                occupancies = [single+double-cancelled]*8
                occupancies[left] = single-at_left+double
                occupancies[right] = at_left+double
                first, second = arb(1), potential
                for position in reversed(range(8)):
                    a, b, c, d = tables[position][occupancies[position]]
                    first, second = a*first+b*second, c*first+d*second
                weight = arb(comb(double, cancelled)*(Q-2)**(double-cancelled))/(Q-1)**double
                mixed_first += weight*first
                mixed_second += weight*second
            maximum = arb_upper_max(maximum, arb_upper_max(mixed_first, mixed_second/potential))
        caps[single, double] = maximum
    return caps


def strata_log_union(regions, groups, z, log_potential):
    """FLOAT full fixed-r homogeneous a=2 union for supplied regional bounds."""
    if not 0 < z < 1 or log_potential > 0:
        raise ValueError('This endpoint normalization requires 0<z<1 and potential<=1')
    start = perf_counter()
    evaluator = MultiplicityCapsFloat(regions, log_potential, groups)
    prefix = log2(comb(REGION_LENGTH, groups))-log2(Q-1)-LIMIT*log2(z)-log_potential
    rows = []
    for mu, count in geometry_stratum_counts(groups).items():
        if count:
            costs = evaluator(mu)
            bound = prefix+log2(count)+joint_log_count(costs, groups)
            rows.append(dict(maximum_pair_multiplicity=mu, geometry_count=count, failure_log2_screen=bound))
    return dict(active_groups=groups, stratum_rows=rows,
                failure_log2_screen=float(np.logaddexp2.reduce([row['failure_log2_screen'] for row in rows])),
                resource_tables=evaluator._envelope.cache_info().currsize,
                elapsed_seconds=perf_counter()-start, directed=False,
                geometry_fixed_across_columns=True, full_homogeneous_group_class=True,
                whole_code_certificate=False)


def scalar_route_screen(groups=17, z=.996548217, log_potential=-63.5):
    from hd_xor_multigroup import region_operators
    result = strata_log_union(region_operators(z, groups), groups, z, log_potential)
    result.update(route_scope='independent-uniform-scalar-regional-permutations',
                  paired_route_certificate=False, mark=z, log2_potential=log_potential)
    return result


def scalar_route_receipt(groups=17, z='996548217/1000000000', log_potential='-127/2', precision=256,
                         *, multiplicities=None):
    """Directed diagnostic, including failed bounds; never a paired-route receipt.

    RS shortened-support dimensions give the same aggregate-floor polynomial
    as hd_xor_joint_columns.  N_r(mu) counts the fixed constituent geometry
    once, outside that 32-column polynomial and the projective normalization.
    """
    from flint import arb, ctx
    from hd_xor_column_potential import arb_region_matrices
    z, log_potential = Fraction(z), Fraction(log_potential)
    if not 0 < z < 1 or log_potential > 0 or precision < 128:
        raise ValueError('Valid marks and at least 128 precision bits required')
    counts = geometry_stratum_counts(groups)
    selected = tuple(mu for mu, count in counts.items() if count) if multiplicities is None else tuple(multiplicities)
    if (not selected or len(set(selected)) != len(selected)
            or any(mu not in counts or counts[mu] == 0 for mu in selected)):
        raise ValueError('A nonempty set of possible multiplicity strata is required')
    start = perf_counter()
    with ctx.workprec(precision):
        zz = arb(z.numerator)/z.denominator
        potential = arb(2)**(arb(log_potential.numerator)/log_potential.denominator)
        evaluator = MultiplicityCapsArb(arb_region_matrices(z, groups), potential, groups)
        prefix = comb(REGION_LENGTH, groups)*zz**(-LIMIT)/(Q-1)/potential
        total, rows = arb(0), []
        def ceiling(probability):
            value = probability.log()/arb(2).log()
            integer = (value.upper()*10**12).ceil().unique_fmpz()
            if integer is None:
                raise ArithmeticError('Could not resolve a rational logarithm ceiling')
            return str(Fraction(int(integer), 10**12))
        for mu in selected:
            count = counts[mu]
            bound = prefix*count*joint_count_arb(evaluator(mu), groups)
            total += bound
            rows.append(dict(maximum_pair_multiplicity=mu, geometry_count=count,
                             failure_log2_upper=ceiling(bound)))
        return dict(active_groups=groups, active_rows_per_group=2, mark=str(z),
                    log2_potential=str(log_potential), failure_log2_upper=ceiling(total),
                    meets_40_bit_target=bool(total < arb(2)**-40), stratum_rows=rows,
                    field_order=Q, forbidden_weight=LIMIT, rs_length=32, rs_dimension=16,
                    total_support_floor=34*groups, polynomial_degree=64*groups,
                    source_pair_assignments=sum(counts[mu] for mu in selected), projective_normalizations=1,
                    route_scope='independent-uniform-scalar-regional-permutations',
                    paired_route_certificate=False, geometry_fixed_across_columns=True,
                    full_homogeneous_group_class=set(selected) == {mu for mu, count in counts.items() if count},
                    whole_code_certificate=False,
                    labelled_region_order=True, cancellation_count_mixture=True,
                    individual_support_floors_relaxed=True, exact_category_multiplicities=True,
                    resource_tables=evaluator._envelope.cache_info().currsize,
                    precision_bits=precision, directed=True, elapsed_seconds=perf_counter()-start)


# Two literal witnesses chosen by a bounded FLOAT diagnostic, not by any
# directed optimizer. Each fixed geometry belongs to exactly one range.
SAVED_SCALAR_R17_SCHEDULE = (
    (1, 15, '996248217/1000000000', '-127/2'),
    (16, 17, '996548217/1000000000', '-127/2'),
)


def replay_scalar_r17(precision=256):
    from flint import arb, ctx
    with ctx.workprec(precision):
        pieces = [scalar_route_receipt(17, z, potential, precision,
                                       multiplicities=range(first, last+1))
                  for first, last, z, potential in SAVED_SCALAR_R17_SCHEDULE]
        rows = [row for piece in pieces for row in piece['stratum_rows']]
        if sorted(row['maximum_pair_multiplicity'] for row in rows) != list(range(1, 18)):
            raise ArithmeticError('The scalar schedule does not partition all geometry strata')
        total = arb(0)
        for piece in pieces:
            ceiling = Fraction(piece['failure_log2_upper'])
            total += arb(2)**(arb(ceiling.numerator)/ceiling.denominator)
        value = total.log()/arb(2).log()
        integer = (value.upper()*10**12).ceil().unique_fmpz()
        if integer is None:
            raise ArithmeticError('Could not resolve the union ceiling')
        return dict(active_groups=17, active_rows_per_group=2,
                    failure_log2_upper=str(Fraction(int(integer), 10**12)),
                    meets_40_bit_target=bool(total < arb(2)**-40), pieces=pieces,
                    route_scope='independent-uniform-scalar-regional-permutations',
                    paired_route_certificate=False, geometry_fixed_across_columns=True,
                    full_homogeneous_group_class=True, whole_code_certificate=False,
                    precision_bits=precision, directed=True)


if __name__ == '__main__':
    import argparse
    import json
    parser = argparse.ArgumentParser(description='Bounded FLOAT scalar-route multiplicity diagnostic')
    parser.add_argument('--groups', type=int, default=17)
    parser.add_argument('--z', type=float, default=.996548217)
    parser.add_argument('--log-potential', type=float, default=-63.5)
    parser.add_argument('--precision', type=int, help='Replay the literal default witness with directed arithmetic')
    parser.add_argument('--schedule', action='store_true', help='Replay the frozen two-witness scalar r17 schedule')
    arguments = parser.parse_args()
    if arguments.schedule:
        if arguments.precision is None or arguments.groups != 17:
            parser.error('The saved r17 schedule requires --precision and --groups 17')
        result = replay_scalar_r17(arguments.precision)
    elif arguments.precision is None:
        result = scalar_route_screen(arguments.groups, arguments.z, arguments.log_potential)
    else:
        if arguments.z != .996548217 or arguments.log_potential != -63.5:
            parser.error('Directed CLI replay uses the literal default witness; call the API for other rationals')
        result = scalar_route_receipt(arguments.groups, precision=arguments.precision)
    print(json.dumps(result, indent=2))
