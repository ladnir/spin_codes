"""Direct blockwise SPIN: exact structures, a failed moment, and capped traps.

The capped route event covers full active source rows in one packet-aligned
window per region, including connected windows. It is not whole-code distance.
"""

from fractions import Fraction as F
from math import comb

from hd_blockwise_geometry import balanced_profile


def _choose(n, k):
    return comb(n, k) if 0 <= k <= n else 0


def moment_obstruction(active=40, state=16, precision=160):
    """Exact routing event, actual outer row space, free terminal state."""
    from flint import arb, ctx
    previous = ctx.prec
    ctx.prec = precision
    try:
        b, k, q = 128, 64, 65536
        profile = balanced_profile(b, 64, F(2, 5))
        rows = profile['rows']
        lengths = [64*v for v in profile['regional_packet_counts']]
        dimension = k*active-(b-1)*state
        if not 1 <= active <= min(lengths) or dimension < 2:
            raise ValueError('Need a positive full-row trap dimension')
        probability = F(1)
        for length in lengths:
            probability *= F(comb(length, active), comb(rows, active))
        lp = (arb(probability.numerator).log()-arb(probability.denominator).log())/arb(2).log()
        lower = lp+16*(dimension-2)
        return dict(active_rows=active, state=state, message_dimension_lower=dimension,
                    low_weight=profile['total_window'], route_event_log2=str(lp),
                    projective_count=f'({q}^{dimension}-1)/{q-1}',
                    actual_projective_moment_log2_lower=str(lower),
                    moment_exceeds_one=bool(lower > 0), distance_failure_not_established=True)
    finally:
        ctx.prec = previous


def connected_terms(regions=128, rows=16384, width=64, state=16, distance=F(2, 5)):
    """Count full-row containment families, with actual connected exit counts.

All regions are used because a full outer row has one coordinate per region.
The formulas are the h=b specialization of connected_window_capacity.py.
"""
    k = regions//2
    blocks = rows//width
    limit = (distance.numerator*(regions*rows)-1)//distance.denominator
    budget = limit//width
    for components in range(1, regions+1):
        for singletons in range(components+1):
            if singletons == components:
                if regions != components:
                    continue
                compositions = 1
            else:
                compositions = _choose(regions-components-1, components-singletons-1)
                if not compositions:
                    continue
            for terminal in ('charged', 'singleton_tail', 'long_tail'):
                if terminal == 'charged':
                    placements = _choose(components, singletons)
                    free_exit, extra_full, anchored = 0, 0, 0
                elif terminal == 'singleton_tail':
                    placements = _choose(components-1, singletons-1)
                    free_exit, extra_full, anchored = 1, 0, 1
                else:
                    placements = _choose(components-1, singletons)
                    free_exit, extra_full, anchored = 1, 1, 0
                full = regions-2*components+singletons+extra_full
                variable = regions-full
                remaining = budget-full*blocks
                if not placements or full < 0 or variable < 1 or remaining < variable:
                    continue
                minimum = (state*(components-free_exit))//k+1
                families = (placements*compositions*_choose(remaining, variable)
                            * blocks**(singletons-anchored))
                p = F(remaining, blocks*variable)
                if not 0 < p < 1:
                    raise ArithmeticError('Invalid variable-window occupancy')
                yield dict(components=components, singletons=singletons, terminal=terminal,
                           full_regions=full, variable_regions=variable,
                           minimum_active_rows=minimum, window_families=families,
                           containment_per_row=p**variable, average_variable_fraction=p)


def capped_full_row_traps(state=16, precision=160):
    """Directed upper for ALL dimension-forced full-row window traps."""
    from flint import arb, ctx
    previous = ctx.prec
    ctx.prec = precision
    try:
        rows = 16384
        total = arb(0)
        count = 0
        worst_float = (-float('inf'), None)
        for term in connected_terms(state=state):
            a = term['minimum_active_rows']
            r = term['containment_per_row']
            ratio = F(rows-a, a+1)*r
            if not ratio < 1:
                raise ArithmeticError('The source-set geometric tail did not contract')
            ar = arb(r.numerator)/arb(r.denominator)
            aratio = arb(ratio.numerator)/arb(ratio.denominator)
            value = arb(term['window_families'])*comb(rows, a)*ar**a/(1-aratio)
            total += value
            logvalue = float(value.log()/arb(2).log())
            if logvalue > worst_float[0]:
                worst_float = (logvalue, {k: v for k, v in term.items()
                                         if k not in ('window_families', 'containment_per_row', 'average_variable_fraction')})
            count += 1
        upper = total.log()/arb(2).log()
        return dict(state=state, grouped_terms=count, capped_family_log2_upper=str(upper),
                    below_2_to_minus_40=bool(upper < -40), worst_group=worst_float[1],
                    scope='full source rows, one aligned interval per region, all connected exits',
                    partial_row_shortening_not_covered=True, whole_code_certificate=False)
    finally:
        ctx.prec = previous


def fixed_window_capacity(mark=F(721, 500), state=16, precision=160):
    """Directed NA Chernoff bound for one profile; union diagnostic is NOT a pass."""
    from flint import arb, ctx
    from window_dimension_mgf import moment
    previous = ctx.prec
    ctx.prec = precision
    try:
        b, rows, blocks = 128, 16384, 256
        mean = moment(b, b//2, F(2, 5), mark)
        mu = arb(mean.numerator)/arb(mean.denominator)
        z = arb(mark.numerator)/arb(mark.denominator)
        logarithm = (rows*mu.log()-((b-1)*state+1)*z.log())/arb(2).log()
        count = (blocks*(blocks+1)//2)**b
        union = logarithm+arb(count).log()/arb(2).log()
        return dict(state=state, rational_mark=str(mark), fixed_profile_log2_upper=str(logarithm),
                    naive_all_window_log2_upper=str(union), naive_union_closes=bool(union < -40),
                    scope='shortened-row dimension, charged127 exits; not all connected windows',
                    whole_code_certificate=False)
    finally:
        ctx.prec = previous


def main():
    import json
    for active in (39, 40, 64):
        print(json.dumps(moment_obstruction(active)))
    for state in (16, 32):
        print(json.dumps(capped_full_row_traps(state)))
    print(json.dumps(fixed_window_capacity()))


if __name__ == '__main__':
    main()
