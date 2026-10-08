"""Directed full-row interval-span bound for outer256/t16/s4.

This does not cover partially shortened rows or multiple intervals per region.
"""

from fractions import Fraction as F


def order_statistic_counts(rows, active):
    """Exact subset counts by anchored length and minimum inclusive span."""
    from math import comb

    if not 1 <= active <= rows:
        raise ValueError('Need a nonempty feasible source set')
    anchored, span = [0]*(rows+1), [0]*(rows+1)
    for length in range(active, rows+1):
        anchored[length] = comb(length-1, active-1)
        if active > 1:
            span[length] = (rows-length+1)*comb(length-2, active-2)
    if active == 1:
        span[1] = rows
    return anchored, span


def new_geometry_moment(active=10, precision=160):
    """Actual-code first moment at outer256/t16/s4; not a failure lower bound."""
    from math import comb
    from flint import arb, ctx
    from hd_blockwise_geometry import balanced_profile

    previous = ctx.prec
    ctx.prec = precision
    try:
        b, k, rows, width, state = 256, 128, 8192, 16, 4
        profile = balanced_profile(b, width, F(2, 5))
        lengths = [width*x for x in profile['regional_packet_counts']]
        dimension = k*active-(b-1)*state
        if not 1 <= active <= min(lengths) or dimension < 2:
            raise ValueError('Need a nontrivial projective-count lower bound')
        probability = F(1)
        for length in lengths:
            probability *= F(comb(length, active), comb(rows, active))
        log_probability = (arb(probability.numerator).log()
                           -arb(probability.denominator).log())/arb(2).log()
        lower = log_probability+16*(dimension-2)
        return dict(outer=[b, k], width=width, state=state, active_rows=active,
                    low_weight=profile['total_window'], dimension_lower=dimension,
                    route_event_log2=str(log_probability),
                    actual_projective_moment_log2_lower=str(lower),
                    moment_exceeds_one=bool(lower > 0),
                    substantial_failure_probability_not_established=True)
    finally:
        ctx.prec = previous


def complete_row_window_cap(mark=F(1023, 1024), precision=160):
    """Directed order-statistic bound for the outer256/t16/s4 full-row family.

    Reuses the already tested component-pattern enumeration, but replaces its
    endpoint union with exact minimum-span/anchored-length distributions.
    All source sets can be trimmed to the minimum dimension-forcing size.
    Partial-row shortening and multiple intervals per region are NOT covered.
    """
    from math import comb
    from flint import arb, ctx
    from hd_blockwise_global import connected_terms

    if not 0 < mark < 1:
        raise ValueError('A lower-tail Chernoff mark must lie between zero and one')
    previous = ctx.prec
    ctx.prec = precision
    try:
        b, rows, width, state, budget = 256, 8192, 16, 4, 838848
        terms = tuple(connected_terms(b, rows, width, state))
        active_sizes = sorted({term['minimum_active_rows'] for term in terms})
        v = arb(mark.numerator)/arb(mark.denominator)
        weights = [arb(1)]
        for _ in range(rows):
            weights.append(weights[-1]*v)
        moment_powers = {}
        for active in active_sizes:
            denominator = comb(rows, active)
            for label, counts in zip(('anchor', 'span'), order_statistic_counts(rows, active)):
                moment = sum((weight*count for weight, count in zip(weights, counts)), arb(0))/denominator
                powers = [arb(1)]
                for _ in range(b):
                    powers.append(powers[-1]*moment)
                moment_powers[active, label] = powers
        fixed_weights = {full: v**(full*rows-budget)
                         for full in {term['full_regions'] for term in terms}}
        total, worst = arb(0), (-float('inf'), None)
        for term in terms:
            c, singletons = term['components'], term['singletons']
            terminal, full = term['terminal'], term['full_regions']
            active = term['minimum_active_rows']
            compositions = 1 if c == singletons else comb(b-c-1, c-singletons-1)
            if terminal == 'charged':
                placements = comb(c, singletons)
            elif terminal == 'singleton_tail':
                placements = comb(c-1, singletons-1)
            else:
                placements = comb(c-1, singletons)
            free_spans = singletons-int(terminal == 'singleton_tail')
            anchors = b-full-free_spans
            value = (placements*compositions*comb(rows, active)*fixed_weights[full]
                     * moment_powers[active, 'anchor'][anchors]
                     * moment_powers[active, 'span'][free_spans])
            total += value
            log_value = float(value.log()/arb(2).log())
            if log_value > worst[0]:
                worst = (log_value, dict(components=c, singletons=singletons,
                                        terminal=terminal, full_regions=full,
                                        anchors=anchors, free_spans=free_spans,
                                        minimum_active_rows=active))
        upper = total.log()/arb(2).log()
        return dict(outer=[256, 128], width=width, state=state, rows=rows,
                    low_weight_budget=budget, rational_mark=str(mark),
                    active_sizes=active_sizes, grouped_terms=len(terms),
                    capped_family_log2_upper=str(upper), below_2_to_minus_40=bool(upper < -40),
                    worst_group=worst[1],
                    scope='dimension-forced complete source rows, one aligned interval per region, all connected exits',
                    partial_row_shortening_not_covered=True, whole_code_certificate=False)
    finally:
        ctx.prec = previous


def main():
    import json
    print(json.dumps(new_geometry_moment()))
    print(json.dumps(complete_row_window_cap()))


if __name__ == '__main__':
    main()
