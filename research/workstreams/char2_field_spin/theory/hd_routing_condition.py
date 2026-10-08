"""A capped routing event for component subspaces, not a full certificate.

Unlike a codeword first moment, this counts each bad setup only once.
"""
from fractions import Fraction as F
from math import comb, prod


def falling(n, k):
    if k > n:
        return 0
    return prod(range(n-k+1, n+1))


def component_intervals(cuts_by_region, rows):
    """Supports of the globally carried free prefix-component coordinates."""
    cuts = sorted(c*rows+x for c, local in enumerate(cuts_by_region) for x in local)
    end = rows*len(cuts_by_region)
    return [tuple(range(left, right)) for left, right in zip(cuts, cuts[1:]+[end])]


def contained_probability(supports, selected, rows, window_sizes):
    """Joint containment under independent uniform regional permutations."""
    counts = [0]*len(window_sizes)
    for j in selected:
        for x in supports[j]:
            counts[x//rows] += 1
    return prod(F(falling(length, h), falling(rows, h))
                for length, h in zip(window_sizes, counts))


def mgf_upper(components, hit_fraction, tilt):
    hit_fraction, tilt = F(hit_fraction), F(tilt)
    if not 0 <= hit_fraction <= 1 or tilt < 1 or components < 0:
        raise ValueError('invalid positive-moment parameters')
    return (1+hit_fraction*(tilt-1))**components


def receipt(first_class=17820, precision=160):
    """All J with |J|>=first_class and all one-window-per-region choices."""
    from flint import arb, ctx
    previous = ctx.prec
    ctx.prec = precision
    try:
        rows, regions, rho = 32768, 64, 32
        # Threshold rho*a+128 is conservative: a positive nullspace from
        # dimension counting needs one more variable than these constraints.
        c = arb(24)/25
        ratio = arb(rows-first_class)/(first_class+1)*c**rho
        if not ratio < 1:
            raise ValueError('class tail is not decreasing geometrically')
        first = (arb(rows)**regions*comb(rows, first_class)*c**(rho*first_class)
                 *(arb(2)/3)**128)
        upper = first/(1-ratio)
        return dict(first_class=first_class, precision=precision,
                    log2_bad_route_upper=str(upper.log()/arb(2).log()),
                    below_2_to_minus_100=bool(upper < arb(2)**-100),
                    max_window_length=13107,
                    includes_all_row_sets_and_window_starts=True,
                    includes_components_crossing_region_boundaries=True,
                    whole_distance_certificate=False)
    finally:
        ctx.prec = previous


def interval_receipt(first_class=18300, max_runs=64, precision=160):
    """Cap a structural obstruction for all low-weight sets with few runs.

    Windows can have unequal lengths or cross regional boundaries. For each
    J, count only the a component-start positions in each region to bound
    complete containment. Their hypergeometric MGFs multiply; AM-GM then
    uses the total output weight, rather than a maximum regional weight.
    """
    from flint import arb, ctx
    previous = ctx.prec
    ctx.prec = precision
    try:
        rows, regions, rho = 32768, 64, 32
        n = rows*regions
        c = arb(24)/25
        ratio = arb(rows-first_class)/(first_class+1)*c**rho
        if not ratio < 1:
            raise ValueError('class tail is not decreasing geometrically')
        class_sum = comb(rows, first_class)*c**(rho*first_class)/(1-ratio)
        # Choosing 2g endpoints counts every g-run subset of [N]. Dropping
        # its weight restriction only increases the number of choices.
        windows = sum(comb(n+1, 2*g)*(arb(2)/3)**(2*g+1)
                      for g in range(1, max_runs+1))
        upper = class_sum*windows
        return dict(first_class=first_class, max_runs=max_runs, precision=precision,
                    log2_bad_route_upper=str(upper.log()/arb(2).log()),
                    below_2_to_minus_100=bool(upper < arb(2)**-100),
                    threshold='32*a + 2*g + 1 entirely contained components',
                    total_output_weight_at_most=838860,
                    no_individual_region_weight_restriction=True,
                    whole_distance_certificate=False)
    finally:
        ctx.prec = previous


if __name__ == '__main__':
    print(receipt())
    print(interval_receipt())
