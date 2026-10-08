"""Exact finite window gate for direct-2K blockwise SPIN, without a pre-stage.

A positive result is a LOWER bound on distance failure over independent
uniform regional permutations. A zero result is not a distance certificate.
The terminal suffix refinement removes one state-exit constraint.
"""

from fractions import Fraction as F
from functools import lru_cache
from math import ceil, log2
import json

from high_distance_window_obstruction import mean_row_dimension

OUTERS = (64, 128, 256, 512)
WIDTHS = (16, 32, 64, 128)
STATES = (4, 8, 16, 32, 64)
DISTANCES = (F(2, 5), F(9, 20))


def failure_lower(mean, maximum, constraints):
    """If integer 0<=S<=maximum, return a lower bound on Pr[S>constraints]."""
    if constraints < 0 or maximum < 0 or not 0 <= mean <= maximum:
        raise ValueError('Invalid dimension mean, maximum, or constraint count')
    if maximum <= constraints or mean <= constraints:
        return F(0)
    return (F(mean)-constraints)/(maximum-constraints)


@lru_cache(maxsize=None)
def window_profile(length, width, distance, k_log=20):
    """No padding: the requested K must be divisible by k and R by width."""
    distance = F(distance)
    if length < 2 or length % 2 or width < 1 or k_log < 1 or not 0 < distance < F(1, 2):
        raise ValueError('Need even outer length and target strictly between zero and one half')
    dimension, requested = length//2, 1 << k_log
    if requested % dimension:
        raise ValueError('The direct-2K menu does not silently pad messages')
    rows = requested//dimension
    if rows % width:
        raise ValueError('Regions must contain complete packets')
    strict_count = (distance.numerator*rows-1)//distance.denominator
    window = (strict_count//width)*width
    mean = rows*mean_row_dimension(length, dimension, F(window, rows))
    return dict(length=length, dimension=dimension, rows=rows, width=width,
                distance=distance, window=window, mean=mean,
                maximum=dimension*window, n=length*rows, k=requested)


def poisson_binomial_coefficients(numerators, denominator):
    """Exact integer coefficients of product((D-u)+u*x), sharing denominator D."""
    coefficients = [1]
    for selected in numerators:
        if not 0 <= selected <= denominator:
            raise ValueError('Invalid selection probability')
        nxt = [0]*(len(coefficients)+1)
        for t, value in enumerate(coefficients):
            nxt[t] += (denominator-selected)*value
            nxt[t+1] += selected*value
        coefficients = nxt
    return coefficients


@lru_cache(maxsize=None)
def balanced_profile(length, width, distance, k_log=20):
    """Use the full strict weight budget with windows differing by one packet."""
    profile = window_profile(length, width, F(distance), k_log).copy()
    distance = F(distance)
    strict_weight = (distance.numerator*profile['n']-1)//distance.denominator
    packets = strict_weight//width
    low, high_regions = divmod(packets, length)
    denominator = profile['rows']//width
    numerators = (low,)*(length-high_regions)+(low+1,)*high_regions
    coefficients = poisson_binomial_coefficients(numerators, denominator)
    rho = length//2
    mean = profile['rows']*F(sum((t-rho)*coefficients[t] for t in range(rho+1, length+1)),
                            denominator**length)
    total = packets*width
    profile.update(window=low*width, high_window=(low+1)*width,
                   high_regions=high_regions, total_window=total, mean=mean,
                   maximum=F(total, 2), regional_packet_counts=numerators)
    return profile


def minimum_unexcluded_state(profile, setup_bits=40, terminal_suffix=True):
    """Smallest integer s for which this particular failure lower is <=2^-bits."""
    if setup_bits < 1:
        raise ValueError('Need a positive setup margin')
    epsilon = F(1, 1 << setup_bits)
    exits = profile['length']-int(terminal_suffix)
    threshold = ((profile['mean']-epsilon*profile['maximum'])
                 /(exits*(1-epsilon)))
    return max(0, ceil(threshold))


def screen(length=128, width=32, state=4, distance=F(2, 5),
           setup_bits=40, terminal_suffix=True, k_log=20, balanced=True):
    if state < 0 or setup_bits < 1:
        raise ValueError('Invalid state or setup margin')
    profile = (balanced_profile if balanced else window_profile)(length, width, F(distance), k_log)
    exits = length-int(terminal_suffix)
    bound = failure_lower(profile['mean'], profile['maximum'], exits*state)
    epsilon = F(1, 1 << setup_bits)
    return dict(outer=[length, length//2], rows=profile['rows'],
                message_k=profile['k'], output_n=profile['n'],
                packet_width=width, state=state, target_distance=str(F(distance)),
                window=profile['window'], window_fraction=str(F(profile['window'], profile['rows'])),
                high_window=profile.get('high_window', profile['window']),
                high_regions=profile.get('high_regions', 0), balanced=balanced,
                output_weight_upper=profile.get('total_window', length*profile['window']),
                state_exit_constraints=exits*state, terminal_suffix=terminal_suffix,
                expected_shortened_dimension=float(profile['mean']),
                failure_probability_lower=float(bound),
                setup_margin_upper_bits=(-log2(bound.numerator)+log2(bound.denominator)) if bound else None,
                excludes_requested_distance_and_margin=bound > epsilon,
                minimum_state_not_excluded=minimum_unexcluded_state(profile, setup_bits, terminal_suffix),
                exact_rational_comparison=True, zero_lower_bound_is_not_a_pass=True,
                probability_space='independent uniform regional permutations',
                no_accumulator_prestage=True)


def menu_summary():
    rows = []
    for length in OUTERS:
        for width in WIDTHS:
            row = dict(outer=[length, length//2], width=width)
            for distance in DISTANCES:
                label = str(distance)
                equal = window_profile(length, width, distance)
                profile = balanced_profile(length, width, distance)
                threshold = minimum_unexcluded_state(profile)
                row[label] = dict(window=profile['window'],
                                  high_window=profile['high_window'], high_regions=profile['high_regions'],
                                  mean=float(profile['mean']),
                                  old_bs_threshold=minimum_unexcluded_state(equal, terminal_suffix=False),
                                  equal_terminal_threshold=minimum_unexcluded_state(equal),
                                  terminal_threshold=threshold,
                                  excluded_states=[s for s in STATES if s < threshold],
                                  states_not_excluded=[s for s in STATES if s >= threshold])
            rows.append(row)
    return rows


def exhaustive_dimension_law(length, dimension, rows, window):
    """Exact small routing oracle, returning a distribution of shortening dimensions."""
    from collections import Counter
    from itertools import combinations, product
    choices = tuple(combinations(range(rows), window))
    total = len(choices)**length
    if total > 1_000_000:
        raise ValueError('This oracle is bounded to small proof tests')
    counts = Counter()
    for regions in product(choices, repeat=length):
        occupancy = [0]*rows
        for region in regions:
            for r in region:
                occupancy[r] += 1
        counts[sum(max(0, t-length+dimension) for t in occupancy)] += 1
    return {size: F(count, total) for size, count in counts.items()}


def main():
    for row in menu_summary():
        print(json.dumps(row))
    for length, width, state, distance in (
        (64, 16, 32, F(2, 5)), (64, 128, 64, F(9, 20)),
        (128, 16, 32, F(9, 20)), (256, 16, 4, F(9, 20))):
        print(json.dumps(screen(length, width, state, distance)))


if __name__ == '__main__':
    main()
