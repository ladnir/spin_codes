"""General finite regional-window capacity gates, not distance certificates.

For a fixed MDS [b,k] outer, h selected coordinate regions leave h-(b-k)
source dimensions per row. A box contains a row when that row's routed
coordinate lies in its prescribed window in every selected region.
Occupancy a with a*(h-(b-k)) > s*h gives a nonzero word confined to the
windows for every stated invertible causal s-state inner. A suffix box
fixes the last region and its final window, saving the last s constraints.

The quick scan uses a Poisson occupancy approximation ONLY to select
parameters. Exact checks use alternating rational occupancy tails and
normalized mixed binomial moments. Passing this necessary gate never
establishes a whole-code distance or setup-failure guarantee.
"""

import argparse
from fractions import Fraction
import json
import math

from scipy.stats import poisson

from rate_design_window_boxes import choose, display, falling, tail_partial


def geometry(k_log=20, length=33, dimension=24, state=4, regions=10,
             windows=3, suffix=False, period=64, width=16, distance=Fraction(1, 10), checksum=False):
    parity = length-dimension
    if not (0 < dimension < length and parity < regions <= length and state >= 1
            and period % width == 0 and windows >= 1 and 0 < distance < 1):
        raise ValueError('invalid regional capacity geometry')
    rows = (((1 << k_log)+dimension*period-1)//(dimension*period))*period
    n = rows*length
    threshold = (n*distance.numerator+distance.denominator-1)//distance.denominator-1
    window = min(rows, threshold//(regions*width)*width)
    if window < 1:
        raise ValueError('empty windows')
    free = regions-int(suffix)
    minimum = state*free//(regions-parity)+1
    if minimum > window:
        raise ValueError('too few positions for the kernel witness')
    if windows == 1:
        starts = (rows-window,)
    else:
        starts = tuple(((rows-window)*i//((windows-1)*width))*width for i in range(windows))
    if len(set(starts)) != windows:
        raise ValueError('window starts coincide')
    return dict(k_log=k_log, rows=rows, k=rows*dimension, n=n, length=length, dimension=dimension,
                rate=Fraction(dimension, length), parity=parity, state=state, regions=regions,
                windows=windows, suffix=suffix, checksum=checksum, period=period, width=width,
                free_regions=free, minimum=minimum, window=window, starts=starts,
                threshold=threshold, total_window=window*regions,
                kernel_dimension=minimum*(regions-parity)-state*free,
                box_count=choose(length-int(suffix), free)*windows**(free-int(checksum)))


def log_falling(n, degree):
    return sum(math.log(n-i) for i in range(degree))


def approximate(record):
    rows, window, regions, minimum = (record[key] for key in ('rows', 'window', 'regions', 'minimum'))
    mean = rows*(window/rows)**regions
    log_tail = float(poisson.logsf(minimum-1, mean))
    # scipy may underflow for the very small tails. The leading Poisson
    # term remains useful for selecting exact checks in that regime.
    if not math.isfinite(log_tail):
        log_tail = -mean+minimum*math.log(mean)-math.lgamma(minimum+1)
    log_first = math.log(record['box_count'])+log_tail
    log_mu = regions*log_falling(window, minimum)-(regions-1)*log_falling(rows, minimum)-math.lgamma(minimum+1)
    return dict(**record, occupancy_mean=mean, poisson_log2_tail=log_tail/math.log(2),
                poisson_log2_first=log_first/math.log(2),
                log2_first_factorial_moment=(math.log(record['box_count'])+log_mu)/math.log(2),
                rigorous=False)


def joint_window_assignments(window, overlap, degree, common):
    """Exact injective assignments for two row subsets into two windows."""
    rest = degree-common
    if common > overlap:
        return 0
    return falling(overlap, common)*sum((-1)**j*choose(rest, j)**2*math.factorial(j)
                                         *falling(overlap-common, j)
                                         *falling(window-common-j, rest-j)**2
                                         for j in range(rest+1))


def cyclic_power_zero(values, power):
    """Constant coefficient of a cyclic polynomial power, exactly."""
    size = len(values)
    current = [1]+[0]*(size-1)
    for _ in range(power):
        next_values = [0]*size
        for a, left in enumerate(current):
            for b, right in enumerate(values):
                next_values[(a+b) % size] += left*right
        current = next_values
    return current[0]


def mixed_pair_upper(record, degree, exact_overlap=True):
    """Exact rational BOX-pair upper bound for arbitrary equal-length windows.

    Different windows use exact inclusion-exclusion counts by default.
    The optional positive relaxation v**c W**(2(m-c)) is easier to audit.
    The suffix case adds one fixed shared region.
    """
    rows, window, minimum, starts = (record[key] for key in ('rows', 'window', 'minimum', 'starts'))
    if not 1 <= degree <= minimum or 2*degree > rows:
        raise ValueError('invalid mixed moment degree')
    free, count, suffix = record['free_regions'], len(starts), int(record['suffix'])
    outer_free = record['length']-suffix
    denominator = falling(rows, 2*degree)
    one = falling(window, degree)*(denominator//falling(rows, degree))
    intersections = [max(0, window-abs(a-b)) for a in starts for b in starts if a != b]
    multiplicities = {overlap: intersections.count(overlap) for overlap in set(intersections)}
    checksum = record.get('checksum', False)
    row_terms = []
    for common in range(degree+1):
        assignments = choose(rows, degree)*choose(degree, common)*choose(rows-degree, degree-common)
        scale = denominator//falling(rows, 2*degree-common)
        same = falling(window, 2*degree-common)*scale
        choices = {overlap: (falling(window, degree)**2 if common == 0 else 0) if overlap == 0
                   else (joint_window_assignments(window, overlap, degree, common) if exact_overlap
                         else overlap**common*window**(2*(degree-common)))
                   for overlap in multiplicities}
        different = sum(multiplicities[overlap]*value for overlap, value in choices.items())*scale
        checked = 0
        if checksum:
            # Both window-index sums are zero modulo count. For each
            # cyclic difference sequence summing to zero, exactly
            # count**(free-1) first-index sequences have checksum zero.
            # Maximize each difference's weight over its first index.
            weights = [same]+[max(choices[max(0, window-abs(starts[i]-starts[(i+delta) % count]))]
                                  for i in range(count))*scale for delta in range(1, count)]
            checked = count**(free-1)*(cyclic_power_zero(weights, free)-same**free)
        row_terms.append((assignments, count*same+different, count*same, same, checked))
    entries = []
    for shared in range(free+1):
        coordinate_pairs = choose(outer_free, free)*choose(free, shared)*choose(outer_free-free, free-shared)
        exclusive = 2*(free-shared)
        row_sum = sum(assignments*same**suffix*
                      (checked if checksum and shared == free else mixed**shared-(identical**shared if shared == free else 0))
                      for assignments, mixed, identical, same, checked in row_terms)
        window_choices = count**(exclusive-2*int(checksum)) if shared < free else 1
        entries.append(Fraction(coordinate_pairs*window_choices*one**exclusive*row_sum,
                                2*choose(minimum, degree)**2*denominator**(exclusive+shared+suffix)))
    return dict(pair_upper=sum(entries), pair_by_shared_free_regions=entries, degree=degree,
                exact_overlap=exact_overlap)


def exact_check(record, degrees=None, extra=None):
    """Finite necessary gate, with both a one-box and a union lower bound."""
    record = {key: value for key, value in record.items() if key != 'rigorous'}
    rows, window, regions, minimum = (record[key] for key in ('rows', 'window', 'regions', 'minimum'))
    mean = rows*(window/rows)**regions
    if extra is None:
        # More than enough for the small-mean cases; exact arithmetic
        # handles the cancellation, and no numeric approximation decides.
        extra = max(31, 2*math.ceil(4*mean+20)+1)
    if extra % 2 != 1:
        raise ValueError('odd extra required for lower-tail truncation')
    tail_lower = max(Fraction(0), tail_partial(rows, window, regions, minimum, extra))
    tail_upper = min(Fraction(1), tail_partial(rows, window, regions, minimum, extra+1))
    if degrees is None:
        degrees = sorted({max(1, minimum-offset) for offset in (0, 4, 8, 12, 16)})
    degrees = list(degrees)
    pair_records = [mixed_pair_upper(record, degree) for degree in degrees]
    # Each overlap class is a sum of intersection probabilities. Taking
    # the best valid bound per class remains a valid summed upper bound.
    best_pairs = [min(entry['pair_by_shared_free_regions'][shared] for entry in pair_records)
                  for shared in range(record['free_regions']+1)]
    pair_upper = sum(best_pairs)
    first_lower, first_upper = record['box_count']*tail_lower, record['box_count']*tail_upper
    union_lower = max(Fraction(0), first_lower-pair_upper)
    lower = max(tail_lower, union_lower)
    return dict(**record, rigorous=True, scope='specified-window-event-family-only',
                extra=extra, moment_degrees=degrees, tail_lower=tail_lower, tail_upper=tail_upper,
                first_lower=first_lower, first_upper=first_upper, pair_upper=pair_upper,
                pair_by_shared_free_regions=best_pairs, union_lower=union_lower,
                failure_lower=lower, rejects_40_bits=lower > Fraction(1, 1 << 40))


def sliding_event_upper(record):
    """All positions of equal-length contiguous windows, not all bad words.

    For a fixed a-subset of row labels, the probability that its positions
    fit into some (not necessarily aligned) length-W interval is exact.
    Unioning row subsets upper bounds this structural event only. It
    does not certify distance and does not cover unequal window lengths.
    """
    rows, window, minimum, regions = (record[key] for key in ('rows', 'window', 'minimum', 'regions'))
    denominator = choose(rows, minimum)
    fitting_subsets = ((rows-window+1)*choose(window-1, minimum-1)
                       + choose(window-1, minimum))
    sliding = Fraction(fitting_subsets, denominator)
    fixed = Fraction(choose(window, minimum), denominator)
    suffix = int(record['suffix'])
    upper = (choose(record['length']-suffix, regions-suffix)*denominator
             *sliding**(regions-suffix)*fixed**suffix)
    return dict(sliding_probability=sliding, event_union_upper=min(Fraction(1), upper),
                unclipped_first_moment=upper,
                scope='any-position-equal-length-window-kernel-event-only')


def all_support_sliding_upper(k_log, length, dimension, state):
    """Union of the stated sliding-window kernel events over support sizes."""
    entries = []
    for regions in range(length-dimension+1, length+1):
        for suffix in (False, True):
            record = geometry(k_log, length, dimension, state, regions, 1, suffix)
            entries.append(dict(regions=regions, suffix=suffix, minimum=record['minimum'],
                                **sliding_event_upper(record)))
    return dict(k_log=k_log, outer=[length, dimension], state=state, entries=entries,
                event_union_upper=min(Fraction(1), sum(entry['event_union_upper'] for entry in entries)),
                scope='all-supports-any-position-equal-length-window-kernel-events-only')


def unequal_window_upper(record):
    """All starts and unequal spans with total length at most threshold.

    For a fixed a-subset, span r has (L-r+1)*C(r-2,a-2)
    placements. Replacing the first factor by L gives generating function
    L*z**a/(1-z)**(a-1). A suffix extent has exact generating function
    z**a/(1-z)**a. Their positive coefficients bound the whole budget.
    """
    rows, h, a, threshold = (record[key] for key in ('rows', 'regions', 'minimum', 'threshold'))
    if a < 2:
        raise ValueError('span generating function requires at least two rows')
    suffix = int(record['suffix'])
    free = h-suffix
    numerator = (choose(record['length']-suffix, free)*rows**free
                 *choose(threshold-h+suffix, h*(a-1)+suffix))
    upper = Fraction(numerator, choose(rows, a)**(h-1))
    return dict(event_union_upper=min(Fraction(1), upper), unclipped_first_moment=upper,
                scope='all-starts-total-length-unequal-window-kernel-event-only')


def all_support_unequal_upper(k_log, length, dimension, state):
    entries = []
    for regions in range(length-dimension+1, length+1):
        for suffix in (False, True):
            record = geometry(k_log, length, dimension, state, regions, 1, suffix)
            entries.append(dict(regions=regions, suffix=suffix, minimum=record['minimum'],
                                **unequal_window_upper(record)))
    return dict(k_log=k_log, outer=[length, dimension], state=state, entries=entries,
                event_union_upper=min(Fraction(1), sum(entry['event_union_upper'] for entry in entries)),
                scope='all-supports-all-starts-total-length-unequal-window-kernel-events-only')


def screen(k_log, length, dimension, state, windows=(1, 2, 3, 4, 5), suffixes=(False, True)):
    records = []
    for regions in range(length-dimension+1, length+1):
        for count in windows:
            for suffix in suffixes:
                records.append(approximate(geometry(k_log, length, dimension, state, regions, count, suffix)))
    return sorted(records, key=lambda item: item['poisson_log2_first'], reverse=True)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--outer', type=int, nargs=2, default=[33, 24])
    parser.add_argument('--state', type=int, default=4)
    parser.add_argument('--k-log', type=int, nargs='+', default=[16, 18, 20])
    parser.add_argument('--top', type=int, default=3)
    parser.add_argument('--exact', action='store_true')
    parser.add_argument('--sliding-upper', action='store_true')
    parser.add_argument('--unequal-upper', action='store_true')
    args = parser.parse_args()
    for k_log in args.k_log:
        if args.sliding_upper or args.unequal_upper:
            method = all_support_unequal_upper if args.unequal_upper else all_support_sliding_upper
            result = method(k_log, *args.outer, args.state)
            result['entries'] = sorted(result['entries'], key=lambda entry: entry['unclipped_first_moment'], reverse=True)[:args.top]
            print(json.dumps(display(result)), flush=True)
            continue
        for record in screen(k_log, *args.outer, args.state)[:args.top]:
            print(json.dumps(display(exact_check(record) if args.exact else record)), flush=True)
