"""Exact regional-window failure lower bound; no samples or native timings.

Independent regional permutations act on L row labels. A box chooses h
regions and one of three equal-length windows per region. Its occupancy
counts row labels routed into all h chosen windows. Occupancy at least a
provides a nonzero shortened-source kernel when a*(h-parity)>state*h.

All decisions use integers and Fractions. Decimal/log output is display
only. The second Bonferroni term uses normalized factorial moments of
degree m, not the much looser degree-a row-subset count.
"""

from fractions import Fraction
import json
import math


def falling(n, k):
    if not 0 <= k <= n:
        return 0
    return math.prod(range(n-k+1, n+1))


def choose(n, k):
    return math.comb(n, k) if 0 <= k <= n else 0


def tail_partial(rows, window, regions, minimum, extra):
    """Alternating occupancy tail: odd extra lower, even extra upper.

    The occupancy X is the intersection size of `regions` independent
    uniform `window`-subsets of `rows`. Its jth binomial moment is
    (window)_j**regions / (j! * (rows)_j**(regions-1)).
    """
    if not 1 <= minimum <= window <= rows or regions < 1 or extra < 0:
        raise ValueError('invalid occupancy parameters')
    stop = min(window, minimum+extra)
    full_falling, full_factorial = falling(rows, stop), math.factorial(stop)
    numerator = 0
    for degree in range(minimum, stop+1):
        term = (choose(degree-1, minimum-1)*falling(window, degree)**regions
                * (full_falling//falling(rows, degree))**(regions-1)
                * (full_factorial//math.factorial(degree)))
        numerator += term if (degree-minimum) % 2 == 0 else -term
    return Fraction(numerator, full_falling**(regions-1)*full_factorial)


def shared_region_numerators(rows, window, overlap, degree, common):
    """Integer numerators over (rows)_(2*degree-common).

    Same/disjoint window probabilities are exact. The adjacent-window
    probability is upper bounded by allowing position collisions: each
    common row has at most `overlap` choices; the remaining rows have at
    most `window` choices. This deliberate relaxation avoids alternating
    inclusion-exclusion in the second-term upper bound.
    """
    if not 0 <= common <= degree or 2*degree-common > rows:
        raise ValueError('invalid row-subset geometry')
    same = falling(window, 2*degree-common)
    adjacent_upper = overlap**common*window**(2*(degree-common))
    disjoint = falling(window, degree)**2 if common == 0 else 0
    return same, adjacent_upper, disjoint


def pair_upper(rows, window, overlap, length, regions, minimum, degree):
    """Sum of upper bounds for all unordered distinct box intersections.

    Three windows have adjacent intersection `overlap` and disjoint end
    windows. Shared-region ordered choices comprise 3 same, 4 adjacent,
    and 2 disjoint pairs. Degree may be less than minimum: divide the
    mixed binomial moment by binom(minimum, degree)**2.
    """
    if not (1 <= degree <= minimum <= window <= rows and 1 <= regions <= length
            and 2*degree <= rows and 0 <= overlap <= window):
        raise ValueError('invalid box-pair parameters')
    denominator = falling(rows, 2*degree)
    one_numerator = falling(window, degree)*(denominator//falling(rows, degree))
    row_terms = []
    for common in range(degree+1):
        multiplicity = choose(rows, degree)*choose(degree, common)*choose(rows-degree, degree-common)
        scale = denominator//falling(rows, 2*degree-common)
        same, adjacent, disjoint = shared_region_numerators(rows, window, overlap, degree, common)
        row_terms.append((multiplicity, scale*(3*same+4*adjacent+2*disjoint), scale*3*same))
    entries = []
    for shared in range(regions+1):
        coordinate_pairs = choose(length, regions)*choose(regions, shared)*choose(length-regions, regions-shared)
        if not coordinate_pairs:
            entries.append(Fraction(0))
            continue
        # When every region is shared, remove identical boxes exactly.
        row_sum = sum(count*(mixed**shared-(same**shared if shared == regions else 0))
                      for count, mixed, same in row_terms)
        exclusive = 2*(regions-shared)
        entries.append(Fraction(coordinate_pairs*3**exclusive*one_numerator**exclusive*row_sum,
                                2*choose(minimum, degree)**2*denominator**(exclusive+shared)))
    return dict(total=sum(entries), by_shared_regions=entries, moment_degree=degree)


def window_box_lower(extra=31, moment_degree=29):
    """Canonical K20 state-four obstruction, with exact rational verdict."""
    rows, length, dimension, state, width = 43712, 32, 24, 4, 16
    regions, minimum, window = 9, 37, 15520
    starts = (0, 14096, 28192)
    overlap = window-starts[1]
    threshold = (rows*length+9)//10-1
    kernel_dimension = minimum*(regions-(length-dimension))-state*regions
    if extra % 2 != 1:
        raise ValueError('odd truncation offset required for a lower bound')
    assert kernel_dimension > 0 and regions*window <= threshold
    assert all(start % width == 0 for start in starts) and window % width == 0
    assert starts[-1]+window == rows and starts[2] >= window
    assert starts[1]-starts[0] == starts[2]-starts[1]
    boxes = choose(length, regions)*3**regions
    tail_lower = max(Fraction(0), tail_partial(rows, window, regions, minimum, extra))
    tail_upper = tail_partial(rows, window, regions, minimum, extra+1)
    pairs = pair_upper(rows, window, overlap, length, regions, minimum, moment_degree)
    first_lower = boxes*tail_lower
    lower = max(Fraction(0), first_lower-pairs['total'])
    return dict(rows=rows, k=rows*dimension, n=rows*length, outer=[length, dimension],
                state=state, selected_regions=regions, minimum_occupancy=minimum,
                window=window, starts=starts, adjacent_overlap=overlap, box_count=boxes,
                total_window=regions*window, bad_weight_threshold=threshold,
                kernel_dimension_lower=kernel_dimension, truncation_extra=extra,
                tail_lower=tail_lower, tail_upper=tail_upper,
                first_bonferroni_lower=first_lower, second_bonferroni_upper=pairs['total'],
                second_by_shared_regions=pairs['by_shared_regions'], moment_degree=moment_degree,
                failure_lower=lower, exceeds_2_minus_38=lower > Fraction(1, 1 << 38),
                exceeds_2_minus_40=lower > Fraction(1, 1 << 40),
                probability_space='independent-uniform-regional-permutations',
                holds_for='every choice in the stated invertible linear causal four-state encoder family',
                rigorous=True)


def display(value):
    if isinstance(value, Fraction):
        return dict(approx=float(value), log2=(math.log2(value.numerator)-math.log2(value.denominator)
                                               if value > 0 else None))
    if isinstance(value, dict):
        return {key: display(item) for key, item in value.items()}
    if isinstance(value, (list, tuple)):
        return [display(item) for item in value]
    return value


if __name__ == '__main__':
    print(json.dumps(display(window_box_lower()), indent=2))
