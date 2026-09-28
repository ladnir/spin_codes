"""Complete shape maxima for the unsplit lazy density column.

For each occupancy j, consider every sorted tuple in {1,2,3,4}^j.
For each shape, intersect available bounds on the same contribution,
multiply by rho^(number of fours) a^(total weight), and only then take
the maximum across shapes. This is an upper bound for every shape, not
an average under a proposed input law. The outer verifier supplies the
reciprocal weights. No odd-weight penalty is supported here.

Feedback records already contain the lazy probability 2^(-rounds).
The conditioned-column bound does not; its uniform-class version also
needs division by that class's size. Only C->C and U_v->C are changed.
"""
from fractions import Fraction as Q
from itertools import combinations_with_replacement
from numbers import Integral

from flint import arb

from column_density import density_caps
from group_rank_one_verify import up
from occupancy_memory import M, C, U
from mature_tail import L48, L56


def shapes(j):
    return combinations_with_replacement(range(1, 5), j)


def _arb(value):
    value = Q(value)
    return arb(value.numerator)/value.denominator


def _feedback_records(data, tilt, spectrum, rounds):
    if data is None:
        return {}, 0
    provenance = data['provenance']
    maximum = provenance['maximum']
    if (not isinstance(maximum, int) or not 1 <= maximum <= 6
            or provenance['rounds'] != rounds
            or provenance['spectrum'] != spectrum):
        raise ValueError('feedback census must match rounds and expansion spectrum')
    records = data['bounds'][str(tilt)]
    expected = {weights for j in range(1, maximum+1) for weights in shapes(j)}
    if set(records) != expected or provenance['checked_shapes'] != len(expected):
        raise ValueError('complete feedback shape census required')
    for record in records.values():
        if set(record['uniform']) != set(spectrum):
            raise ValueError('all uniform expansion classes are required')
        for value in (record['density'], *record['uniform'].values()):
            if not value.is_finite() or not value >= 0:
                raise ValueError('finite nonnegative feedback coefficients required')
    return records, maximum


def coefficients(maximum, spectrum, tilt, penalty='1', input_penalty='1', *,
                 window_averages=None, feedback=None, rounds=2):
    """Outward (C,U_48,...,U_80)->C coefficients for covered occupancies.

One-window records and feedback data must be computed for these same
maps and tilt. Feedback data come from feedback_density.build(), never
from a partial list or a floating-point diagnostic snapshot.
"""
    if not isinstance(maximum, int) or not 1 <= maximum <= 32:
        raise ValueError('one through 32 distinct packet windows required')
    if not isinstance(rounds, int) or isinstance(rounds, bool) or rounds < 1:
        raise ValueError('positive integer transvection count required')
    sizes = [spectrum[v] for v in sorted(spectrum)]
    if len(sizes) != 5 or any(not isinstance(h, Integral) or isinstance(h, bool) or h < 1 for h in sizes):
        raise ValueError('five positive expansion-class sizes required')
    if Q(tilt) < 0 or not 0 < Q(penalty) <= 1 or Q(input_penalty) <= 0:
        raise ValueError('nonnegative tilt, 0 < rho <= 1, and a > 0 required')
    records, feedback_maximum = _feedback_records(feedback, tilt, spectrum, rounds)
    caps = density_caps(window_averages) if window_averages is not None else None
    if caps is not None and set(caps) != {1, 2, 3, 4}:
        raise ValueError('all four one-window weight classes are required')
    lam, rho, a = _arb(tilt), _arb(penalty), _arb(input_penalty)
    lazy = arb(2)**(-rounds)
    # Shared powers avoid exponentials in the complete shape loop.
    powers = [up((lam*w).exp()) for w in range(4*maximum+1)] if caps else []
    cap_values = {b:up(_arb(cap)) for b, cap in (caps or {}).items()}
    scales = [[up(rho**f*a**w) for w in range(4*maximum+1)] for f in range(maximum+1)]
    result = {}
    for j in range(1, maximum+1):
        use_column = caps is not None and j >= 2
        if not use_column and j > feedback_maximum:
            continue
        upper = [arb(0) for _ in range(6)]
        factor = arb(32)/(33-j)
        for weights in shapes(j):
            weight = sum(weights)
            candidates = None
            if use_column:
                column = min(arb(1), *(up(factor*powers[weight-b]*cap_values[b]) for b in set(weights)))
                column = up(lazy*column)
                candidates = [column, *(up(column/h) for h in sizes)]
            if j <= feedback_maximum:
                record = records[weights]
                exact = [up(record['density']), *(up(record['uniform'][v]) for v in sorted(spectrum))]
                candidates = exact if candidates is None else [min(x, y) for x, y in zip(candidates, exact)]
            scale = scales[weights.count(4)][weight]
            for i, value in enumerate(candidates):
                upper[i] = max(upper[i], up(value*scale))
        result[j] = tuple(upper)
    return result


def refine(base, spectrum, tilt, penalty='1', input_penalty='1', *,
           window_averages=None, feedback=None, rounds=2):
    """Copy a two-update envelope and intersect complete density bounds."""
    if not 2 <= len(base) <= 33 or any(t.nrows() != 11 or t.ncols() != 11 for t in base):
        raise ValueError('eleven-coordinate epoch operators through at most 32 packets required')
    if any(t[source, C] != 0 for t in base for source in (M, L48, L56)):
        raise ValueError('the mature density column must be unsplit')
    candidates = coefficients(len(base)-1, spectrum, tilt, penalty, input_penalty,
                              window_averages=window_averages, feedback=feedback, rounds=rounds)
    result = [t*1 for t in base]
    for j, values in candidates.items():
        for source, value in zip((C, *range(U, U+5)), values):
            result[j][source, C] = min(result[j][source, C], value)
    return result
