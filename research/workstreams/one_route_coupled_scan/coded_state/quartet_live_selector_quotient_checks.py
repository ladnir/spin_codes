"""Exact tiny checks for QUARTET_LIVE_SELECTOR_QUOTIENT.md.

These tests use GF(2) and GF(4), integer weight polynomials, exhaustive
amplitude fibres, coupled homogeneous residual equations, and support
histories. They check the quotient identities, not the global geometry
sum or the construction's distance. No numerical experiment is launched
on import. The parent runs --checks under the standard bounded guard.
"""

import argparse
from collections import Counter
from fractions import Fraction
from itertools import combinations, product
import os
import shutil
import threading
import time


_LQ_GF4 = ((0, 0, 0, 0), (0, 1, 2, 3),
           (0, 2, 3, 1), (0, 3, 1, 2))


def _lq_mul(a, b, q):
    if q == 2:
        return a * b
    if q == 4:
        return _LQ_GF4[a][b]
    raise ValueError('The exact toys use only GF(2) or GF(4)')


def _lq_trim(poly):
    poly = list(poly)
    while len(poly) > 1 and poly[-1] == 0:
        poly.pop()
    return tuple(poly)


def _lq_convolve(first, second):
    result = [0] * (len(first) + len(second) - 1)
    for i, a in enumerate(first):
        for j, b in enumerate(second):
            result[i + j] += a * b
    return _lq_trim(result)


def _lq_power(poly, exponent):
    result = (1,)
    for _ in range(exponent):
        result = _lq_convolve(result, poly)
    return result


def _lq_eval(poly, x):
    value = Fraction(0)
    for coefficient in reversed(poly):
        value = value * x + coefficient
    return value


def quotient_run_formula(q, length, forced, closed):
    """Four exact integer coefficient formulas; forced start is fixed."""
    if q not in (2, 4) or length < 1:
        raise ValueError('Invalid tiny run')
    a, c = (1, q - 1), (1, -1)
    if not closed:
        return ((0,) + _lq_power(a, length - 1) if forced
                else _lq_power(a, length))
    exponent = length if forced else length + 1
    apow, cpow = _lq_power(a, exponent), _lq_power(c, exponent)
    numerator = [left + (-1 if forced else q - 1) * right
                 for left, right in zip(apow, cpow)]
    assert all(value % q == 0 for value in numerator)
    result = tuple(value // q for value in numerator)
    assert all(value >= 0 for value in result)
    return _lq_trim((0,) + result if forced else result)


def quotient_run_enumeration(q, delta, length, forced, closed, birth_value=1):
    """Independent enumeration via Y_t=Q_t+delta*Q_(t-1)."""
    if not 0 < delta < q or (forced and not 0 < birth_value < q):
        raise ValueError('Nonzero delta and forced value required')
    counts = [0] * (length + int(closed) + 1)
    tail = product(range(q), repeat=length - int(forced))
    for values in tail:
        states = (birth_value,) + values if forced else values
        previous, weight = 0, 0
        for state in states:
            weight += (state ^ _lq_mul(delta, previous, q)) != 0
            previous = state
        if closed:
            weight += _lq_mul(delta, previous, q) != 0
        counts[weight] += 1
    return _lq_trim(counts)


def quotient_allowed_runs(width, allowed, births):
    """Return(length,forced,closed) for coordinate runs; validate births."""
    if len(allowed) != len(births):
        raise ValueError('Allowed and birth masks must have the same length')
    end = len(allowed)
    full = (1 << width) - 1
    for t, (a, b) in enumerate(zip(allowed, births)):
        if a & ~full or b & ~a or (t and b & allowed[t - 1]):
            raise ValueError('Each declared birth must start an allowed run')
    runs = []
    for coordinate in range(width):
        bit, t = 1 << coordinate, 0
        while t < end:
            if not allowed[t] & bit:
                t += 1
                continue
            start = t
            while t < end and allowed[t] & bit:
                if t > start and births[t] & bit:
                    raise ValueError('A run cannot contain an interior birth')
                t += 1
            runs.append((t - start, bool(births[start] & bit), t < end))
    assert sum(forced for _, forced, _ in runs) == sum(b.bit_count() for b in births)
    return runs


def quotient_fibre_histograms(width, allowed, births, q, delta):
    """Enumerate whole allowed histories, grouped by their nonzero B values."""
    times = len(allowed)
    cells = [(t, i) for t, mask in enumerate(allowed)
             for i in range(width) if mask >> i & 1]
    birth_indices = [k for k, (t, i) in enumerate(cells) if births[t] >> i & 1]
    histograms = {}
    for values in product(range(q), repeat=len(cells)):
        b = tuple(values[k] for k in birth_indices)
        if not all(b):
            continue
        states = [[0] * width for _ in range(times)]
        for value, (t, i) in zip(values, cells):
            states[t][i] = value
        weight = 0
        previous = [0] * width
        for state in states:
            weight += sum((value ^ _lq_mul(delta, old, q)) != 0
                          for value, old in zip(state, previous))
            previous = state
        histograms.setdefault(b, Counter())[weight] += 1
    return {b: _lq_trim([histogram.get(w, 0) for w in range(width * times + 1)])
            for b, histogram in histograms.items()}


def _lq_residual_zero(matrix, b, q):
    for row in matrix:
        value = 0
        for coefficient, coordinate in zip(row, b):
            value ^= _lq_mul(coefficient, coordinate, q)
        if value:
            return False
    return True


def _lq_rank(matrix, columns, q):
    rows = [list(row) for row in matrix]
    rank = 0
    for column in range(columns):
        pivot = next((i for i in range(rank, len(rows)) if rows[i][column]), None)
        if pivot is None:
            continue
        rows[rank], rows[pivot] = rows[pivot], rows[rank]
        inverse = next(value for value in range(1, q)
                       if _lq_mul(value, rows[rank][column], q) == 1)
        rows[rank] = [_lq_mul(inverse, value, q) for value in rows[rank]]
        for i in range(len(rows)):
            if i != rank and rows[i][column]:
                multiplier = rows[i][column]
                rows[i] = [a ^ _lq_mul(multiplier, b, q)
                           for a, b in zip(rows[i], rows[rank])]
        rank += 1
    return rank


def _lq_gamma(mask, partition):
    union = 0
    for quartet in partition:
        if mask & quartet:
            union |= quartet
    return union


def quotient_support_containment(width, supports, partitions):
    """Check actual S ->(U,B,A), free terminal, and selector-run identities."""
    times, full = len(supports), (1 << width) - 1
    live = [_lq_gamma(0 if t == 0 else supports[t - 1], partitions[t])
            for t in range(times)] + [full]
    births = [s & ~live[t] for t, s in enumerate(supports)]
    allowed = [(live[t] & live[t + 1]) | births[t] for t in range(times)]
    assert live[0] == 0 and live[-1] == full
    for t, support in enumerate(supports):
        assert support & ~allowed[t] == 0
        assert births[t] & ~(live[t + 1] & ~live[t]) == 0
        assert ((allowed[t - 1] if t else 0) & ~live[t]) == 0
        assert live[t] | support == live[t] | births[t]
    f = sum(mask.bit_count() for mask in births)
    assert bool(f) == any(supports)
    runs = quotient_allowed_runs(width, allowed, births)
    drops = sum((live[t] & ~live[t + 1]).bit_count() for t in range(times))
    G = sum(mask.bit_count() for mask in live[:-1])
    assert sum(mask.bit_count() for mask in allowed) == G + f - drops
    closed_from_live = 0
    for coordinate in range(width):
        bit, t = 1 << coordinate, 0
        while t <= times:
            if not live[t] & bit:
                t += 1
                continue
            start = t
            while t <= times and live[t] & bit:
                t += 1
            if t <= times:  # Closed U run; the formal final U is never closed.
                assert start > 0
                closed_from_live += t - start >= 2 or bool(births[start - 1] & bit)
    assert sum(closed for _, _, closed in runs) == closed_from_live
    return len(runs)


def live_selector_quotient_checks():
    run_checks = envelope_checks = fibre_checks = residual_checks = containment_checks = 0
    for q in (2, 4):
        for delta in range(1, q):
            for length in range(1, 6):
                for forced, closed in product((False, True), repeat=2):
                    formula = quotient_run_formula(q, length, forced, closed)
                    for birth in (range(1, q) if forced else (1,)):
                        assert quotient_run_enumeration(q, delta, length, forced, closed, birth) == formula
                        run_checks += 1
                    assert sum(formula) == q ** (length - int(forced))
                    for x in (Fraction(1, 4), Fraction(1, 2), Fraction(1)):
                        a = 1 + (q - 1) * x
                        d = (1 + (q - 1) * x * x) / a
                        assert d >= a / q
                        envelope = a ** length * (x / a if forced else 1) * (d if closed else 1)
                        assert _lq_eval(formula, x) <= envelope
                        envelope_checks += 1
        assert quotient_run_formula(q, 1, True, False) == (0, 1)
        assert quotient_run_formula(q, 1, True, True) == (0, 0, 1)
        assert quotient_run_formula(q, 1, False, True) == (1, 0, q - 1)
    cases = ((2, (3, 1, 2, 2), (1, 0, 2, 0)),
             (2, (1, 2, 3, 1, 1), (1, 2, 0, 0, 0)),
             (3, (5, 3, 6, 2), (5, 2, 0, 0)))
    for q in (2, 4):
        for delta in range(1, q):
            for width, allowed, births in cases:
                runs = quotient_allowed_runs(width, allowed, births)
                predicted = (1,)
                for length, forced, closed in runs:
                    predicted = _lq_convolve(predicted, quotient_run_formula(q, length, forced, closed))
                fibres = quotient_fibre_histograms(width, allowed, births, q, delta)
                f = sum(mask.bit_count() for mask in births)
                assert len(fibres) == (q - 1) ** f
                assert all(histogram == predicted for histogram in fibres.values())
                fibre_checks += len(fibres)
                matrices = ((), ((1,) + (0,) * (f - 1),), ((1,) * f,),
                            tuple(tuple(int(i == j) for j in range(f)) for i in range(f)),
                            ((1,) * f, (1,) * f))
                if q == 4 and f == 3:
                    matrices += (((1, 2, 0), (0, 1, 3)),)
                for matrix in matrices:
                    compatible = [b for b in fibres if _lq_residual_zero(matrix, b, q)]
                    rank = _lq_rank(matrix, f, q)
                    all_solutions = sum(_lq_residual_zero(matrix, b, q)
                                        for b in product(range(q), repeat=f))
                    assert all_solutions == q ** (f - rank)
                    assert len(compatible) <= all_solutions
                    weighted = [0] * len(predicted)
                    representatives = [0] * len(predicted)
                    for b in compatible:
                        for i, value in enumerate(fibres[b]):
                            weighted[i] += value
                            if b[0] == 1:
                                representatives[i] += value
                    assert tuple(weighted) == tuple(len(compatible) * value for value in predicted)
                    assert all(a == (q - 1) * b for a, b in zip(weighted, representatives))
                    if f == 2 and matrix == ((1, 1),):
                        assert len(compatible) == q - 1
                    if rank == f:
                        assert not compatible
                    residual_checks += 1
    for supports in product(range(16), repeat=3):
        quotient_support_containment(4, supports, ((15,),) * 3)
        containment_checks += 1
    partitions8 = [(first, 255 ^ first) for rest in combinations(range(1, 8), 3)
                   for first in (1 + sum(1 << i for i in rest),)]
    assert len(partitions8) == 35
    for partition in partitions8:
        for population in range(9):
            for final in range(256):
                quotient_support_containment(8, ((1 << population) - 1, final),
                                             ((15, 240), partition))
                containment_checks += 1
    # T=1: the formal last live set permits every terminal birth, without reset.
    for support in range(16):
        quotient_support_containment(4, (support,), ((15,),))
        containment_checks += 1
    # A relaxed selector can be live despite zero incoming state. Its random
    # multiplier output is then concentrated, so field rank must come first.
    for q in (2, 4):
        zero = Counter(_lq_mul(multiplier, 0, q) for multiplier in range(q))
        nonzero = Counter(_lq_mul(multiplier, 1, q) for multiplier in range(q))
        assert zero == Counter({0: q})
        assert nonzero == Counter({value: 1 for value in range(q)})
    artificial_live, artificial_births = (0, 15, 0, 15), (0, 0, 1)
    artificial_allowed = tuple((artificial_live[t] & artificial_live[t + 1])
                               | artificial_births[t] for t in range(3))
    assert artificial_allowed == (0, 0, 1)
    assert quotient_allowed_runs(4, artificial_allowed, artificial_births) == [(1, True, False)]
    assert _lq_gamma(0, (15,)) != artificial_live[1]
    print(dict(exact_run_enumerators=run_checks, run_envelope_checks=envelope_checks,
               fixed_birth_fibres=fibre_checks, residual_linear_checks=residual_checks,
               support_containment_checks=containment_checks,
               nontrivial_GF4_delta=True, free_terminal_checked=True,
               zero_incoming_negative_control=True,
               global_geometry_sum_bounded=False, full_code_certificate=False), flush=True)


def live_selector_quotient_main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--checks', action='store_true', required=True)
    parser.parse_args()
    import resource
    if shutil.disk_usage('/tmp').free < 20 * (1 << 30):
        raise RuntimeError('20GiB /tmp free-space guard')
    resource.setrlimit(resource.RLIMIT_AS, (2 * (1 << 30), 2 * (1 << 30)))
    timer = threading.Timer(175, lambda: os._exit(124))
    timer.daemon = True
    timer.start()
    started = time.monotonic()
    try:
        live_selector_quotient_checks()
        print(dict(seconds=round(time.monotonic() - started, 3)), flush=True)
    finally:
        timer.cancel()


if __name__ == '__main__':
    live_selector_quotient_main()
