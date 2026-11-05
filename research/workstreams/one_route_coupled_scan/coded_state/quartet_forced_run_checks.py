"""Exact bounded checks of the whole-run packing and amplitude interfaces.

Research verification only. Run file-free with Python -B on the designated
numerical host. No production probability certificate is computed here.
"""

from collections import Counter
from fractions import Fraction
from itertools import combinations_with_replacement, product


def runs(mask, length):
    result = []
    start = 0
    while start < length:
        if not (mask >> start) & 1:
            start += 1
            continue
        end = start + 1
        while end < length and (mask >> end) & 1:
            end += 1
        result.append((start, end-start, end == length))
        start = end
    return result


def subsets(mask):
    sub = mask
    while True:
        yield sub
        if not sub:
            break
        sub = (sub-1) & mask


def deaths(mask, length):
    return (mask & ~(mask >> 1) & ((1 << (length-1))-1)).bit_count()


def endpoints(parts, budget):
    fc = [ell for ell, forced, opened in parts if forced and not opened]
    fo = [ell for ell, forced, opened in parts if forced and opened]
    uc = [ell for ell, forced, opened in parts if not forced and not opened]
    uo = [ell for ell, forced, opened in parts if not forced and opened]
    remaining = budget-len(fc)
    if remaining < 0:
        return None
    low = len(fc)+sum(fo)-sum(sorted((ell-1 for ell in fo), reverse=True)[:remaining])
    high = sum(fc)+sum(fo)+sum(uo)+sum(sorted(uc, reverse=True)[:remaining])
    return low, high


def support_polynomial(parts, z, v, e):
    value = Fraction(1)
    for ell, forced, opened in parts:
        row = (Fraction(0), z) if forced else (Fraction(1), Fraction(0))
        for _ in range(ell-int(forced)):
            row = (row[0]+v*row[1], z*(row[0]+row[1]))
        value *= row[0]+(e if opened else v)*row[1]
    return value


def checks():
    families = interval_checks = inequalities = polynomials = amplitude_checks = 0
    multi_run_families = multi_run_intervals = 0
    z, v, e = Fraction(3, 2), Fraction(2, 3), Fraction(5, 4)
    for length in range(1, 9):
        for allowed in range(1 << length):
            rr = runs(allowed, length)
            starts = sum(1 << start for start, _, _ in rr)
            for forced_mask in subsets(starts):
                parts = [(ell, bool((forced_mask >> start) & 1), opened)
                         for start, ell, opened in rr]
                counts = Counter()
                exact_poly = Fraction(0)
                for optional in subsets(allowed & ~forced_mask):
                    support = optional | forced_mask
                    vol, d = support.bit_count(), deaths(support, length)
                    terminal = (support >> (length-1)) & 1
                    counts[(vol, d)] += 1
                    exact_poly += z**vol*v**d*e**terminal
                    holes = allowed.bit_count()-vol
                    for cap in range(1, length+1):
                        target = sum(cap if forced else min(cap, ell)
                                     for ell, forced, opened in parts if not opened)
                        assert holes+cap*d >= target
                        inequalities += 1
                assert exact_poly == support_polynomial(parts, z, v, e)
                polynomials += 1
                for budget in range(length+1):
                    observed = {vol for vol, d in counts if d <= budget}
                    bounds = endpoints(parts, budget)
                    expected = set() if bounds is None else set(range(bounds[0], bounds[1]+1))
                    assert observed == expected, (length, allowed, forced_mask, budget)
                    interval_checks += 1
                families += 1

    # Independent coordinate runs allow several open forced runs at once.
    # This checks the Top_b lower endpoint beyond a single coordinate.
    run_types = [(ell, forced, opened) for ell in range(1, 5)
                 for forced in (False, True) for opened in (False, True)]
    for count in range(1, 4):
        for parts in combinations_with_replacement(run_types, count):
            feasible = {(0, 0)}
            for ell, forced, opened in parts:
                local = {(mask.bit_count(), deaths(mask, ell)+
                          (0 if opened else ((mask >> (ell-1)) & 1)))
                         for mask in range(1 << ell) if not forced or mask & 1}
                feasible = {(l1+l2, d1+d2) for l1, d1 in feasible for l2, d2 in local}
            for budget in range(sum(ell for ell, _, _ in parts)+1):
                observed = {vol for vol, d in feasible if d <= budget}
                bounds = endpoints(parts, budget)
                expected = set() if bounds is None else set(range(bounds[0], bounds[1]+1))
                assert observed == expected, (parts, budget)
                multi_run_intervals += 1
            multi_run_families += 1

    # A strict improvement: four forced singleton runs and an unforced
    # closed length-four run. The old one-cap inequality accepts L7,D4.
    parts = [(1, True, False)]*4+[(4, False, False)]
    holes, d = 1, 4
    assert holes+d >= 5
    assert holes+2*d < 10
    assert endpoints(parts, d) == (4, 4)

    # The amplitude test independently enumerates actual GF(4) values.
    # Addition is XOR and delta=1; births have fixed, alternating nonzero
    # values. Every output is charged, with no terminal flush.
    q, x = 4, Fraction(2, 3)
    cx = 1+(q-2)*x
    vx, ex = (q-1)*x*x/cx, (q-1)*x/cx
    for length in range(1, 6):
        for allowed in range(1 << length):
            rr = runs(allowed, length)
            starts = sum(1 << start for start, _, _ in rr)
            for forced_mask in subsets(starts):
                parts = [(ell, bool((forced_mask >> start) & 1), opened)
                         for start, ell, opened in rr]
                options = [(1+(j % 2),) if (forced_mask >> j) & 1 else
                           (tuple(range(q)) if (allowed >> j) & 1 else (0,))
                           for j in range(length)]
                exact = Fraction(0)
                for values in product(*options):
                    previous = weight = 0
                    for current in values:
                        weight += bool(current ^ previous)
                        previous = current
                    exact += x**weight
                expected = support_polynomial(parts, cx, vx, ex)/(q-1)**forced_mask.bit_count()
                assert exact == expected, (length, allowed, forced_mask)
                amplitude_checks += 1
    print(dict(support_families=families, exact_volume_intervals=interval_checks,
               capped_run_inequalities=inequalities, support_polynomials=polynomials,
               fixed_birth_amplitude_counts=amplitude_checks,
               independent_run_families=multi_run_families,
               independent_run_volume_intervals=multi_run_intervals,
               strict_two_cap_example=True, free_terminal=True,
               full_code_certificate=False), flush=True)


if __name__ == '__main__':
    checks()
