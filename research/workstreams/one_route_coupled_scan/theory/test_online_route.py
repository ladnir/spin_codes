"""Small exact route-law and local-operator checks; no numerical certificate.

Enumerates at most six pair slots and two labeled active pairs. Generic
linear operators test chronological without-replacement averaging. This
does not assert that a nonlinear Bellman operator commutes with averaging.
No files, large-field enumeration, or benchmark; run serially with -B.
"""

from collections import Counter
from fractions import Fraction as F
from functools import lru_cache
from itertools import permutations, product

from check_endpoint_tail import active_terminals
from online_route import VALID, active, combine, counts, region, zeros


def require(condition, detail):
    if not condition:
        raise AssertionError(detail)


MASKS = ((0,), (1, 2, 4, 8), (5, 6, 9, 10))
TERMINAL = (F(2, 3), F(7, 5))


def matrix(mask):
    # Deliberately noncommuting positive operators, unrelated to the proof.
    return ((F(2 + mask, 17), F(1 + (mask & 3), 19)),
            (F(1 + mask.bit_count(), 23), F(3 + mask * mask % 11, 29)))


def apply(mask, value):
    operator = matrix(mask)
    return tuple(sum((operator[i][j] * value[j] for j in range(2)), F(0))
                 for i in range(2))


def average_apply(masks, value):
    children = [apply(mask, value) for mask in masks]
    return tuple(sum((child[i] for child in children), F(0)) / len(children)
                 for i in range(2))


@lru_cache(maxsize=None)
def linear_route(steps, remaining):
    if not steps:
        require(remaining == 0, "unplaced active pair")
        return TERMINAL
    result = [F(0), F(0)]
    for consumed, numerator, denominator in counts(steps, remaining):
        child = average_apply(MASKS[consumed], linear_route(steps - 1, remaining - consumed))
        for i in range(2):
            result[i] += F(numerator, denominator) * child[i]
    return tuple(result)


@lru_cache(maxsize=None)
def linear_full_pair(steps):
    empty_suffix = TERMINAL
    for _ in range(steps - 1):
        empty_suffix = apply(0, empty_suffix)
    placed = average_apply((3, 12), empty_suffix)
    if steps == 1:
        return placed
    later = apply(0, linear_full_pair(steps - 1))
    return tuple(placed[i] / steps + F(steps - 1, steps) * later[i] for i in range(2))


def evaluate_sequence(masks):
    value = TERMINAL
    for mask in reversed(masks):
        value = apply(mask, value)
    return value


def check_routes():
    laws = samples = 0
    for steps in (1, 2, 3):
        for remaining in (0, 1, 2):
            first = Counter()
            accumulated = [F(0), F(0)]
            total = 0
            for slots in permutations(range(2 * steps), remaining):
                for swaps in product((0, 1), repeat=remaining):
                    masks = [0] * steps
                    for slot, swap in zip(slots, swaps):
                        masks[slot // 2] |= 1 << (2 * (slot % 2) + swap)
                    first[masks[0]] += 1
                    value = evaluate_sequence(masks)
                    for i in range(2):
                        accumulated[i] += value[i]
                    total += 1
            predicted = {mask: F(0) for mask in range(16)}
            for consumed, numerator, denominator in counts(steps, remaining):
                for mask in MASKS[consumed]:
                    predicted[mask] += F(numerator, denominator * len(MASKS[consumed]))
            require(sum(predicted.values()) == 1, ("normalization", steps, remaining))
            for mask in range(16):
                require(F(first[mask], total) == predicted[mask],
                        ("first-batch law", steps, remaining, mask))
            require(tuple(value / total for value in accumulated) == linear_route(steps, remaining),
                    ("ordered singleton-pair recurrence", steps, remaining))
            laws += 1
            samples += total
        full_values = []
        for slot in range(2 * steps):
            masks = [0] * steps
            masks[slot // 2] = 3 << (2 * (slot % 2))
            full_values.append(evaluate_sequence(masks))
        expected = tuple(sum((value[i] for value in full_values), F(0)) / len(full_values)
                         for i in range(2))
        require(expected == linear_full_pair(steps), ("ordered full-pair recurrence", steps))
        samples += len(full_values)
        laws += 1
    require(laws == 12, "route fixture count")
    return laws, samples


def check_local_operators():
    z, tau = F(3, 5), F(1024)
    live = zeros(z)
    quiet = zeros(z)
    quiet[0][0] = F(1)
    for rank, mask in VALID[1:]:
        live[rank][mask] = F(1 + 3 * rank + 2 * mask, 11)
    comparisons = 0
    for fresh, outer_randomized in ((False, False), (True, False), (False, True)):
        mode = dict(fresh=fresh, outer_randomized=outer_randomized)
        for mask in (1, 3, 5):
            reset_expected, live_expected = active_terminals(live, mask, z, tau, **mode)
            require(active(quiet, mask, z, tau, **mode) == reset_expected,
                    ("local reset terminal", fresh, mask))
            require(active(live, mask, z, tau, **mode) == live_expected,
                    ("local live terminal", fresh, mask))
            comparisons += 2
        for terminal in (quiet, live):
            pair, split = region(terminal, z, tau, steps=1, **mode)
            expected_pair = combine(z, [(F(1, 2), active(terminal, mask, z, tau, **mode))
                                        for mask in (3, 12)])
            expected_split = combine(z, [(F(1, 4), active(terminal, mask, z, tau, **mode))
                                         for mask in (5, 6, 9, 10)])
            require(pair == expected_pair and split == expected_split, ("one-batch boundary", fresh))
            comparisons += 2
    return comparisons


def check_directed_arithmetic():
    from decimal import Decimal as D, ROUND_CEILING, localcontext
    z, tau = F(3, 5), F(1024)
    terminal = zeros(z)
    for rank, mask in VALID:
        terminal[rank][mask] = F(1 + rank + mask, 11)
    comparisons = 0
    for fresh, outer_randomized in ((False, False), (True, False), (False, True)):
        mode = dict(fresh=fresh, outer_randomized=outer_randomized)
        exact = region(terminal, z, tau, steps=3, **mode)
        for precision in (80, 120):
            with localcontext() as context:
                context.prec = precision
                context.rounding = ROUND_CEILING
                up = lambda x: D(x.numerator) / D(x.denominator)
                rounded = region([[up(x) for x in row] for row in terminal], up(z), up(tau), steps=3, **mode)
                for expected, actual in zip(exact, rounded):
                    for rank, mask in VALID:
                        require(F(actual[rank][mask]) >= expected[rank][mask],
                                ("directed comparison", mode, precision, rank, mask))
                        comparisons += 1
    return comparisons


def main():
    laws, samples = check_routes()
    operators = check_local_operators()
    directed = check_directed_arithmetic()
    print(f"PASS: {laws} exact small-route laws/linear recurrences; {samples} labeled placement/swap samples.")
    print(f"PASS: {operators} local terminal/boundary table comparisons, ordinary/fresh/lazy, tau>q.")
    print(f"PASS: {directed} directed80/120 comparisons against exact three-batch nonlinear tables.")
    print("No claim that nonlinear averaging commutes; no whole-code certificate.")


if __name__ == "__main__":
    main()
