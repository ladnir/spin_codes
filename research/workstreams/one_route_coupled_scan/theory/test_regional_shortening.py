#!/usr/bin/env python3
"""Bounded exact tests of regional_shortening.py; not a code certificate.

Exhaust all binary incidence matrices in four tiny geometries, with at most
2^16 matrices per geometry. Compare the maximum in every column-count class
against the formula and an explicit nested-support witness. Check the old
secant, complement identity, monotonicity, and virtual completions separately.
All arithmetic is integer. No dependencies, field sweeps, or output files.
The root agent owns execution; run serially with python -B.
"""

from itertools import product
from time import monotonic

from regional_shortening import (
    largest_half_regions,
    regional_shortening_cap,
    secant_shortening_cap,
    virtual_completion_cap,
)


def require(condition, message):
    if not condition:
        raise AssertionError(message)


def guard(deadline):
    if monotonic() > deadline:
        raise TimeoutError("regional shortening tests exceeded the 30-second cap")


def excess(columns, rows):
    half = len(columns) // 2
    return sum(max(0, sum((column >> i) & 1 for column in columns) - half)
               for i in range(rows))


def check_geometry(regions, rows, deadline):
    require(regions * rows <= 16, "exhaustive matrix cap exceeded")
    maxima = {}
    matrix_count = 0
    for columns in product(range(1 << rows), repeat=regions):
        counts = tuple(column.bit_count() for column in columns)
        value = excess(columns, rows)
        maxima[counts] = max(maxima.get(counts, 0), value)
        matrix_count += 1
        if matrix_count % 1024 == 0:
            guard(deadline)

    require(len(maxima) == (rows + 1) ** regions, "missing count class")
    virtual_checks = monotone_checks = 0
    for counts, maximum in maxima.items():
        cap = regional_shortening_cap(counts, rows)
        require(maximum == cap, f"exact maximum failed for {counts}")
        nested = tuple((1 << a) - 1 for a in counts)
        require(excess(nested, rows) == cap, "nested witness does not attain cap")
        require(cap <= secant_shortening_cap(counts, rows), "new cap lost to old secant")
        complement = tuple(rows - a for a in counts)
        require(cap - regional_shortening_cap(complement, rows)
                == sum(counts) - rows * regions // 2, "complement identity failed")
        largest = largest_half_regions(counts, rows)
        require(virtual_completion_cap(counts, rows, largest) == cap,
                "largest-half completion did not attain cap")

        for j, a in enumerate(counts):
            if a < rows:
                increased = list(counts)
                increased[j] += 1
                changed = regional_shortening_cap(increased, rows)
                require(cap <= changed <= cap + 1, "one-coordinate monotonicity failed")
                monotone_checks += 1

        if regions <= 4:
            best = None
            for subset in range(1 << regions):
                selected = tuple(j for j in range(regions) if (subset >> j) & 1)
                candidate = virtual_completion_cap(counts, rows, selected)
                require(candidate >= cap, "virtual completion beat exact maximum")
                best = candidate if best is None else min(best, candidate)
                virtual_checks += 1
            require(best == cap, "optimal virtual completion differs from cap")
        guard(deadline)
    return matrix_count, len(maxima), virtual_checks, monotone_checks


def check_large_count_fixtures(deadline):
    """Test large integers only; never enumerate the large row population."""
    rows, regions = 32768, 64
    seed = 0x5EC4A7
    for fixture in range(32):
        counts = []
        for j in range(regions):
            seed = (1664525 * seed + 1013904223) & 0xFFFFFFFF
            value = seed % (rows + 1)
            if fixture % 4 == 0:
                value = rows - 1 if j < regions // 2 else 1
            elif fixture % 4 == 1:
                value = rows if j < fixture else 0
            counts.append(value)
        cap = regional_shortening_cap(counts, rows)
        require(cap <= secant_shortening_cap(counts, rows), "large secant comparison failed")
        require(virtual_completion_cap(counts, rows, largest_half_regions(counts, rows)) == cap,
                "large completion comparison failed")
        require(cap - regional_shortening_cap((rows - a for a in counts), rows)
                == sum(counts) - rows * regions // 2, "large complement comparison failed")
        guard(deadline)
    # No region is full, but its exact cardinality still improves the bound.
    require(regional_shortening_cap((3, 3, 1, 1), 4) == 2,
            "near-full fixture lost its exact cap")
    require(secant_shortening_cap((3, 3, 1, 1), 4) == 4,
            "near-full fixture must separate the bounds")
    return 32


def check_rejections():
    for counts, rows in (((), 1), ((0,), 1), ((0, 0), 0), ((-1, 0), 2),
                         ((0, 3), 2), ((0, 1.5), 2), ((0, 0), 1.5)):
        for function in (regional_shortening_cap, secant_shortening_cap, largest_half_regions):
            try:
                function(counts, rows)
            except ValueError:
                pass
            else:
                raise AssertionError("invalid geometry was accepted")
    for index in (-1, 2, 0.5):
        try:
            virtual_completion_cap((0, 1), 2, (index,))
        except ValueError:
            pass
        else:
            raise AssertionError("invalid completion index was accepted")


def main():
    deadline = monotonic() + 30.0
    totals = [0, 0, 0, 0]
    for regions, rows in ((2, 4), (4, 3), (4, 4), (8, 2)):
        counts = check_geometry(regions, rows, deadline)
        totals = [a + b for a, b in zip(totals, counts)]
    large = check_large_count_fixtures(deadline)
    check_rejections()
    guard(deadline)
    print(f"PASS regional shortening: {totals[0]} binary matrices, {totals[1]} count classes, "
          f"{totals[2]} virtual completions, {totals[3]} monotonicity checks, "
          f"{large} large-integer fixtures.")
    print("Exact local combinatorial tests only; no packet-route attainability or distance certificate.")


if __name__ == "__main__":
    main()
