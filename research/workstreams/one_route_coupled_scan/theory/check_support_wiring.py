#!/usr/bin/env python3
"""Exact necessary support patterns for resets in the coupled scan.

For nonzero diagonal entries, Q_next=0 requires X=(R+D)Q. Write H for
the nonzero state coordinates and I for the nonzero input coordinates.
Necessarily H xor R(H) is contained in I, and I is contained in H union
R(H). This script enumerates these MASK conditions, not field values.

Natural and pair-antipodal cycles are compared at widths 2, 4, and 8.
Physical pairs are {0,1}, {2,3}, and so on. A cycle lists the successive
physical positions visited by R. These fixed local wires add no route.

The mask condition is not an algebraic sufficiency claim: H=all positions
and I=empty passes it, but an invertible R+D forbids that reset. No distance,
probability, implementation-cost, or coefficient-realizability claim is made.
Uses only the standard library, writes no files, and should run with -B.
"""

from collections import Counter, defaultdict


def require(condition, message):
    if not condition:
        raise AssertionError(message)


def rotate_mask(mask, cycle):
    result = 0
    for index, position in enumerate(cycle):
        if mask & (1 << position):
            result |= 1 << cycle[(index + 1) % len(cycle)]
    return result


def mask_condition(state, inputs, cycle):
    rotated = rotate_mask(state, cycle)
    boundary = state ^ rotated
    union = state | rotated
    return (boundary & ~inputs) == 0 and (inputs & ~union) == 0


def physical_pair_unions(width):
    """Return (selected pair mask, allowed physical-coordinate mask)."""
    result = []
    for selected in range(1 << (width // 2)):
        allowed = 0
        for pair in range(width // 2):
            if selected & (1 << pair):
                allowed |= 3 << (2 * pair)
        result.append((selected, allowed))
    return tuple(result)


def enumerate_patterns(cycle):
    width = len(cycle)
    patterns = []
    by_state = defaultdict(list)
    for state in range(1, 1 << width):
        for inputs in range(1 << width):
            if mask_condition(state, inputs, cycle):
                patterns.append((state, inputs))
                by_state[state].append(inputs)
    return tuple(patterns), by_state


def minimum_for_allowed(allowed, cycle, by_state):
    """Enumerate exact I masks, then minimize over nonempty state masks H."""
    minimum = len(cycle) + 1
    count = 0
    for state, input_masks in by_state.items():
        admissible = [inputs for inputs in input_masks if (inputs & ~allowed) == 0]
        boundary = state ^ rotate_mask(state, cycle)
        require(bool(admissible) == ((boundary & ~allowed) == 0),
                ("allowed-input relaxation mismatch", cycle, state, allowed))
        if admissible:
            minimum = min(minimum, state.bit_count())
            count += len(admissible)
    require(minimum <= len(cycle), ("full support mask disappeared", cycle, allowed))
    return minimum, count


def shortest_allowed_gap(allowed, cycle):
    """Independent arc characterization of the minimum in the mask model.

With fewer than two allowed boundary positions, only the full state mask
survives. Otherwise a shortest arc gives the minimum proper state support.
"""
    positions = [index for index, physical in enumerate(cycle)
                 if allowed & (1 << physical)]
    width = len(cycle)
    if len(positions) < 2:
        return width
    return min((positions[(i + 1) % len(positions)] - position) % width
               for i, position in enumerate(positions))


def expected_pair_histogram(width, name, selected_count):
    # Counts refer to distinct unions of physical pairs, not probabilities.
    if selected_count == 0:
        return Counter({width: 1})
    if name == "natural" or width == 2:
        counts = {2: (1, 1), 4: (1, 2, 1), 8: (1, 4, 6, 4, 1)}
        return Counter({1: counts[width][selected_count]})
    if width == 4:
        return Counter({2: 2}) if selected_count == 1 else Counter({1: 1})
    if selected_count == 1:
        return Counter({4: 4})
    if selected_count == 2:
        return Counter({1: 4, 2: 2})
    return Counter({1: 4 if selected_count == 3 else 1})


def check_cycle(width, name, cycle):
    require(sorted(cycle) == list(range(width)), ("not one full cycle", cycle))
    full = (1 << width) - 1
    patterns, by_state = enumerate_patterns(cycle)
    expected_patterns = {2: 6, 4: 46, 8: 2206}
    require(len(patterns) == expected_patterns[width],
            ("unexpected global mask count", width, name, len(patterns)))
    require(mask_condition(full, 0, cycle),
            "necessary-mask model must retain its full-support false positive")

    for state in range(1, full):
        boundary = state ^ rotate_mask(state, cycle)
        require(boundary.bit_count() >= 2 and boundary.bit_count() % 2 == 0,
                ("proper support has invalid cyclic boundary", cycle, state))

    single_coordinate_cases = 0
    for position in range(width):
        allowed = 1 << position
        minimum, _ = minimum_for_allowed(allowed, cycle, by_state)
        require(minimum == width, ("single-coordinate reset mask", cycle, position, minimum))
        require(minimum == shortest_allowed_gap(allowed, cycle),
                "single-coordinate arc check")
        single_coordinate_cases += 1

    summaries = defaultdict(Counter)
    pattern_counts = Counter()
    union_cases = 0
    for selected, allowed in physical_pair_unions(width):
        minimum, count = minimum_for_allowed(allowed, cycle, by_state)
        require(minimum == shortest_allowed_gap(allowed, cycle),
                ("independent shortest-arc check", cycle, allowed, minimum))
        selected_count = selected.bit_count()
        if selected_count == 1:
            expected_minimum = 1 if name == "natural" else width // 2
            require(minimum == expected_minimum,
                    ("one physical pair", width, name, allowed, minimum))
        summaries[selected_count][minimum] += 1
        pattern_counts[selected_count] += count
        union_cases += 1

    for selected_count in range(width // 2 + 1):
        require(summaries[selected_count] == expected_pair_histogram(width, name, selected_count),
                ("physical-pair union histogram", width, name, selected_count,
                 summaries[selected_count]))

    return {
        "width": width,
        "name": name,
        "cycle": cycle,
        "global_patterns": len(patterns),
        "single_coordinate_cases": single_coordinate_cases,
        "pair_union_cases": union_cases,
        "summaries": summaries,
        "pattern_counts": pattern_counts,
    }


def main():
    results = []
    for width in (2, 4, 8):
        cycles = {
            "natural": tuple(range(width)),
            "pair-antipodal": tuple(range(0, width, 2)) + tuple(range(1, width, 2)),
        }
        for name, cycle in cycles.items():
            results.append(check_cycle(width, name, cycle))

    require(sum(item["single_coordinate_cases"] for item in results) == 28,
            "unexpected single-coordinate context count")
    require(sum(item["pair_union_cases"] for item in results) == 44,
            "unexpected physical-pair union context count")
    print("PASS: exact necessary support masks; no coefficient-realizability or distance claim")
    print("single-coordinate contexts=28; physical-pair union contexts=44")
    for item in results:
        print(f"w={item['width']} {item['name']} cycle={item['cycle']}; "
              f"global (H,I) patterns={item['global_patterns']}")
        for selected_count, histogram in sorted(item["summaries"].items()):
            formatted = ", ".join(f"{minimum}:{count}" for minimum, count in sorted(histogram.items()))
            note = " [empty input: full-support false positive only]" if selected_count == 0 else ""
            print(f"  allowed pairs={selected_count}; min-|H| histogram {{{formatted}}}; "
                  f"(allowed,H,I) patterns={item['pattern_counts'][selected_count]}{note}")


if __name__ == "__main__":
    main()
