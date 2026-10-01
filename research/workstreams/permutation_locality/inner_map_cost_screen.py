"""Bounded, exact local screen of proposed t64/t128, s16 inner maps.

Stdlib only; prints JSON to stdout and does not write reports or map files.
This is a structural/cost screen, not a whole-code distance certificate.
The declared ordering of aligned four-coordinate packets is preserved.
"""

from collections import Counter, defaultdict
from hashlib import sha256
from itertools import combinations
import json
from math import comb
from pathlib import Path

from inner_packet_codegen import cse


HERE = Path(__file__).resolve().parent
MAPS = HERE.parent / "rate_quarter_bch/inner_calibration/maps"


def rank(words):
    pivots = {}
    for word in words:
        while word:
            bit = word.bit_length() - 1
            if bit not in pivots:
                pivots[bit] = word
                break
            word ^= pivots[bit]
    return len(pivots)


def anf(row, length):
    values = [(row >> p) & 1 for p in range(length)]
    for bit in (1 << j for j in range(length.bit_length() - 1)):
        for p in range(length):
            if p & bit:
                values[p] ^= values[p ^ bit]
    if any(value and p.bit_count() > 2 for p, value in enumerate(values)):
        raise ArithmeticError("map is not a quadratic RM submap")
    return sum(value << p for p, value in enumerate(values))


def low_dual(columns):
    """Enumerate distinct physical supports, not repeated-index tuples."""
    by_value = defaultdict(list)
    for i, column in enumerate(columns):
        by_value[column].append(i)
    pairs = defaultdict(list)
    for i, j in combinations(range(len(columns)), 2):
        pairs[columns[i] ^ columns[j]].append((i, j))
    triples = sum(1 for i, j in combinations(range(len(columns)), 2)
                  for k in by_value[columns[i] ^ columns[j]] if k > j)
    quads = set()
    for bucket in pairs.values():
        for left, right in combinations(bucket, 2):
            support = tuple(sorted(left + right))
            if len(set(support)) == 4:
                quads.add(support)
    shapes = Counter()
    for support in quads:
        occupancy = Counter(p // 4 for p in support)
        shapes[",".join(map(str, sorted(occupancy.values())))] += 1
    return {
        "weight_1": len(by_value[0]),
        "weight_2": len(pairs[0]),
        "weight_3": triples,
        "weight_4": len(quads),
        "weight_4_packet_occupancies": dict(sorted(shapes.items())),
        "first_weight_4_witness": list(min(quads)) if quads else None,
    }


def dual_count(spectrum, length, dimension, weight):
    numerator = sum(multiplicity * sum(
        (-1) ** j * comb(i, j) * comb(length - i, weight - j)
        for j in range(max(0, weight - length + i), min(weight, i) + 1))
        for i, multiplicity in spectrum.items())
    count, remainder = divmod(numerator, 1 << dimension)
    if remainder or count < 0:
        raise ArithmeticError("nonintegral or negative MacWilliams count")
    return count


def screen(name, length, rows, source=None, selected_monomials=None):
    if len(rows) != 16 or any(not 0 <= r < 1 << length for r in rows):
        raise ValueError("sixteen in-range expansion rows required")
    dimension = rank(rows)
    if dimension != 16:
        raise ArithmeticError("screen candidates must have rank sixteen")
    columns = [sum(((row >> p) & 1) << j for j, row in enumerate(rows))
               for p in range(length)]
    packet_ranks = [rank(columns[p:p + 4]) for p in range(0, length, 4)]
    pair_packet_ranks = Counter(rank(columns[4 * a:4 * a + 4] +
                                   columns[4 * b:4 * b + 4])
                               for a, b in combinations(range(length // 4), 2))
    spectrum = Counter({0: 1})
    packet_spectrum = Counter({0: 1})
    images = [0] * (1 << 16)
    for state in range(1, 1 << 16):
        low = state & -state
        image = images[state ^ low] ^ rows[low.bit_length() - 1]
        images[state] = image
        spectrum[image.bit_count()] += 1
        active = sum(bool((image >> p) & 15) for p in range(0, length, 4))
        packet_spectrum[active] += 1
    if len(set(images)) != 1 << 16:
        raise ArithmeticError("expansion enumeration is not injective")
    dual = low_dual(columns)
    for weight in range(1, 5):
        if dual[f"weight_{weight}"] != dual_count(spectrum, length, 16, weight):
            raise ArithmeticError("support enumeration and MacWilliams disagree")
    dual_minimum = next(weight for weight in range(1, length + 1)
                        if dual_count(spectrum, length, 16, weight))

    monomials = [m for m in range(length) if m.bit_count() <= 2]
    anf_rows = [anf(row, length) for row in rows]
    coefficients = [sum(((row >> m) & 1) << j for j, row in enumerate(anf_rows))
                    for m in monomials]
    feedback = [sum(((row >> m) & 1) << j for j, m in enumerate(monomials))
                for row in anf_rows]
    # These are basis-dependent straight-line XOR counts, not instructions,
    # cycle estimates, or complete costs including packet evaluation/packing.
    coefficient_cost = cse(coefficients, 16)[2]
    feedback_cost = cse(feedback, len(monomials))[2]
    return dict(
        name=name, step_bits=length, state_bits=16, expansion_rank=dimension,
        feedback_rank=rank(columns),
        feedback_times_expansion_zero=all(not ((a & b).bit_count() & 1)
                                          for a in rows for b in rows),
        zero_columns=columns.count(0), distinct_columns=len(set(columns)),
        aligned_four_packet_ranks=packet_ranks,
        single_packet_nonzero_labels_with_zero_feedback=[(1 << (4 - r)) - 1 for r in packet_ranks],
        two_packet_union_rank_histogram=dict(sorted(pair_packet_ranks.items())),
        minimum_expansion_weight=min(weight for weight in spectrum if weight),
        expansion_weight_spectrum=dict(sorted(spectrum.items())),
        minimum_active_output_packets=min(weight for weight in packet_spectrum if weight),
        active_output_packet_spectrum=dict(sorted(packet_spectrum.items())),
        minimum_feedback_kernel_weight=dual_minimum, low_weight_feedback_kernel=dual,
        native_basis_anf_coefficient_xors=coefficient_cost,
        native_basis_feedback_finish_xors=feedback_cost,
        selected_monomials=selected_monomials, source=source,
        whole_code_certificate=False,
    )


def load_declared(name):
    path = MAPS / (name + ".json")
    raw = path.read_bytes()
    record = json.loads(raw)
    rows = [int(word, 16) for word in record["generator_rows_hex"]]
    length = record["step_bits"]
    columns = [sum(((row >> p) & 1) << j for j, row in enumerate(rows))
               for p in range(length)]
    if record["state_bits"] != 16 or columns != record["columns"]:
        raise ArithmeticError("declared expansion/feedback geometry disagrees")
    result = screen(name, length, rows, dict(path=str(path), sha256=sha256(raw).hexdigest()))
    if path.read_bytes() != raw:
        raise ArithmeticError("map source changed while screening")
    return result


def monomial_map(name, variables, edges):
    monomials = [0] + [1 << j for j in range(variables)]
    monomials += [(1 << a) | (1 << b) for a, b in edges]
    if len(monomials) != 16 or len(set(monomials)) != 16:
        raise ValueError("sixteen distinct RM2 monomials required")
    length = 1 << variables
    rows = [sum(1 << p for p in range(length) if p & m == m) for m in monomials]
    return screen(name, length, rows, selected_monomials=monomials)


def selected_neighbors():
    """Bounded 9*15 single-monomial mutations; fully enumerate six finalists.

    Reject basis-only changes and all weight-four feedback cancellations
    before considering cost. This is a proposal search, not certification.
    """
    path = MAPS / "t64_s16_selected.json"
    raw = path.read_bytes()
    original = [int(word, 16) for word in json.loads(raw)["generator_rows_hex"]]
    monomials = [m for m in range(64) if m.bit_count() <= 2]
    proposals = []
    counts = Counter()
    for row_index in range(7, 16):
        for monomial in (m for m in monomials if m.bit_count() == 2):
            counts["tested"] += 1
            rows = original.copy()
            rows[row_index] ^= sum(1 << p for p in range(64) if p & monomial == monomial)
            if rank(rows) != 16:
                counts["rank_failure"] += 1
                continue
            if rank(original + rows) == 16:
                counts["basis_only_change"] += 1
                continue
            columns = [sum(((row >> p) & 1) << j for j, row in enumerate(rows)) for p in range(64)]
            if len({columns[i] ^ columns[j] for i, j in combinations(range(64), 2)}) != comb(64, 2):
                counts["weight_four_cancellation"] += 1
                continue
            if any(rank(columns[p:p + 4]) != 4 for p in range(0, 64, 4)):
                counts["packet_rank_failure"] += 1
                continue
            if any(rank(columns[4*a:4*a+4] + columns[4*b:4*b+4]) != 7
                   for a, b in combinations(range(16), 2)):
                counts["two_packet_rank_failure"] += 1
                continue
            anf_rows = [anf(row, 64) for row in rows]
            coefficients = [sum(((row >> m) & 1) << j for j, row in enumerate(anf_rows)) for m in monomials]
            feedback = [sum(((row >> m) & 1) << j for j, m in enumerate(monomials)) for row in anf_rows]
            expansion_cost, feedback_cost = cse(coefficients, 16)[2], cse(feedback, len(monomials))[2]
            counts["local_filter_pass"] += 1
            proposals.append((expansion_cost + feedback_cost, expansion_cost,
                              feedback_cost, row_index, monomial, rows))
    finalists = []
    for _, _, _, row_index, monomial, rows in sorted(proposals)[:6]:
        finalists.append(screen(
            f"t64_selected_toggle_row{row_index}_monomial{monomial}", 64, rows,
            source=dict(parent_path=str(path), parent_sha256=sha256(raw).hexdigest(),
                        modified_row=row_index, xor_monomial=monomial,
                        generator_rows_hex=[hex(row) for row in rows])))
    if path.read_bytes() != raw:
        raise ArithmeticError("selected parent changed during neighbor screen")
    return dict(counts), finalists


def main():
    results = [load_declared(name) for name in (
        "t64_s16_selected", "t64_s16_nested", "t128_s16_nested")]
    pairs6 = list(combinations(range(6), 2))
    results += [
        monomial_map("t64_monomial_lex9", 6, pairs6[:9]),
        monomial_map("t64_monomial_balanced_k33", 6,
                     [(a, b) for a in (0, 2, 4) for b in (1, 3, 5)]),
        monomial_map("t64_monomial_prism", 6,
                     [(0,1),(1,2),(0,2),(3,4),(4,5),(3,5),(0,3),(1,4),(2,5)]),
        monomial_map("t64_monomial_no_x0x1", 6, [p for p in pairs6 if p != (0, 1)][:9]),
        monomial_map("t128_monomial_cycle_chord", 7,
                     [(0,1),(1,2),(2,3),(3,4),(4,5),(5,6),(0,6),(0,3)]),
    ]
    neighbor_counts, neighbors = selected_neighbors()
    results += neighbors
    print(json.dumps(dict(
        schema="inner-map-static-cost-screen-1", whole_code_certificate=False,
        limitations=[
            "Local algebra and exact finite spectra only; no multi-step distance or security claim.",
            "Native-basis XOR CSE counts omit SIMD layout, matrix refresh, packet evaluation, routing and outer work.",
            "Changing the map or physical step length invalidates reuse of a selected-map whole-code certificate.",
            "Physical t128 is not two independent t64 updates; state refresh frequency changes.",
            "Aligned four-coordinate packets use the declared coordinate order without hidden permutation.",
        ], bounded_single_monomial_search=neighbor_counts, candidates=results), indent=2))


if __name__ == "__main__":
    main()
