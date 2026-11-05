"""Directed cell census with ratio-binned two-batch fractionalization.

For positive raw local matrices M,N, write (MN)[a,b]=x+y. Then
(x+y)^alpha=c(x/y)*(x^alpha+y^alpha), where
c(r)=(1+r)^alpha/(1+r^alpha). The ratio factors into one row ratio of M
and one column ratio of N. Exact rational comparisons assign every local
occupancy to a log2 ratio bin of width1 or1/4. The correction increases
with |log2(r)|, so the farther endpoint of a summed bin interval bounds it.
Positive Arb bin sums and precomputed correction matrices bound the IID
two-batch expectation without enumerating all4845^2 occupancy pairs.

The IID law is UNCONDITIONED. Its two-batch upper matrix W gives W^64
for128 batches. Divide once by a certified lower bound on the full4096-pair
histogram probability, not its alpha power. Endpoint minima cover complete
full-pair cells. Every ordered cell pair retains its original chronology.
An interval-maximum sweep supplies the bound for each occupied-run class.
Exact cycle-run multiplicities remain outside alpha and sum to binom(8192,g).

The route-only census divided by binom(8192,g) is an averaged raw anchor.
Multiply it by the sum of the joint source factors at every actual size in
the requested interval, using the same z and alpha. This is not transport
of individual run classes or an anchor-size source-density substitution.

All raw entries must be positive; zero-ratio cases are rejected. Compact
zero-radius Arb upper endpoints are used for comparisons and census storage;
huge-exponent operators are never expanded into fmpq denominators. Partial
runs fail closed. No FLOAT helper is imported as a proof dependency.
No files are written. Root owns all runs; use Python -B and external limits.
"""

import argparse
from collections import Counter
from fractions import Fraction
import hashlib
from itertools import combinations_with_replacement
import json
import math
from pathlib import Path
import shutil
import sys
import time

sys.dont_write_bytecode = True
import verify_small_state_dense_anchor as exact


reference = exact.reference
GROUPS, PAIRS = exact.GROUPS, exact.PAIRS
OUTPUT_LIMIT = 256 * 1024
MAX_CORRECTION_ENTRIES = 2_000_000


def cells_for_group(groups, maximum, tail_counts):
    """Disjoint complete integer cells, including singleton endpoint cells."""
    low, high = max(0, groups - PAIRS), groups // 2
    tail = min(tail_counts, high - low + 1)
    bulk_high, cells = high - tail, []
    if bulk_high >= low:
        width = (bulk_high - low + maximum) // maximum
        cells.extend((start, min(bulk_high, start + width - 1))
                     for start in range(low, bulk_high + 1, width))
    cells.extend((full, full) for full in range(bulk_high + 1, high + 1))
    if (not 1 <= len(cells) <= maximum + tail_counts
            or cells[0][0] != low or cells[-1][1] != high
            or any(left[1] + 1 != right[0] for left, right in zip(cells, cells[1:]))):
        raise ArithmeticError("Incomplete or overlapping matching-count partition")
    return cells


def ratio_bin(left, right, subdivisions):
    """Exact floor(subdivisions*log2(left/right)); no floating assignment."""
    if left <= 0 or right <= 0:
        raise ArithmeticError("A raw ratio requires strictly positive entries")
    numerator = (left.numerator * right.denominator)**subdivisions
    denominator = (left.denominator * right.numerator)**subdivisions
    exponent = numerator.bit_length() - denominator.bit_length()
    below = (numerator < (denominator << exponent) if exponent >= 0
             else (numerator << -exponent) < denominator)
    return exponent - int(below)


def prepare_kernel(z, alpha, bin_width, local_moments, matrix, arb, arb_mat,
                   rational, guard):
    """Precompute local powers, exact bin membership, and correction matrices."""
    guard()
    profiles = tuple(combinations_with_replacement(range(17), 4))
    if len(profiles) != 4845:
        raise ArithmeticError("Incomplete sorted occupancy domain")
    divisors = [math.prod(math.factorial(n) for n in Counter(row).values())
                for row in profiles]
    subdivisions, exponent = bin_width.denominator, rational(alpha)
    local, epsilon = local_moments(z, q=256, h=4), Fraction(1, 256**16)
    powered_rows, memberships = [], [{}, {}, {}, {}]
    # Families0,1 are left row ratios; families2,3 are right column ratios.
    columns = ((0, 1), (2, 3), (0, 2), (1, 3))
    for index, occupancy in enumerate(profiles):
        if index % 32 == 0:
            guard()
        quiet, birth, ret, live = matrix(occupancy, z, local, q=256, h=4)
        total = quiet + birth
        entries = (quiet, total, epsilon * quiet + (1 - epsilon) * ret,
                   epsilon * total + (1 - epsilon) * live)
        if min(entries) <= 0 or max(entries) > 1:
            raise ArithmeticError("A raw rational local entry left(0,1]")
        if index == 0 and entries[:3] != (Fraction(1), Fraction(1), epsilon):
            raise ArithmeticError("Empty input lost its raw transitions")
        powered = [rational(value) if alpha == 1 else rational(value)**exponent
                   for value in entries]
        if any(not value.is_finite() or not value > 0 for value in powered):
            raise ArithmeticError("Invalid positive directed local moment")
        powered_rows.append(powered)
        for family, (first, second) in enumerate(columns):
            key = ratio_bin(entries[first], entries[second], subdivisions)
            memberships[family].setdefault(key, []).append(index)

    families = []
    for family, membership in enumerate(memberships):
        guard()
        keys = sorted(membership)
        indices = [tuple(membership[key]) for key in keys]
        if sorted(index for group in indices for index in group) != list(range(len(profiles))):
            raise ArithmeticError("Ratio bins did not partition all occupancies")
        first, second = columns[family]
        blocks = []
        for group in indices:
            guard()
            blocks.append(arb_mat([[powered_rows[index][first], powered_rows[index][second]]
                                   for index in group]))
        families.append(dict(keys=keys, indices=indices, blocks=blocks))

    total_entries = sum(len(families[a]["keys"]) * len(families[2 + b]["keys"])
                        for a in range(2) for b in range(2))
    if total_entries > MAX_CORRECTION_ENTRIES:
        raise MemoryError("ratio-bin-pair-cap")
    correction_by_sum, correction_by_distance = {}, {}

    def correction(total_bin):
        if total_bin not in correction_by_sum:
            # The closed log-ratio interval is[total_bin/q,(total_bin+2)/q].
            distance = max(abs(total_bin), abs(total_bin + 2))
            if distance not in correction_by_distance:
                guard()
                if alpha == 1:
                    upper = arb(1)
                else:
                    small = arb(2)**rational(Fraction(-distance, subdivisions))
                    value = (1 + small)**exponent / (1 + small**exponent)
                    if not value.is_finite() or not value > 0:
                        raise ArithmeticError("Invalid directed ratio correction")
                    upper = value.upper()
                    # c(r)<=1 analytically, so clipping an upper endpoint is safe.
                    if upper > 1:
                        upper = arb(1)
                if not upper > 0 or not upper <= 1:
                    raise ArithmeticError("Correction upper endpoint left(0,1]")
                correction_by_distance[distance] = upper
            correction_by_sum[total_bin] = correction_by_distance[distance]
        return correction_by_sum[total_bin]

    corrections = [[None, None], [None, None]]
    for a in range(2):
        for b in range(2):
            rows = []
            for left in families[a]["keys"]:
                guard()
                rows.append([correction(left + right) for right in families[2 + b]["keys"]])
            corrections[a][b] = arb_mat(rows)
            guard()
    receipt = dict(occupancies=len(profiles), raw_positive_entries=4 * len(profiles),
                   exact_rational_bin_assignment=True, bin_width=str(bin_width),
                   left_bin_counts=[len(families[a]["keys"]) for a in range(2)],
                   right_bin_counts=[len(families[2 + b]["keys"]) for b in range(2)],
                   correction_entries=total_entries, correction_entry_cap=MAX_CORRECTION_ENTRIES,
                   distinct_bin_sums=len(correction_by_sum),
                   distinct_correction_distances=len(correction_by_distance),
                   correction_rule="compact Arb upper at farther log-ratio endpoint, clipped at1")
    return dict(profiles=profiles, divisors=divisors, families=families,
                corrections=corrections), receipt


def regional_cell(prepared, groups, cell, tilt, pseudocount, arb, arb_mat,
                  rational, log_two, guard, details):
    """Bound the unconditioned grouped expectation, then condition once."""
    guard()
    probabilities = exact.proposal(groups, cell, tilt, pseudocount)
    weights, mass = exact.integer_occupancy_weights(
        probabilities, prepared["profiles"], prepared["divisors"], guard)
    aggregates = []
    for family_index, family in enumerate(prepared["families"]):
        values = []
        for indices, block in zip(family["indices"], family["blocks"]):
            guard()
            weighted = arb_mat([[weights[index] for index in indices]]) * block
            pair = (weighted[0, 0], weighted[0, 1])
            if any(not value.is_finite() or not value > 0 for value in pair):
                raise ArithmeticError("Invalid positive ratio-bin aggregate")
            values.append(pair)
        # Left matrices are2-by-bins; right matrices arebins-by2.
        aggregates.append(arb_mat([[row[0] for row in values], [row[1] for row in values]])
                          if family_index < 2 else arb_mat([[row[0], row[1]] for row in values]))

    old_average = arb_mat(2, 2)
    for a in range(2):
        for state in range(2):
            total = arb(0)
            for index in range(len(prepared["families"][a]["keys"])):
                total += aggregates[a][state, index]
            old_average[a, state] = total / arb(mass)
    old_window = old_average * old_average
    window = arb_mat(2, 2)
    intersections = 0
    for a in range(2):
        for b in range(2):
            guard()
            # Both paths are positive. No gain is subtracted from another bound.
            products = aggregates[a] * (prepared["corrections"][a][b] * aggregates[2 + b])
            value = (products[0, 0] + products[1, 1]) / arb(mass * mass)
            if not value.is_finite() or not value > 0:
                raise ArithmeticError("Invalid positive grouped IID entry")
            upper, baseline = value.upper(), old_window[a, b].upper()
            if not baseline.is_finite() or not baseline > 0:
                raise ArithmeticError("Invalid positive batchwise upper endpoint")
            if baseline < upper:
                upper = baseline
                intersections += 1
            window[a, b] = upper
    region = window
    for _ in range(6):
        guard()
        region = region * region

    lowers, endpoint_logs = [], []
    for full in cell:
        guard()
        histogram = exact.pair_histogram(groups, full)
        count = math.comb(PAIRS, histogram[0]) * math.comb(PAIRS - histogram[0], histogram[1])
        probability = arb(count)
        for p, n in zip(probabilities, histogram):
            probability *= rational(p)**n
        if not probability.is_finite() or not probability > 0:
            raise ArithmeticError("Nonpositive histogram probability enclosure")
        lower = probability.lower()
        if not lower > 0:
            raise ArithmeticError("Nonpositive histogram lower endpoint")
        lowers.append(lower)
        if details:
            endpoint_logs.append(str(probability.log() / log_two))
    divisor = min(lowers)
    region /= divisor
    upper_region = arb_mat(2, 2)
    for a in range(2):
        for b in range(2):
            value = region[a, b]
            if not value.is_finite() or not value > 0:
                raise ArithmeticError("Invalid normalized regional matrix")
            upper_region[a, b] = value.upper()
    receipt = dict(full_pair_low=cell[0], full_pair_high=cell[1],
                   covered_counts=cell[1] - cell[0] + 1,
                   exact_occupancy_mass_checked=True, proposal=list(map(str, probabilities)),
                   baseline_upper_intersections=intersections)
    if details:
        receipt.update(endpoint_probability_log2_intervals=endpoint_logs,
                       grouped_iid_log2_intervals=[[str(window[a, b].log() / log_two)
                                                   for b in range(2)] for a in range(2)],
                       regional_log2_intervals=[[str(upper_region[a, b].log() / log_two)
                                                for b in range(2)] for a in range(2)])
    guard()
    return upper_region, receipt


def cell_census(groups, cells, regions, arb, guard, progress):
    """Exact maxima of compact upper endpoints; exact integer run counts."""
    maximum_runs, count = min(groups, GROUPS - groups), len(cells)
    candidates = []
    expected_pairs = sum(max(1, groups - cells[i][1] - cells[j][1])
                         <= min(maximum_runs, groups - cells[i][0] - cells[j][0])
                         for i in range(count) for j in range(count))
    progress.update(requested_compatible_cell_pairs=expected_pairs,
                    computed_ordered_cell_pairs=0, assigned_run_classes=0,
                    counted_run_classes=0, complete_run_cover=False)
    for even, left in enumerate(regions):
        for odd, right in enumerate(regions):
            first = max(1, groups - cells[even][1] - cells[odd][1])
            last = min(maximum_runs, groups - cells[even][0] - cells[odd][0])
            if first > last:
                continue
            if len(candidates) % 32 == 0:
                guard()
            run = left * right
            for _ in range(7):
                run = run * run
            value = run[0, 0] + run[0, 1]
            if not value.is_finite() or not value > 0:
                raise ArithmeticError("Invalid complete ordered cell-pair operator")
            candidates.append((value.upper(), even, odd, first, last))
            progress["computed_ordered_cell_pairs"] = len(candidates)
    if len(candidates) != expected_pairs:
        raise ArithmeticError("Incomplete compatible cell-pair coverage")
    guard()
    # upper() endpoints are compact exact dyadics, so comparisons need no fmpq.
    candidates.sort(key=lambda row: row[0], reverse=True)
    successor = list(range(maximum_runs + 2))
    maxima, maximizing = [None] * (maximum_runs + 1), [None] * (maximum_runs + 1)

    def next_unfilled(runs):
        while successor[runs] != runs:
            successor[runs] = successor[successor[runs]]
            runs = successor[runs]
        return runs

    assigned = 0
    for index, (upper, even, odd, first, last) in enumerate(candidates):
        if index % 64 == 0:
            guard()
        runs = next_unfilled(first)
        while runs <= last:
            maxima[runs], maximizing[runs] = upper, (even, odd)
            assigned += 1
            successor[runs] = next_unfilled(runs + 1)
            runs = successor[runs]
        progress["assigned_run_classes"] = assigned
        if assigned == maximum_runs:
            break
    if assigned != maximum_runs or any(value is None for value in maxima[1:]):
        raise ArithmeticError("Incomplete source-run coverage")
    multiplicity, count_sum, total, terms = GROUPS, 0, arb(0), []
    for runs in range(1, maximum_runs + 1):
        if runs % 64 == 1:
            guard()
        if multiplicity <= 0:
            raise ArithmeticError("Nonpositive exact cycle-run count")
        term = maxima[runs] * multiplicity
        total += term
        terms.append((term.upper(), runs, multiplicity, maximizing[runs]))
        count_sum += multiplicity
        progress["counted_run_classes"] = runs
        if runs < maximum_runs:
            numerator = multiplicity * (groups - runs) * (GROUPS - groups - runs)
            multiplicity, remainder = divmod(numerator, runs * (runs + 1))
            if remainder:
                raise ArithmeticError("Nonintegral exact cycle-run recurrence")
    expected = math.comb(GROUPS, groups)
    if count_sum != expected:
        raise ArithmeticError("Cycle-run counts did not partition all source sets")
    if not total.is_finite() or not total > 0:
        raise ArithmeticError("Invalid positive route census")
    progress.update(exact_source_count_partition=True, complete_run_cover=True,
                    source_set_count_bits=expected.bit_length())
    terms.sort(key=lambda row: row[0], reverse=True)
    guard()
    return total, total / arb(expected), terms[:8]


def source_interval(groups, end_groups, z, alpha, arb, rational, guard, progress):
    """Joint source factors at the actual sizes; no source count is powered."""
    if (reference.LENGTH // 64 != 128 * 256
            or not 2 * (16 * (reference.LENGTH // 64)) < (1 << 32) - 1):
        raise ArithmeticError("The original full-refresh Gamma<2 receipt is unavailable")
    exponent = rational(alpha)
    beta = arb(256)**160 / arb(255)**32
    kappa = arb(256)**256 / beta
    growth = beta * kappa**rational(Fraction(1, 256))
    powered_growth = growth if alpha == 1 else growth**exponent
    common = 2 * rational(z)**(-reference.WEIGHT_LIMIT) / 255
    common = common if alpha == 1 else common**exponent
    anchor_density = common * powered_growth**groups
    term = arb(math.comb(GROUPS, groups)) * powered_growth**groups
    total = arb(0)
    for actual in range(groups, end_groups + 1):
        guard()
        total += term
        progress["completed_source_sizes"] += 1
        if actual < end_groups:
            term *= powered_growth * (GROUPS - actual)
            term /= actual + 1
    total *= common
    if not total.is_finite() or not total > 0:
        raise ArithmeticError("Invalid actual-size source-factor sum")
    return total, anchor_density


def verify(groups, z, alpha, end_groups=None, max_cells=128, tail_counts=32,
           tilt=Fraction(4, 3), pseudocount=Fraction(1, 65536),
           bin_width=Fraction(1, 4), precision=128, target_bits=100,
           seconds=50., memory_mib=1024, details=False):
    end_groups = groups if end_groups is None else end_groups
    z, alpha, tilt, pseudocount, bin_width = map(Fraction, (z, alpha, tilt, pseudocount, bin_width))
    if (not isinstance(groups, int) or not isinstance(end_groups, int)
            or not 1 <= groups < GROUPS or not groups <= end_groups <= GROUPS
            or not isinstance(max_cells, int) or not 1 <= max_cells <= 128
            or not isinstance(tail_counts, int) or not 0 <= tail_counts <= 32
            or not 0 < z < 1 or not 0 < alpha <= 1 or tilt <= 0 or pseudocount <= 0
            or bin_width not in (Fraction(1), Fraction(1, 4))
            or precision not in (128, 256) or not 1 <= target_bits <= 10000
            or not 0 < seconds <= 50 or not 64 <= memory_mib <= 1024):
        raise ValueError("Require proper anchor1..8191, anchor<=end<=8192,1..128 bulk cells,"
                         "0..32 tail cells, valid rational markers, width1 or1/4,"
                         "precision128/256, bounded target, <=50seconds and64..1024MiB")
    free = shutil.disk_usage(Path(__file__).resolve().parent).free
    if free < 20 * 1024**3:
        raise RuntimeError("At least20 GiB free space is required")
    started, next_memory = time.monotonic(), 0.

    def guard():
        nonlocal next_memory
        now = time.monotonic()
        if now - started >= seconds:
            raise TimeoutError("time-cap")
        if now >= next_memory:
            if reference.resident_bytes() > memory_mib * 1024**2:
                raise MemoryError("memory-cap")
            next_memory = now + .05

    cells = cells_for_group(groups, max_cells, tail_counts)
    regions, reports, dependencies, result = [], [], {}, {}
    progress = dict(computed_ordered_cell_pairs=0, assigned_run_classes=0,
                    counted_run_classes=0, complete_run_cover=False, completed_source_sizes=0)
    status, error, complete, passed = "directed-complete-grouped-route-replay", None, False, False
    try:
        guard()
        directory = (reference.OLD_THEORY if reference.OLD_THEORY.is_dir()
                     else Path(__file__).resolve().parent)
        sys.path.insert(0, str(directory))
        from hd_common_refresh_packets import local_moments, matrix
        from flint import arb, arb_mat, ctx
        ctx.prec = precision

        def rational(value):
            value = Fraction(value)
            return arb(value.numerator) / value.denominator

        log_two = arb(2).log()
        prepared, receipt = prepare_kernel(z, alpha, bin_width, local_moments, matrix,
                                           arb, arb_mat, rational, guard)
        result["local_kernel"] = receipt
        for index, cell in enumerate(cells):
            region, receipt = regional_cell(prepared, groups, cell, tilt, pseudocount,
                                              arb, arb_mat, rational, log_two, guard, details)
            regions.append(region)
            reports.append(dict(index=index, **receipt))
        if len(regions) != len(cells):
            raise ArithmeticError("Incomplete regional cell cover")
        census, average, top_terms = cell_census(groups, cells, regions, arb, guard, progress)
        source, anchor_density = source_interval(groups, end_groups, z, alpha,
                                                  arb, rational, guard, progress)
        bound = average * source
        if not bound.is_finite() or not bound > 0:
            raise ArithmeticError("Invalid final source-interval bound")
        top = []
        for upper, runs, count, (even, odd) in top_terms:
            top.append(dict(runs=runs, source_sets=str(count),
                            counted_route_log2_upper_interval=str(upper.log() / log_two),
                            anchor_log2_upper_interval=str((upper * anchor_density).log() / log_two),
                            even_cell_index=even, odd_cell_index=odd,
                            even_full_pair_cell=list(cells[even]), odd_full_pair_cell=list(cells[odd])))
        result.update(log2_bound_interval=str(bound.log() / log_two),
                      anchor_log2_bound_interval=str((census * anchor_density).log() / log_two),
                      raw_anchor_census_operator_sum_log2_interval=str(census.log() / log_two),
                      raw_anchor_average_operator_log2_interval=str(average.log() / log_two),
                      interval_source_factor_log2_interval=str(source.log() / log_two),
                      top_run_terms=top, worst_run_term=top[0])
        complete = (len(regions) == len(cells) and progress["complete_run_cover"]
                    and progress["completed_source_sizes"] == end_groups - groups + 1)
        passed = complete and bool(bound < arb(2)**(-target_bits))
        for name in reference.DEPENDENCIES:
            module = sys.modules.get(name)
            if module is None:
                raise RuntimeError("Missing loaded proof dependency: " + name)
            dependencies[name + ".py"] = hashlib.sha256(Path(module.__file__).read_bytes()).hexdigest()
        guard()
    except (TimeoutError, MemoryError, ArithmeticError, RuntimeError) as failure:
        status = str(failure) if isinstance(failure, (TimeoutError, MemoryError)) else "verification-error"
        error, complete, passed = str(failure), False, False
    for path in (Path(__file__), Path(exact.__file__), Path(reference.__file__)):
        dependencies[path.name] = hashlib.sha256(path.read_bytes()).hexdigest()
    result.update(schema="small-state-grouped-route-arb-v1", status=status, error=error,
                  passed=passed, complete_profile_cell_cover=complete,
                  complete_source_interval_cover=complete,
                  anchor_groups=groups, end_groups=end_groups, z=str(z), alpha=str(alpha),
                  tilt=str(tilt), pseudocount=str(pseudocount), bin_width=str(bin_width),
                  maximum_bulk_cells=max_cells, requested_tail_counts=tail_counts,
                  requested_cells=len(cells), completed_cells=len(regions),
                  unprocessed_cells=len(cells) - len(regions), cell_partition=[list(cell) for cell in cells],
                  requested_run_classes=min(groups, GROUPS - groups),
                  unprocessed_run_classes=min(groups, GROUPS - groups) - progress["counted_run_classes"],
                  requested_source_sizes=end_groups - groups + 1,
                  completed_source_sizes=progress["completed_source_sizes"],
                  unprocessed_source_sizes=end_groups - groups + 1 - progress["completed_source_sizes"],
                  census_engine=progress, precision_bits=precision, target_bits=target_bits,
                  acceptance="complete coverage and Arb bound < exact2^(-target_bits)",
                  source_mode="joint", gamma_upper="2",
                  grouping="two ordered raw batches before entrywise alpha",
                  regional_rule="unconditioned grouped IID W^64 / full4096-pair histogram probability",
                  histogram_cell_rule="minimum compact exact lower endpoint; divisor not powered by alpha",
                  correction_rule="exact rational ratio bins; positive Arb correction aggregates",
                  baseline_intersection="entrywise minimum of grouped and batchwise upper endpoints",
                  profile_accounting="exact run counts times maximum compatible ordered cell-pair bound",
                  anchor_accounting="route-only census divided by binom(8192,anchor)",
                  source_accounting="actual-r joint factors and source binomials; counts outside alpha",
                  endpoint_storage="compact exact dyadic Arb upper endpoints; no operator fmpq conversion",
                  interval_rule="exact grouped raw inclusion monotonicity; common z and alpha",
                  empty_input_rule="raw[[1,1],[256^-16,V_empty]]; omega=1; positive entries required",
                  terminal="positive f<=M*ones; original Gamma; no final refresh or reset charge",
                  seconds=time.monotonic() - started, peak_resident_bytes=reference.resident_bytes(),
                  time_cap_seconds=seconds, memory_cap_bytes=memory_mib * 1024**2,
                  free_bytes_at_start=free, output_cap_bytes=OUTPUT_LIMIT, dependencies=dependencies,
                  directed_arithmetic=True, directed=complete,
                  whole_code_certificate=False, generated_file_bytes=0)
    if details:
        result["cells"] = reports
    return result


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--groups", type=int, default=2048)
    parser.add_argument("--z", type=Fraction, required=True)
    parser.add_argument("--alpha", type=Fraction, required=True)
    parser.add_argument("--end-groups", type=int)
    parser.add_argument("--max-cells", type=int, default=128)
    parser.add_argument("--tail-counts", type=int, default=32)
    parser.add_argument("--tilt", type=Fraction, default=Fraction(4, 3))
    parser.add_argument("--pseudocount", type=Fraction, default=Fraction(1, 65536))
    parser.add_argument("--bin-width", type=Fraction, default=Fraction(1, 4))
    parser.add_argument("--precision", type=int, choices=(128, 256), default=128)
    parser.add_argument("--target-bits", type=int, default=100)
    parser.add_argument("--seconds", type=float, default=50.)
    parser.add_argument("--memory-mib", type=int, default=1024)
    parser.add_argument("--details", action="store_true")
    args = parser.parse_args()
    result = verify(**vars(args))
    encoded = json.dumps(result, allow_nan=False, sort_keys=True)
    if len(encoded.encode("utf-8")) > OUTPUT_LIMIT:
        result.pop("cells", None)
        result["details_omitted_for_output_cap"] = True
        encoded = json.dumps(result, allow_nan=False, sort_keys=True)
    if len(encoded.encode("utf-8")) > OUTPUT_LIMIT:
        raise RuntimeError("Compact grouped receipt exceeded its output cap")
    print(encoded, flush=True)
    if not result["passed"]:
        raise SystemExit(1)


if __name__ == "__main__":
    main()
