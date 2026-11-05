"""Directed actual-size interval with fixed parity IID proposal families.

For each parity's first size, cache positive unconditioned regional upper
matrices K_i=W(p_i)^64 and all ordered C_ij=e0*(K_i*K_j)^128*ones. Neither
depends on the later actual source size. Translate cells by one full pair
per two added groups, extending the lowest cell to the actual feasible
lower bound. Check the complete disjoint partition separately at every g.

At each actual g, bound both endpoint histogram probabilities from below.
Their minimum h_i(g) covers that entire cell. The ordered route bound is
C_ij/(h_i(g)*h_j(g))^128. Rebuild all compatible pair bounds and their order
at every g, then apply exact cycle-run multiplicities N_g(c). Check their
integer sum against binom(8192,g). Multiply by the actual-g joint source
factor and sum the requested sizes. No averaged-anchor transport, anchor
normalization, or extra source binomial occurs.

All proof arithmetic is positive Arb arithmetic at128 or256 bits. Compact
upper/lower endpoints avoid huge fmpq denominators. After endpoint sorting,
certify the complete descending order before using its interval-max sweep;
an ambiguous comparison fails closed. Partial coverage never passes.
No FLOAT proof helper is imported. No files are written. Root owns all
runs; use Python -B and external50-second/1-GiB limits. Internal caps also
require20 GiB free disk and at most256 actual sizes.
"""

import argparse
from fractions import Fraction
import hashlib
import json
import math
from pathlib import Path
import shutil
import sys
import time

sys.dont_write_bytecode = True
import verify_small_state_grouped_route as grouped


exact, reference = grouped.exact, grouped.reference
GROUPS, PAIRS = grouped.GROUPS, grouped.PAIRS
OUTPUT_LIMIT, MAX_SOURCE_SIZES = 256 * 1024, 256


def moving_cells(groups, anchor, anchor_cells):
    """Translate fixed-deficit cells, extend the lowest, and check coverage."""
    if groups < anchor or (groups-anchor) % 2:
        raise ValueError("A fixed proposal family requires its anchor parity")
    shift = (groups-anchor)//2
    feasible_low, feasible_high = max(0, groups-PAIRS), groups//2
    active, cells = [], []
    for index, (low, high) in enumerate(anchor_cells):
        low, high = low+shift, high+shift
        if index == 0:
            low = feasible_low
        low, high = max(feasible_low, low), min(feasible_high, high)
        if low <= high:
            active.append(index)
            cells.append((low, high))
    if (not cells or cells[0][0] != feasible_low or cells[-1][1] != feasible_high
            or any(left[1]+1 != right[0] for left, right in zip(cells, cells[1:]))
            or sum(high-low+1 for low, high in cells) != feasible_high-feasible_low+1):
        raise ArithmeticError("Incomplete translated actual-size partition")
    return active, cells


def positive_upper(value, context):
    if not value.is_finite() or not value > 0:
        raise ArithmeticError("Invalid positive enclosure: " + context)
    return value.upper()


def regional_numerator(prepared, probabilities, arb, arb_mat, guard):
    """The existing positive grouped IID contraction, without conditioning."""
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
                raise ArithmeticError("Invalid positive grouped occupancy aggregate")
            values.append(pair)
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
    window, intersections = arb_mat(2, 2), 0
    for a in range(2):
        for b in range(2):
            guard()
            products = aggregates[a] * (prepared["corrections"][a][b] * aggregates[2+b])
            value = (products[0, 0]+products[1, 1]) / arb(mass*mass)
            upper = positive_upper(value, "grouped IID entry")
            baseline = positive_upper(old_window[a, b], "batchwise IID entry")
            if baseline <= upper:
                upper = baseline
                intersections += 1
            elif not upper <= baseline:
                raise ArithmeticError("Ambiguous compact IID endpoint comparison")
            window[a, b] = upper
    region = window
    for _ in range(6):
        guard()
        region = region * region
    result = arb_mat(2, 2)
    for a in range(2):
        for b in range(2):
            result[a, b] = positive_upper(region[a, b], "unconditioned regional numerator")
    guard()
    return result, intersections


def prepare_family(prepared, anchor, cells, tilt, pseudocount, arb, arb_mat,
                   log_two, guard, progress, details):
    """Cache all ordered numerator operators, not only anchor-feasible ones."""
    numerators, proposals, reports = [], [], []
    for index, cell in enumerate(cells):
        guard()
        probabilities = exact.proposal(anchor, cell, tilt, pseudocount)
        numerator, intersections = regional_numerator(prepared, probabilities, arb, arb_mat, guard)
        proposals.append(probabilities)
        numerators.append(numerator)
        if details:
            reports.append(dict(index=index, anchor_cell=list(cell), proposal=list(map(str, probabilities)),
                                batchwise_intersections=intersections,
                                numerator_log2_intervals=[[str(numerator[a, b].log()/log_two)
                                                          for b in range(2)] for a in range(2)]))
        progress["completed_family_cells"] += 1
    count = len(cells)
    if len(numerators) != count:
        raise ArithmeticError("Incomplete parity-family numerator cover")
    operators = []
    for even, left in enumerate(numerators):
        row = []
        for odd, right in enumerate(numerators):
            if (odd & 31) == 0:
                guard()
            run = left * right
            for _ in range(7):
                run = run * run
            row.append(positive_upper(run[0, 0]+run[0, 1], "ordered numerator operator"))
            progress["cached_ordered_cell_pairs"] += 1
        operators.append(row)
    if len(operators) != count or any(len(row) != count for row in operators):
        raise ArithmeticError("Incomplete ordered numerator cache")
    guard()
    return dict(anchor=anchor, cells=cells, proposals=proposals, numerators=numerators,
                operators=operators, reports=reports)


def histogram_lower(groups, full, probabilities, arb, rational):
    histogram = exact.pair_histogram(groups, full)
    count = math.comb(PAIRS, histogram[0]) * math.comb(PAIRS-histogram[0], histogram[1])
    probability = arb(count)
    for p, n in zip(probabilities, histogram):
        probability *= rational(p)**n
    if not probability.is_finite() or not probability > 0:
        raise ArithmeticError("Invalid actual-size histogram enclosure")
    lower = probability.lower()
    if not lower > 0:
        raise ArithmeticError("Nonpositive compact histogram lower endpoint")
    return lower


def actual_census(groups, family, arb, rational, log_two, guard, progress):
    """Recondition every cell and rebuild the exact actual-g run census."""
    guard()
    active, cells = moving_cells(groups, family["anchor"], family["cells"])
    maximum_runs, count = min(groups, GROUPS-groups), len(cells)
    progress.update(current_groups=groups, current_requested_cells=count,
                    current_completed_cells=0, current_computed_cell_pairs=0,
                    current_assigned_run_classes=0, current_counted_run_classes=0,
                    current_complete_run_cover=False)
    scales = []
    for family_index, cell in zip(active, cells):
        guard()
        probabilities = family["proposals"][family_index]
        left = histogram_lower(groups, cell[0], probabilities, arb, rational)
        right = histogram_lower(groups, cell[1], probabilities, arb, rational)
        if left <= right:
            divisor = left
        elif right <= left:
            divisor = right
        else:
            raise ArithmeticError("Ambiguous compact histogram endpoint comparison")
        # Each regional normalizer occurs128 times. Precompute its positive
        # scalar power once per cell, rather than once per ordered pair.
        scales.append(positive_upper(divisor**(-128), "inverse histogram power"))
        progress["current_completed_cells"] += 1
    candidates, expected_pairs = [], 0
    for even, (left_low, left_high) in enumerate(cells):
        for odd, (right_low, right_high) in enumerate(cells):
            first = max(1, groups-left_high-right_high)
            last = min(maximum_runs, groups-left_low-right_low)
            if first > last:
                continue
            expected_pairs += 1
            if (expected_pairs & 63) == 1:
                guard()
            numerator = family["operators"][active[even]][active[odd]]
            upper = positive_upper(numerator*scales[even]*scales[odd], "conditioned cell-pair bound")
            candidates.append((upper, even, odd, first, last))
            progress["current_computed_cell_pairs"] = len(candidates)
    if not candidates or len(candidates) != expected_pairs:
        raise ArithmeticError("Incomplete actual-size compatible-pair cover")
    guard()
    candidates.sort(key=lambda row: row[0], reverse=True)
    # upper() returns compact exact dyadics. Explicitly certify the whole
    # order; overlapping or otherwise ambiguous comparisons cannot pass.
    for index in range(1, len(candidates)):
        if (index & 255) == 0:
            guard()
        if not candidates[index-1][0] >= candidates[index][0]:
            raise ArithmeticError("Uncertified descending upper-endpoint order")
    successor = list(range(maximum_runs+2))
    maxima, maximizing = [None]*(maximum_runs+1), [None]*(maximum_runs+1)

    def next_unfilled(runs):
        while successor[runs] != runs:
            successor[runs] = successor[successor[runs]]
            runs = successor[runs]
        return runs

    assigned = 0
    for index, (upper, even, odd, first, last) in enumerate(candidates):
        if (index & 63) == 0:
            guard()
        runs = next_unfilled(first)
        while runs <= last:
            if (assigned & 63) == 0:
                guard()
            maxima[runs], maximizing[runs] = upper, (even, odd)
            assigned += 1
            successor[runs] = next_unfilled(runs+1)
            runs = successor[runs]
        progress["current_assigned_run_classes"] = assigned
        if assigned == maximum_runs:
            break
    if assigned != maximum_runs or any(value is None for value in maxima[1:]):
        raise ArithmeticError("Incomplete actual-size run-class cover")
    multiplicity, count_sum, total, worst = GROUPS, 0, arb(0), None
    for runs in range(1, maximum_runs+1):
        if (runs & 63) == 1:
            guard()
        if multiplicity <= 0:
            raise ArithmeticError("Nonpositive exact cycle-run multiplicity")
        term = maxima[runs] * multiplicity
        upper = positive_upper(term, "counted run term")
        total += term
        if worst is None or upper > worst[0]:
            worst = (upper, runs, maximizing[runs])
        count_sum += multiplicity
        progress["current_counted_run_classes"] = runs
        if runs < maximum_runs:
            numerator = multiplicity*(groups-runs)*(GROUPS-groups-runs)
            multiplicity, remainder = divmod(numerator, runs*(runs+1))
            if remainder:
                raise ArithmeticError("Nonintegral exact cycle-run recurrence")
    expected = math.comb(GROUPS, groups)
    if count_sum != expected:
        raise ArithmeticError("Exact actual-size census did not partition all source sets")
    positive_upper(total, "complete actual-size route census")
    _, runs, (even, odd) = worst
    row = dict(groups=groups, family_anchor=family["anchor"], active_cells=count,
               covered_full_pair_counts=sum(high-low+1 for low, high in cells),
               compatible_ordered_cell_pairs=expected_pairs, run_classes=maximum_runs,
               exact_source_count_partition=True, source_set_count_bits=expected.bit_length(),
               route_census_log2_interval=str(total.log()/log_two),
               worst_run=runs, worst_even_cell=list(cells[even]), worst_odd_cell=list(cells[odd]))
    progress["current_complete_run_cover"] = True
    guard()
    return total, row


def source_growth(z, alpha, arb, rational):
    """Return the original joint prefactor common*growth^actual_g."""
    if (reference.LENGTH//64 != 128*256
            or not 2*(16*(reference.LENGTH//64)) < (1 << 32)-1):
        raise ArithmeticError("The original full-refresh Gamma<2 receipt is unavailable")
    exponent = rational(alpha)
    beta = arb(256)**160 / arb(255)**32
    kappa = arb(256)**256 / beta
    growth = beta * kappa**rational(Fraction(1, 256))
    common = 2 * rational(z)**(-reference.WEIGHT_LIMIT) / 255
    if alpha != 1:
        common, growth = common**exponent, growth**exponent
    positive_upper(common, "common joint source factor")
    positive_upper(growth, "actual-size source growth")
    return common, growth


def verify(groups=2048, end_groups=2073, z=Fraction(23, 100), alpha=Fraction(19, 200),
           max_cells=128, tail_counts=32, tilt=Fraction(33, 25),
           pseudocount=Fraction(1, 65536), bin_width=Fraction(1, 4),
           precision=128, target_bits=100, seconds=50., memory_mib=1024, details=False):
    end_groups = groups if end_groups is None else end_groups
    z, alpha, tilt, pseudocount, bin_width = map(Fraction, (z, alpha, tilt, pseudocount, bin_width))
    if (not isinstance(groups, int) or not isinstance(end_groups, int)
            or not 1 <= groups <= end_groups < GROUPS or end_groups-groups+1 > MAX_SOURCE_SIZES
            or not isinstance(max_cells, int) or not 1 <= max_cells <= 128
            or not isinstance(tail_counts, int) or not 0 <= tail_counts <= 32
            or not 0 < z < 1 or not 0 < alpha <= 1 or tilt <= 0 or pseudocount <= 0
            or bin_width not in (Fraction(1), Fraction(1, 4)) or precision not in (128, 256)
            or not isinstance(target_bits, int) or not 1 <= target_bits <= 10000
            or not 0 < seconds <= 50 or not 64 <= memory_mib <= 1024):
        raise ValueError("Require <=256 proper actual sizes, <=128 bulk cells, <=32 tail cells,"
                         "positive rational markers, width1 or1/4, precision128/256, <=50s/1GiB")
    free = shutil.disk_usage(Path(__file__).resolve().parent).free
    if free < 20*1024**3:
        raise RuntimeError("At least20 GiB free space is required")
    started, next_memory = time.monotonic(), 0.

    def guard():
        nonlocal next_memory
        now = time.monotonic()
        if now-started >= seconds:
            raise TimeoutError("time-cap")
        if now >= next_memory:
            if reference.resident_bytes() > memory_mib*1024**2:
                raise MemoryError("memory-cap")
            next_memory = now+.05

    anchors = [groups] if groups == end_groups else [groups, groups+1]
    partitions = {anchor: grouped.cells_for_group(anchor, max_cells, tail_counts) for anchor in anchors}
    progress = dict(completed_family_cells=0, cached_ordered_cell_pairs=0,
                    current_groups=None, current_requested_cells=0, current_completed_cells=0,
                    current_computed_cell_pairs=0, current_assigned_run_classes=0,
                    current_counted_run_classes=0, current_complete_run_cover=False)
    families, rows, result, dependencies = {}, [], {}, {}
    status, error, complete, passed = "directed-complete-grouped-interval-replay", None, False, False
    try:
        guard()
        directory = reference.OLD_THEORY if reference.OLD_THEORY.is_dir() else Path(__file__).resolve().parent
        sys.path.insert(0, str(directory))
        from hd_common_refresh_packets import local_moments, matrix
        from flint import arb, arb_mat, ctx
        ctx.prec = precision

        def rational(value):
            value = Fraction(value)
            return arb(value.numerator)/value.denominator

        log_two = arb(2).log()
        prepared, receipt = grouped.prepare_kernel(z, alpha, bin_width, local_moments, matrix,
                                                   arb, arb_mat, rational, guard)
        result["local_kernel"] = receipt
        for anchor in anchors:
            families[anchor % 2] = prepare_family(prepared, anchor, partitions[anchor], tilt,
                                                  pseudocount, arb, arb_mat, log_two,
                                                  guard, progress, details)
        common, growth = source_growth(z, alpha, arb, rational)
        density, union = common * growth**groups, arb(0)
        for actual in range(groups, end_groups+1):
            guard()
            census, row = actual_census(actual, families[actual % 2], arb, rational,
                                         log_two, guard, progress)
            bound = census * density
            positive_upper(bound, "actual-size union bound")
            union += bound
            row.update(source_factor_log2_interval=str(density.log()/log_two),
                       log2_bound_interval=str(bound.log()/log_two))
            rows.append(row)
            if actual < end_groups:
                density *= growth
        positive_upper(union, "complete source-interval union")
        if len(rows) != end_groups-groups+1:
            raise ArithmeticError("Incomplete direct source-size interval")
        for name in reference.DEPENDENCIES:
            guard()
            module = sys.modules.get(name)
            if module is None:
                raise RuntimeError("Missing loaded proof dependency: " + name)
            dependencies[name+".py"] = hashlib.sha256(Path(module.__file__).read_bytes()).hexdigest()
        guard()
        result["log2_bound_interval"] = str(union.log()/log_two)
        complete = True
        passed = bool(union < arb(2)**(-target_bits))
    except (TimeoutError, MemoryError, ArithmeticError, RuntimeError) as failure:
        status = str(failure) if isinstance(failure, (TimeoutError, MemoryError)) else "verification-error"
        error, complete, passed = str(failure), False, False
    for path in (Path(__file__), Path(grouped.__file__), Path(exact.__file__), Path(reference.__file__)):
        dependencies[path.name] = hashlib.sha256(path.read_bytes()).hexdigest()
    requested_cells = sum(len(cells) for cells in partitions.values())
    requested_pairs = sum(len(cells)**2 for cells in partitions.values())
    requested_runs = sum(min(g, GROUPS-g) for g in range(groups, end_groups+1))
    completed_runs = sum(row["run_classes"] for row in rows)
    result.update(schema="small-state-grouped-interval-arb-v1", status=status, error=error,
                  passed=passed, complete_profile_cell_cover=complete, complete_source_interval_cover=complete,
                  first_groups=groups, end_groups=end_groups, z=str(z), alpha=str(alpha),
                  tilt=str(tilt), pseudocount=str(pseudocount), bin_width=str(bin_width),
                  maximum_bulk_cells=max_cells, requested_tail_counts=tail_counts,
                  parity_anchors=anchors, completed_families=len(families),
                  family_cell_counts={str(anchor): len(cells) for anchor, cells in partitions.items()},
                  requested_family_cells=requested_cells,
                  unprocessed_family_cells=requested_cells-progress["completed_family_cells"],
                  requested_cached_ordered_pairs=requested_pairs,
                  unprocessed_cached_ordered_pairs=requested_pairs-progress["cached_ordered_cell_pairs"],
                  requested_source_sizes=end_groups-groups+1, completed_source_sizes=len(rows),
                  unprocessed_source_sizes=end_groups-groups+1-len(rows),
                  requested_run_classes=requested_runs, completed_run_classes=completed_runs,
                  unprocessed_run_classes=requested_runs-completed_runs,
                  progress=progress, source_sizes=rows, precision_bits=precision, target_bits=target_bits,
                  acceptance="complete coverage and Arb union < exact2^(-target_bits)",
                  interval_rule="direct actual-g conditioning and census; no subset transport",
                  cache_rule="fixed parity proposals; unconditioned K_i=W_i^64; every ordered C_ij cached",
                  cell_rule="shift high-f deficit; extend lowest to feasible lower; clip and check each size",
                  histogram_rule="minimum certified lower endpoint; each inverse divisor powered128, not alpha",
                  maxima_rule="compact upper endpoints; complete descending order certified before max sweep",
                  profile_accounting="exact N_g(c) and complete compatible ordered cell pairs at every actual g",
                  source_accounting="actual-g joint source factor; no extra binomial or anchor normalization",
                  grouping="two ordered raw batches before entrywise alpha", source_mode="joint", gamma_upper="2",
                  empty_input_rule="raw[[1,1],[256^-16,V_empty]]; omega=1; strictly positive entries",
                  terminal="positive f<=M*ones; original Gamma; no terminal reset or final refresh",
                  endpoint_storage="compact Arb endpoints; no operator fmpq conversion",
                  seconds=time.monotonic()-started, peak_resident_bytes=reference.resident_bytes(),
                  time_cap_seconds=seconds, memory_cap_bytes=memory_mib*1024**2,
                  free_bytes_at_start=free, output_cap_bytes=OUTPUT_LIMIT, dependencies=dependencies,
                  directed_arithmetic=True, directed=complete, whole_code_certificate=False, generated_file_bytes=0)
    if details:
        result["families"] = [dict(anchor=family["anchor"], cells=family["reports"])
                              for family in families.values()]
    return result


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--groups", type=int, default=2048)
    parser.add_argument("--end-groups", type=int, default=2073)
    parser.add_argument("--z", type=Fraction, default=Fraction(23, 100))
    parser.add_argument("--alpha", type=Fraction, default=Fraction(19, 200))
    parser.add_argument("--max-cells", type=int, default=128)
    parser.add_argument("--tail-counts", type=int, default=32)
    parser.add_argument("--tilt", type=Fraction, default=Fraction(33, 25))
    parser.add_argument("--pseudocount", type=Fraction, default=Fraction(1, 65536))
    parser.add_argument("--bin-width", type=Fraction, default=Fraction(1, 4))
    parser.add_argument("--precision", type=int, choices=(128, 256), default=128)
    parser.add_argument("--target-bits", type=int, default=100)
    parser.add_argument("--seconds", type=float, default=50.)
    parser.add_argument("--memory-mib", type=int, default=1024)
    parser.add_argument("--details", action="store_true")
    result = verify(**vars(parser.parse_args()))
    encoded = json.dumps(result, allow_nan=False, sort_keys=True)
    if len(encoded.encode("utf-8")) > OUTPUT_LIMIT:
        result.pop("families", None)
        result["details_omitted_for_output_cap"] = True
        encoded = json.dumps(result, allow_nan=False, sort_keys=True)
    if len(encoded.encode("utf-8")) > OUTPUT_LIMIT:
        raise RuntimeError("Compact directed-interval receipt exceeded its output cap")
    print(encoded, flush=True)
    if not result["passed"]:
        raise SystemExit(1)


if __name__ == "__main__":
    main()
