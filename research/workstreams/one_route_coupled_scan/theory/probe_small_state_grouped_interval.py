"""FLOAT direct source-size interval using fixed parity proposal families.

At each parity's first size, prepare the existing complete grouped cells and
their positive IID proposals. K_i=W(p_i)^64 is an UNCONDITIONED regional
numerator, independent of the actual size g. Cache every ordered
C_ij=[1,0](K_i*K_j)^128*ones, including pairs impossible at the anchor.

As g increases by two, translate every cell by one full pair. This holds
the high-f deficit fixed. Extend the lowest cell to the actual feasible
lower bound, and clip cells when necessary. Check complete integer coverage
at EVERY g. A fixed positive proposal remains valid on an expanded cell.
Separate parity families avoid using a zero-singleton proposal for one
singleton. No regional override or selected-profile shortcut is accepted.

For each actual g, h_i(g) is the smaller endpoint histogram probability.
The ordered cell-pair bound is C_ij/(h_i(g)*h_j(g))^128. Rebuild the run-line
maxima after these scalar corrections; their sorting can change with g.
Multiply by exact N_g(c), sum all runs, then use the actual-g joint source
factor. Sum these direct size bounds, with no averaged-anchor transport,
anchor binomial divisor, or additional source/profile count.

All arithmetic except bin assignment and integer census counts is FLOAT.
Finite log arithmetic avoids materializing tiny positive operators. This
is not a directed certificate. Caps: 256 sizes, 50 seconds, 1 GiB, 20 GiB
free disk, compact stdout only. No files are written. Root owns all runs.
"""

import argparse
from fractions import Fraction
import hashlib
import json
import math
from pathlib import Path
import sys

sys.dont_write_bytecode = True
import probe_small_state_grouped_route as grouped


census, dense, bridge, np = grouped.census, grouped.dense, grouped.bridge, grouped.np
OUTPUT_LIMIT = 256 * 1024
MAX_SOURCE_SIZES = 256


def moving_cells(groups, anchor, anchor_cells):
    """Return active family indices and a complete translated partition."""
    if groups < anchor or (groups-anchor) % 2:
        raise ValueError("A proposal family requires its anchor parity and later sizes")
    shift = (groups-anchor)//2
    feasible_low, feasible_high = max(0, groups-bridge.PAIRS), groups//2
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
        raise ArithmeticError("Incomplete translated actual-size cell partition")
    return active, cells


def prepare_family(prepared, anchor, anchor_cells, alpha, tilt, pseudocount,
                   guard, progress):
    """Cache log2 K_i and log2 C_ij once; expose both for root cross-checks."""
    numerators, proposals, reports = [], [], []
    for index, cell in enumerate(anchor_cells):
        guard()
        _, receipt = grouped.regional_cell(prepared, anchor, cell, alpha, tilt,
                                           pseudocount, guard)
        # Rebuild K from W, rather than undoing a large conditioning division.
        value = np.asarray(receipt["grouped_iid_log2"], dtype=float)
        with np.errstate(over="raise", invalid="raise", under="ignore"):
            for _ in range(6):
                guard()
                value = dense.log_matrix_product(value, value)
        probabilities = np.asarray([float(Fraction(p)) for p in receipt["proposal"]])
        if (value.shape != (2, 2) or np.any(~np.isfinite(value))
                or np.any(~np.isfinite(probabilities)) or np.any(probabilities <= 0)
                or abs(float(probabilities.sum())-1.) > 1e-12):
            raise ArithmeticError("Invalid fixed positive regional numerator or proposal")
        numerators.append(value)
        proposals.append(probabilities)
        reports.append(dict(index=index, anchor_cell=list(cell), proposal=receipt["proposal"],
                            numerator_log2=value.tolist()))
        progress["completed_family_cells"] += 1
    numerators = np.asarray(numerators)
    count = len(anchor_cells)
    if numerators.shape != (count, 2, 2):
        raise ArithmeticError("Incomplete parity-family regional numerators")
    with np.errstate(over="raise", invalid="raise", under="ignore"):
        run = dense.log_matrix_product(numerators[:, None], numerators[None, :])
        for _ in range(7):
            guard()
            run = dense.log_matrix_product(run, run)
        operators = np.logaddexp2(run[..., 0, 0], run[..., 0, 1])
    if operators.shape != (count, count) or np.any(~np.isfinite(operators)):
        raise ArithmeticError("Invalid cached ordered numerator operators")
    progress["cached_ordered_cell_pairs"] += count*count
    guard()
    return dict(anchor=anchor, cells=anchor_cells, proposals=np.asarray(proposals),
                numerator_log2=numerators, operator_log2=operators, reports=reports)


def actual_census(groups, family, guard, progress, details=False):
    """Apply actual-g scalar denominators, then count every actual source set."""
    guard()
    active, cells = moving_cells(groups, family["anchor"], family["cells"])
    count, maximum_runs = len(cells), min(groups, bridge.GROUPS-groups)
    divisors = np.empty(count)
    progress.update(current_groups=groups, current_completed_cells=0,
                    current_requested_cells=count, current_assigned_run_classes=0,
                    current_counted_run_classes=0, current_complete_run_cover=False)
    for index, (family_index, cell) in enumerate(zip(active, cells)):
        if (index & 15) == 0:
            guard()
        probabilities = family["proposals"][family_index]
        endpoints = [bridge.histogram_log_probability(bridge.pair_histogram(groups, full),
                                                       probabilities) for full in cell]
        divisors[index] = min(endpoints)
        if not all(math.isfinite(value) for value in endpoints):
            raise ArithmeticError("Nonfinite actual-size endpoint histogram probability")
        progress["current_completed_cells"] += 1
    # Scalar normalizers commute with matrices; the two regional divisors
    # each occur exactly128 times. No alpha acts on either divisor.
    with np.errstate(over="raise", invalid="raise", under="ignore"):
        operators = family["operator_log2"][np.ix_(active, active)].copy()
        operators -= 128.*(divisors[:, None]+divisors[None, :])
    if np.any(~np.isfinite(operators)):
        raise ArithmeticError("Nonfinite actual-size ordered cell-pair bound")
    lower = np.asarray([cell[0] for cell in cells], dtype=int)
    upper = np.asarray([cell[1] for cell in cells], dtype=int)
    first = np.maximum(1, groups-upper[:, None]-upper[None, :])
    last = np.minimum(maximum_runs, groups-lower[:, None]-lower[None, :])
    candidates = np.flatnonzero(first <= last)
    # Re-sort for every g: the scalar denominator can change pair ordering.
    order = candidates[np.argsort(operators.ravel()[candidates])[::-1]]
    successor = list(range(maximum_runs+2))
    maxima = np.full(maximum_runs+1, -math.inf)
    maximizing = np.full((maximum_runs+1, 2), -1, dtype=int)

    def next_unfilled(runs):
        while successor[runs] != runs:
            successor[runs] = successor[successor[runs]]
            runs = successor[runs]
        return runs

    assigned = 0
    for index, flat in enumerate(order):
        if (index & 63) == 0:
            guard()
        even, odd = divmod(int(flat), count)
        runs, end = next_unfilled(int(first[even, odd])), int(last[even, odd])
        while runs <= end:
            if (assigned & 63) == 0:
                guard()
            maxima[runs], maximizing[runs] = operators[even, odd], (even, odd)
            assigned += 1
            successor[runs] = next_unfilled(runs+1)
            runs = successor[runs]
        progress["current_assigned_run_classes"] = assigned
        if assigned == maximum_runs:
            break
    if assigned != maximum_runs or np.any(~np.isfinite(maxima[1:])):
        raise ArithmeticError("Incomplete actual-size source-run cover")

    multiplicity, count_sum = bridge.GROUPS, 0
    terms = np.empty(maximum_runs)
    for runs in range(1, maximum_runs+1):
        if (runs & 63) == 1:
            guard()
        if multiplicity <= 0:
            raise ArithmeticError("Nonpositive exact cycle-run multiplicity")
        terms[runs-1] = math.log2(multiplicity)+maxima[runs]
        count_sum += multiplicity
        progress["current_counted_run_classes"] = runs
        if runs < maximum_runs:
            numerator = multiplicity*(groups-runs)*(bridge.GROUPS-groups-runs)
            multiplicity, remainder = divmod(numerator, runs*(runs+1))
            if remainder:
                raise ArithmeticError("Nonintegral exact cycle-run recurrence")
    expected = math.comb(bridge.GROUPS, groups)
    if count_sum != expected:
        raise ArithmeticError("Actual-size census did not partition all source sets")
    with np.errstate(over="raise", invalid="raise", under="ignore"):
        total = float(np.logaddexp2.reduce(terms))
    if not math.isfinite(total) or np.any(~np.isfinite(terms)):
        raise ArithmeticError("Nonfinite complete actual-size route census")
    worst = int(np.argmax(terms))+1
    even, odd = map(int, maximizing[worst])
    row = dict(groups=groups, family_anchor=family["anchor"], active_cells=count,
               covered_full_pair_counts=sum(high-low+1 for low, high in cells),
               compatible_ordered_cell_pairs=len(candidates), run_classes=maximum_runs,
               exact_source_count_partition=True, source_set_count_bits=expected.bit_length(),
               route_census_log2=total,
               worst_run_term=dict(runs=worst, counted_route_log2=float(terms[worst-1]),
                                   even_family_index=active[even], odd_family_index=active[odd],
                                   even_cell=list(cells[even]), odd_cell=list(cells[odd])))
    if details:
        row["minimum_histogram_log2"] = float(divisors.min())
        row["maximum_histogram_log2"] = float(divisors.max())
    progress["current_complete_run_cover"] = True
    guard()
    return row


def probe(groups=2048, end_groups=2051, z=Fraction(23, 100), alpha=Fraction(19, 200),
          max_cells=128, tail_counts=32, tilt=Fraction(33, 25),
          pseudocount=Fraction(1, 65536), bin_width=Fraction(1, 4),
          seconds=45., memory_mib=1024, target_bits=100, details=False):
    end_groups = groups if end_groups is None else end_groups
    z, alpha, tilt, pseudocount, bin_width = map(Fraction, (z, alpha, tilt, pseudocount, bin_width))
    if (not isinstance(groups, int) or not isinstance(end_groups, int)
            or not 1 <= groups <= end_groups < bridge.GROUPS
            or end_groups-groups+1 > MAX_SOURCE_SIZES
            or not isinstance(max_cells, int) or not 1 <= max_cells <= 128
            or not isinstance(tail_counts, int) or not 0 <= tail_counts <= 32
            or not 0 < z < 1 or not 0 < alpha <= 1 or tilt <= 0 or pseudocount <= 0
            or not 0 < float(z) < 1 or not 0 < float(alpha) <= 1
            or bin_width not in (Fraction(1), Fraction(1, 2), Fraction(1, 4))
            or not isinstance(target_bits, int) or not 1 <= target_bits <= 10000):
        raise ValueError("Require a proper source interval of at most256 sizes, <=128 bulk cells,"
                         "<=32 tail cells, valid positive rational markers, width1,1/2,1/4")
    budget = dense.Budget(seconds, memory_mib)
    anchors = [groups] if groups == end_groups else [groups, groups+1]
    partitions = {anchor: census.cells_for_group(anchor, max_cells, tail_counts) for anchor in anchors}
    progress = dict(completed_family_cells=0, cached_ordered_cell_pairs=0,
                    current_groups=None, current_requested_cells=0, current_completed_cells=0,
                    current_assigned_run_classes=0, current_counted_run_classes=0,
                    current_complete_run_cover=False)
    families, rows, result = {}, [], {}
    status, error, complete = "complete-direct-grouped-interval-float", None, False
    try:
        prepared, receipt = grouped.local_bins(z, alpha, bin_width, budget.guard)
        result["local_kernel"] = receipt
        for anchor in anchors:
            families[anchor % 2] = prepare_family(prepared, anchor, partitions[anchor], alpha,
                                                  tilt, pseudocount, budget.guard, progress)
        for actual in range(groups, end_groups+1):
            budget.guard()
            row = actual_census(actual, families[actual % 2], budget.guard, progress, details)
            source = bridge.source_factor(actual, "joint")
            prefix = float(alpha)*(source["log2_prefactor"]-bridge.LIMIT*math.log2(float(z)))
            bound = row["route_census_log2"]+prefix
            if not math.isfinite(prefix) or not math.isfinite(bound):
                raise ArithmeticError("Nonfinite actual-size source factor or total bound")
            row.update(source_factor_log2=prefix, log2_uncapped_bound=bound)
            rows.append(row)
        with np.errstate(over="raise", invalid="raise", under="ignore"):
            total = bridge.log2_sum(row["log2_uncapped_bound"] for row in rows)
        if len(rows) != end_groups-groups+1 or not math.isfinite(total):
            raise ArithmeticError("Incomplete source interval or nonfinite union bound")
        budget.guard()
        result.update(log2_uncapped_bound=total,
                      worst_source_size=max(rows, key=lambda row: row["log2_uncapped_bound"])["groups"])
        complete = True
    except (TimeoutError, MemoryError, ArithmeticError, RuntimeError) as failure:
        status = str(failure) if isinstance(failure, (TimeoutError, MemoryError)) else "arithmetic-error"
        error = str(failure)
    requested_runs = sum(min(g, bridge.GROUPS-g) for g in range(groups, end_groups+1))
    completed_runs = sum(row["run_classes"] for row in rows)
    requested_cells = sum(len(cells) for cells in partitions.values())
    result.update(schema="small-state-grouped-interval-float-v1", status=status, error=error,
                  first_groups=groups, end_groups=end_groups, z=str(z), alpha=str(alpha),
                  tilt=str(tilt), pseudocount=str(pseudocount), bin_width=str(bin_width),
                  maximum_bulk_cells=max_cells, requested_tail_counts=tail_counts,
                  parity_anchors=anchors, completed_families=len(families),
                  requested_family_cells=requested_cells,
                  unprocessed_family_cells=requested_cells-progress["completed_family_cells"],
                  family_cell_counts={str(anchor): len(partitions[anchor]) for anchor in anchors},
                  requested_source_sizes=end_groups-groups+1, completed_source_sizes=len(rows),
                  unprocessed_source_sizes=end_groups-groups+1-len(rows),
                  requested_run_classes=requested_runs, completed_run_classes=completed_runs,
                  unprocessed_run_classes=requested_runs-completed_runs,
                  progress=progress, source_sizes=rows,
                  complete_profile_cell_cover=complete, complete_source_interval_cover=complete,
                  target_bits=target_bits,
                  diagnostic_within_target=complete and result["log2_uncapped_bound"] < -target_bits,
                  interval_rule="direct actual-g conditioning and census; no subset transport",
                  cache_rule="fixed parity proposals; K_i=W_i^64; all ordered C_ij cached once",
                  cell_rule="shift high-f deficit; extend lowest cell to feasible lower bound; clip and check",
                  histogram_rule="actual-g minimum of two endpoint probabilities, outside alpha",
                  profile_accounting="exact N_g(c) times maximum compatible ordered cell-pair bound",
                  source_accounting="actual-g joint source factor; no extra binomial or anchor normalization",
                  grouping="two ordered raw batches before entrywise alpha",
                  empty_input_rule="raw[[1,1],[256^-16,V_empty]]; omega=1",
                  underflow_rule="positive occupancy mass checked; all operators retained in finite log form",
                  source_mode="joint", directed=False, whole_code_certificate=False,
                  generated_file_bytes=0, output_cap_bytes=OUTPUT_LIMIT, **budget.receipt())
    if details:
        result["families"] = [dict(anchor=family["anchor"], cells=family["reports"])
                              for family in families.values()]
    paths = {Path(module.__file__).resolve() for module in (grouped, census, dense, bridge)}
    paths.add(Path(__file__).resolve())
    provider = sys.modules.get("hd_common_refresh_packets")
    if provider is not None:
        paths.add(Path(provider.__file__).resolve())
    result["dependencies"] = {path.name: hashlib.sha256(path.read_bytes()).hexdigest()
                              for path in sorted(paths)}
    return result


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--groups", type=int, default=2048)
    parser.add_argument("--end-groups", type=int, default=2051)
    parser.add_argument("--z", type=Fraction, default=Fraction(23, 100))
    parser.add_argument("--alpha", type=Fraction, default=Fraction(19, 200))
    parser.add_argument("--max-cells", type=int, default=128)
    parser.add_argument("--tail-counts", type=int, default=32)
    parser.add_argument("--tilt", type=Fraction, default=Fraction(33, 25))
    parser.add_argument("--pseudocount", type=Fraction, default=Fraction(1, 65536))
    parser.add_argument("--bin-width", type=Fraction, default=Fraction(1, 4))
    parser.add_argument("--seconds", type=float, default=45.)
    parser.add_argument("--memory-mib", type=int, default=1024)
    parser.add_argument("--target-bits", type=int, default=100)
    parser.add_argument("--details", action="store_true")
    result = probe(**vars(parser.parse_args()))
    encoded = json.dumps(result, allow_nan=False, sort_keys=True)
    if len(encoded.encode("utf-8")) > OUTPUT_LIMIT:
        result.pop("families", None)
        result["details_omitted_for_output_cap"] = True
        encoded = json.dumps(result, allow_nan=False, sort_keys=True)
    if len(encoded.encode("utf-8")) > OUTPUT_LIMIT:
        raise RuntimeError("Compact direct-interval receipt exceeded its output cap")
    print(encoded, flush=True)
    if not result["complete_source_interval_cover"]:
        raise SystemExit(1)


if __name__ == "__main__":
    main()
