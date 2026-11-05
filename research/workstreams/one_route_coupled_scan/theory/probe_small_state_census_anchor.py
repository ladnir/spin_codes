"""Complete FLOAT cell census and averaged raw anchor, with optional exact tails.

Partition every feasible full-pair count into at most128 bulk cells plus
at most32 singleton tail cells. A cell's regional matrix bounds every count
in that cell, using the existing rational IID proposals and endpoint rule.
Keep both orders of every cell pair; do not replace all regions by one U.

For cells[i]=[l_i,u_i], their integer count pairs meet f_even+f_odd=g-c iff
g-u_i-u_j<=c<=g-l_i-l_j. Intersect this with1..min(g,8192-g); including
cycle-unrealizable pairs is conservative. Compute each ordered cell-pair
operator once, then maximize over pairs covering each c. Exact cycle-run
counts multiply these maxima OUTSIDE alpha. Their exact integer sum must
equal binom(8192,g).

Divide this route-only census sum by binom(8192,g) to obtain an averaged
raw anchor. Raw inclusion monotonicity then extends it to actual sizes
g..end_groups, with the SAME z/alpha and each actual size's source factor.
No source-density factor at the anchor is substituted for a later size.

regional_overrides is an optional dictionary{actual_full_pair_count:2x2
log2_matrix}. Only singleton cells may be overridden. The caller must supply
RAW4096-slot regional upper matrices for the SAME g,z,alpha and local kernel.
This FLOAT interface cannot certify externally supplied values; its receipt
records every override. Default: no overrides. A future directed verifier
must verify their provenance separately.

No numerical proof or whole-code certificate is claimed. Root owns all
runs; use Python -B and external limits. No files are written.
"""

import argparse
from fractions import Fraction
import hashlib
import json
import math
from pathlib import Path
import sys

sys.dont_write_bytecode = True
import probe_small_state_dense_anchor as dense


np, bridge = dense.np, dense.bridge
OUTPUT_LIMIT = 256 * 1024


def cells_for_group(groups, max_cells, tail_counts):
    low, high = max(0, groups-bridge.PAIRS), groups//2
    tail = min(tail_counts, high-low+1)
    bulk_high = high-tail
    cells = []
    if bulk_high >= low:
        width = (bulk_high-low+max_cells)//max_cells
        cells.extend((start, min(bulk_high, start+width-1))
                     for start in range(low, bulk_high+1, width))
    cells.extend((full, full) for full in range(bulk_high+1, high+1))
    if (not 1 <= len(cells) <= max_cells+tail_counts
            or cells[0][0] != low or cells[-1][1] != high
            or any(a[1]+1 != b[0] for a, b in zip(cells, cells[1:]))):
        raise ArithmeticError("Incomplete or overlapping matching-count partition")
    return cells


def prepare_overrides(overrides, cells):
    if overrides is None:
        return {}
    if not isinstance(overrides, dict):
        raise ValueError("regional_overrides must be a dictionary")
    singleton_counts = {low for low, high in cells if low == high}
    result = {}
    for full, matrix in overrides.items():
        if not isinstance(full, int) or full not in singleton_counts:
            raise ValueError("Regional overrides require a singleton cell at the actual full-pair count")
        value = np.asarray(matrix, dtype=float)
        if value.shape != (2, 2) or np.any(~np.isfinite(value)):
            raise ValueError("Override must be a finite2x2 raw regional log matrix")
        result[full] = value.copy()
    return result


def cell_census(groups, cells, regions, guard, progress, details=False):
    """Maximize ordered cell-pair bounds along every integer run-count line."""
    guard()
    count, maximum_runs = len(cells), min(groups, bridge.GROUPS-groups)
    regional = np.asarray(regions)
    if regional.shape != (count, 2, 2) or np.any(~np.isfinite(regional)):
        raise ArithmeticError("Incomplete or nonfinite regional cell table")
    with np.errstate(over="raise", invalid="raise", under="ignore"):
        run = dense.log_matrix_product(regional[:, None], regional[None, :])
        for _ in range(7):
            guard()
            run = dense.log_matrix_product(run, run)
        operators = np.logaddexp2(run[..., 0, 0], run[..., 0, 1])
    if np.any(~np.isfinite(operators)):
        raise ArithmeticError("Nonfinite ordered cell-pair operator")
    progress["computed_ordered_cell_pairs"] = count*count
    lower = np.asarray([cell[0] for cell in cells], dtype=int)
    upper = np.asarray([cell[1] for cell in cells], dtype=int)
    first_run = np.maximum(1, groups-upper[:, None]-upper[None, :])
    last_run = np.minimum(maximum_runs, groups-lower[:, None]-lower[None, :])
    allowed = first_run <= last_run
    candidates = np.flatnonzero(allowed)
    # Each candidate covers a whole integer interval of c. Process candidates
    # in decreasing bound order; assign each c only at its first covering pair.
    order = candidates[np.argsort(operators.ravel()[candidates])[::-1]]
    successor = list(range(maximum_runs+2))
    maxima = np.full(maximum_runs+1, -math.inf)
    maximizing_pair = np.full((maximum_runs+1, 2), -1, dtype=int)

    def next_unfilled(runs):
        while successor[runs] != runs:
            successor[runs] = successor[successor[runs]]
            runs = successor[runs]
        return runs

    assigned = 0
    progress.update(compatible_ordered_cell_pairs=len(candidates), assigned_run_classes=0)
    for index, flat in enumerate(order):
        if (index & 63) == 0:
            guard()
        even, odd = divmod(int(flat), count)
        runs = next_unfilled(int(first_run[even, odd]))
        end = int(last_run[even, odd])
        while runs <= end:
            if (assigned & 63) == 0:
                guard()
            maxima[runs] = operators[even, odd]
            maximizing_pair[runs] = (even, odd)
            assigned += 1
            successor[runs] = next_unfilled(runs+1)
            runs = successor[runs]
        progress["assigned_run_classes"] = assigned
        if assigned == maximum_runs:
            break
    if assigned != maximum_runs or np.any(~np.isfinite(maxima[1:])):
        raise ArithmeticError("Incomplete source-run coverage")
    guard()

    multiplicity, count_sum = bridge.GROUPS, 0
    terms = []
    for runs in range(1, maximum_runs+1):
        if (runs & 63) == 1:
            guard()
        if multiplicity <= 0:
            raise ArithmeticError("Nonpositive exact source-run multiplicity")
        even, odd = map(int, maximizing_pair[runs])
        terms.append(dict(runs=runs, source_sets=str(multiplicity),
                          log2_source_sets=math.log2(multiplicity),
                          log2_route_operator=float(maxima[runs]),
                          log2_counted_route_operator=math.log2(multiplicity)+float(maxima[runs]),
                          even_cell_index=even, odd_cell_index=odd,
                          even_full_pair_cell=list(cells[even]), odd_full_pair_cell=list(cells[odd])))
        count_sum += multiplicity
        if runs < maximum_runs:
            numerator = multiplicity*(groups-runs)*(bridge.GROUPS-groups-runs)
            multiplicity, remainder = divmod(numerator, runs*(runs+1))
            if remainder:
                raise ArithmeticError("Nonintegral cycle-run recurrence")
    expected = math.comb(bridge.GROUPS, groups)
    if count_sum != expected:
        raise ArithmeticError("Exact source-run counts failed to partition source sets")
    total = bridge.log2_sum(row["log2_counted_route_operator"] for row in terms)
    average = total-math.log2(expected)
    if not math.isfinite(total) or not math.isfinite(average):
        raise ArithmeticError("Nonfinite averaged raw anchor")
    guard()
    progress.update(counted_run_classes=maximum_runs, exact_source_count_partition=True,
                    complete_run_cover=True)
    dominant = sorted(terms, key=lambda row: row["log2_counted_route_operator"], reverse=True)[:8]
    result = dict(raw_anchor_census_operator_sum_log2=total,
                  raw_anchor_average_operator_log2=average,
                  exact_source_sets=str(expected), counted_source_sets=str(count_sum),
                  run_classes=maximum_runs, worst_run_term=dominant[0], top_run_terms=dominant)
    if details:
        result["run_terms"] = terms
    return result


def probe(groups, z, alpha, end_groups=None, max_cells=128, tail_counts=0,
          tilts=(Fraction(1),), pseudocount=Fraction(1, 65536),
          seconds=45., memory_mib=1024, details=False, regional_overrides=None):
    end_groups = groups if end_groups is None else end_groups
    z, alpha, pseudocount = map(Fraction, (z, alpha, pseudocount))
    tilts = tuple(dict.fromkeys(map(Fraction, tilts)))
    if (not isinstance(groups, int) or not isinstance(end_groups, int)
            or not 1 <= groups <= end_groups <= bridge.GROUPS or groups == bridge.GROUPS
            or not isinstance(max_cells, int) or not 1 <= max_cells <= 128
            or not isinstance(tail_counts, int) or not 0 <= tail_counts <= 32
            or not 0 < z < 1 or not 0 < alpha <= 1 or pseudocount <= 0
            or not 1 <= len(tilts) <= 5 or any(tilt <= 0 for tilt in tilts)
            or not 0 < float(z) < 1 or not 0 < float(alpha) <= 1):
        raise ValueError("Require proper anchor1..8191, anchor<=end<=8192,1..128 bulk cells,"
                         "0..32 tail singletons, rational markers, and1..5 positive tilts")
    budget = dense.Budget(seconds, memory_mib)
    cells = cells_for_group(groups, max_cells, tail_counts)
    overrides = prepare_overrides(regional_overrides, cells)
    regions, cell_reports, used_overrides, result = [], [], [], {}
    progress = dict(computed_ordered_cell_pairs=0, assigned_run_classes=0,
                    counted_run_classes=0, complete_run_cover=False)
    status, error, complete = "complete-census-anchor-float", None, False
    completed_sources = 0
    try:
        table = None
        for index, cell in enumerate(cells):
            budget.guard()
            if cell[0] == cell[1] and cell[0] in overrides:
                region = overrides[cell[0]]
                receipt = dict(full_pair_low=cell[0], full_pair_high=cell[1], covered_counts=1,
                               regional_override=True, log2_matrix=region.tolist(),
                               provenance="caller-supplied raw regional FLOAT override")
                used_overrides.append(cell[0])
            else:
                if table is None:
                    table = dense.raw_local_table(z, budget.guard)
                region, receipt = dense.regional_cell(table, groups, cell, alpha, tilts[0],
                                                      pseudocount, budget.guard, tilts[1:])
                receipt["regional_override"] = False
            if region.shape != (2, 2) or np.any(~np.isfinite(region)):
                raise ArithmeticError("Invalid complete regional cell majorant")
            regions.append(region)
            cell_reports.append(dict(index=index, **receipt))
        result.update(cell_census(groups, cells, regions, budget.guard, progress, details))
        source, source_receipt = dense.source_interval_log(groups, end_groups, z, alpha, budget.guard)
        completed_sources = end_groups-groups+1
        bound = result["raw_anchor_average_operator_log2"]+source
        if not math.isfinite(bound):
            raise ArithmeticError("Nonfinite source-interval diagnostic bound")
        anchor_source = bridge.source_factor(groups, "joint")
        anchor_prefix = float(alpha)*(anchor_source["log2_prefactor"]-bridge.LIMIT*math.log2(float(z)))
        result.update(log2_uncapped_bound=bound, interval_source_factor_log2=source,
                      source_interval=source_receipt,
                      anchor_log2_uncapped_bound=result["raw_anchor_census_operator_sum_log2"]+anchor_prefix)
        for term in result["top_run_terms"]:
            term["anchor_log2_bound"] = term["log2_counted_route_operator"]+anchor_prefix
        selected = sorted({term[key] for term in result["top_run_terms"]
                           for key in ("even_cell_index", "odd_cell_index")})
        result["dominant_cell_reports"] = [cell_reports[index] for index in selected]
        budget.guard()
        complete = True
    except (TimeoutError, MemoryError, ArithmeticError, RuntimeError) as failure:
        status = str(failure) if isinstance(failure, (TimeoutError, MemoryError)) else "arithmetic-error"
        error = str(failure)
    result.update(schema="small-state-census-anchor-float-v1", status=status, error=error,
                  anchor_groups=groups, end_groups=end_groups, z=str(z), alpha=str(alpha),
                  tilts=list(map(str, tilts)), pseudocount=str(pseudocount),
                  maximum_bulk_cells=max_cells, requested_tail_counts=tail_counts,
                  requested_cells=len(cells), completed_cells=len(regions),
                  unprocessed_cells=len(cells)-len(regions), cell_partition=[list(cell) for cell in cells],
                  requested_run_classes=min(groups, bridge.GROUPS-groups),
                  unprocessed_run_classes=min(groups, bridge.GROUPS-groups)-progress["counted_run_classes"],
                  requested_source_sizes=end_groups-groups+1, completed_source_sizes=completed_sources,
                  unprocessed_source_sizes=end_groups-groups+1-completed_sources,
                  complete_profile_cell_cover=complete, complete_source_interval_cover=complete,
                  census_engine=progress, source_mode="joint",
                  requested_override_counts=sorted(overrides), used_override_counts=used_overrides,
                  override_contract="caller asserts raw4096-slot matrices at this anchor,z,alpha,kernel; FLOAT only",
                  overrides=[dict(full_pairs=full, log2_matrix=overrides[full].tolist())
                             for full in sorted(overrides)],
                  profile_accounting="exact run counts times maximum over compatible ordered cell pairs; no extra binomial",
                  anchor_accounting="route-only census divided by binom(8192,anchor)",
                  source_accounting="sum actual-r binomials and joint density factors outside route average",
                  interval_rule="averaged raw-anchor inclusion monotonicity; common z and alpha",
                  empty_input_rule="raw[[1,1],[256^-16,V_empty]]; omega=1",
                  terminal="positive f<=M*ones envelope; no final reset or refresh charge",
                  directed=False, whole_code_certificate=False, generated_file_bytes=0,
                  output_cap_bytes=OUTPUT_LIMIT, **budget.receipt())
    if details:
        result["cells"] = cell_reports
    paths = {Path(__file__).resolve(), Path(dense.__file__).resolve(), Path(bridge.__file__).resolve()}
    provider = sys.modules.get("hd_common_refresh_packets")
    if provider is not None:
        paths.add(Path(provider.__file__).resolve())
    result["dependencies"] = {path.name: hashlib.sha256(path.read_bytes()).hexdigest()
                              for path in sorted(paths)}
    return result


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--groups", type=int, required=True)
    parser.add_argument("--z", type=Fraction, required=True)
    parser.add_argument("--alpha", type=Fraction, required=True)
    parser.add_argument("--end-groups", type=int)
    parser.add_argument("--max-cells", type=int, default=128)
    parser.add_argument("--tail-counts", type=int, default=0)
    parser.add_argument("--tilts", type=Fraction, nargs="+", default=(Fraction(1),))
    parser.add_argument("--pseudocount", type=Fraction, default=Fraction(1, 65536))
    parser.add_argument("--seconds", type=float, default=45.)
    parser.add_argument("--memory-mib", type=int, default=1024)
    parser.add_argument("--details", action="store_true")
    args = parser.parse_args()
    result = probe(args.groups, args.z, args.alpha, args.end_groups, args.max_cells,
                   args.tail_counts, args.tilts, args.pseudocount,
                   args.seconds, args.memory_mib, args.details)
    encoded = json.dumps(result, allow_nan=False, sort_keys=True)
    if len(encoded.encode("utf-8")) > OUTPUT_LIMIT:
        result.pop("cells", None)
        result.pop("run_terms", None)
        result["details_omitted_for_output_cap"] = True
        encoded = json.dumps(result, allow_nan=False, sort_keys=True)
    if len(encoded.encode("utf-8")) > OUTPUT_LIMIT:
        raise RuntimeError("Compact summary exceeded its output cap")
    print(encoded, flush=True)
    if not result["complete_profile_cell_cover"]:
        raise SystemExit(1)


if __name__ == "__main__":
    main()
