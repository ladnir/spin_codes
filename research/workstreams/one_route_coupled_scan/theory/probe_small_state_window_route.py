"""Bounded FLOAT census with projectively rounded power-of-two batch windows.

Represent a positive row as s*(r,1), and round r UP to a rational geometric
grid. This increases the row coordinatewise. For grid point q and local raw
matrix M, put b=q*M01+M11 and r'=(q*M00+M10)/b. One occupancy sends one grid
state to the ceiling of r', with positive fractional-scale weight b^alpha.
There is no splitting over intermediate quiet/live paths and no duplicated
probability. The first batch supplies two initial rows; the final batch is
not rounded. W_k=nu*T^(k-2)*F bounds the powered k-batch IID expectation.

Precompute only grid-by-4845 destination and coefficient planes. For each
proposal, positive bincount aggregation builds T; propagate two rows, not a
full matrix power. The region is W_k^(128/k), divided ONCE by the complete
4096-pair histogram probability. Reuse the complete source-run census.

The grid points are rational, but kernel arithmetic and ceiling assignment
use binary64. This is NOT a directed certificate. All actual occupancies
are included; no profiles are sampled. No files are written. Root owns all
runs; use Python -B and external limits.
"""

import argparse
from fractions import Fraction
import hashlib
import json
import math
from pathlib import Path
import sys

sys.dont_write_bytecode = True
import probe_small_state_census_anchor as census


dense, np, bridge = census.dense, census.np, census.bridge
OUTPUT_LIMIT = 256 * 1024


def prepare(z, alpha, rho, max_grid, guard):
    guard()
    table = dense.raw_local_table(z, guard)
    if table.shape != (4845, 2, 2) or np.any(~np.isfinite(table)):
        raise ArithmeticError("Incomplete positive raw local table")
    ratios = table[:, :, 0]-table[:, :, 1]
    # An extra lower dyadic and an extra upper grid point keep floating
    # boundary noise away from the finite grid's outer endpoints.
    exponent = math.floor(float(ratios.min()))-1
    point = Fraction(2**exponent) if exponent >= 0 else Fraction(1, 2**(-exponent))
    grid, logs = [], []
    upper_target = float(ratios.max())+math.log2(float(rho))
    while True:
        guard()
        if len(grid) >= max_grid:
            raise MemoryError("projective-grid-cap")
        grid.append(point)
        logs.append(math.log2(point.numerator)-math.log2(point.denominator))
        if logs[-1] >= upper_target:
            break
        point *= rho
    log_grid = np.asarray(logs)
    if np.any(np.diff(log_grid) <= 0):
        raise ArithmeticError("Floating grid lost strict order")
    count, occupancies = len(grid), len(table)
    # Every local entry is <=1. This common shift makes both endpoint
    # coefficient planes <=1 while preserving their common scale exactly
    # in logarithmic bookkeeping.
    shift = float(alpha)*float(np.logaddexp2(0., log_grid[-1]))
    destination = np.empty((count, occupancies), dtype=np.int32)
    scale = np.empty((count, occupancies))
    final_zero = np.empty_like(scale)
    with np.errstate(over="raise", invalid="raise", under="raise"):
        for start in range(0, count, 32):
            guard()
            end = min(count, start+32)
            q = log_grid[start:end, None]
            first = np.logaddexp2(q+table[None, :, 0, 0], table[None, :, 1, 0])
            second = np.logaddexp2(q+table[None, :, 0, 1], table[None, :, 1, 1])
            target = np.searchsorted(log_grid, first-second, side="left")
            if np.any(target < 0) or np.any(target >= count):
                raise ArithmeticError("Projective image escaped the complete grid")
            destination[start:end] = target+np.arange(start, end)[:, None]*count
            np.exp2(float(alpha)*second-shift, out=scale[start:end])
            np.exp2(float(alpha)*first-shift, out=final_zero[start:end])
        initial_power = np.exp2(float(alpha)*table[:, :, 1]).T.copy()
    initial_destination = np.searchsorted(log_grid, ratios.T, side="left")
    if (np.any(initial_destination >= count) or np.any(scale <= 0)
            or np.any(final_zero <= 0) or np.any(initial_power <= 0)
            or np.any(~np.isfinite(scale)) or np.any(~np.isfinite(final_zero))):
        raise ArithmeticError("A positive precomputed coefficient vanished")
    guard()
    receipt = dict(occupancies=occupancies, grid_states=count, grid_cap=max_grid,
                   rho=str(rho), first_grid=str(grid[0]), last_grid=str(grid[-1]),
                   first_grid_log2=float(log_grid[0]), last_grid_log2=float(log_grid[-1]),
                   raw_row_ratio_log2_range=[float(ratios.min()), float(ratios.max())],
                   kernel_log2_shift=shift,
                   precomputed_array_bytes=sum(value.nbytes for value in
                       (destination, scale, final_zero, initial_power, initial_destination)),
                   rational_grid=True, ceiling_assignment="FLOAT; not directed")
    return dict(destination=destination.ravel(), scale=scale, final_zero=final_zero,
                initial_power=initial_power, initial_destination=initial_destination,
                shift=shift, count=count, weighted=np.empty_like(scale)), receipt


def normalize_rows(rows, log_scale):
    maxima = rows.max(axis=1)
    if np.any(maxima <= 0) or np.any(~np.isfinite(maxima)):
        raise ArithmeticError("A positive projective row vanished")
    rows /= maxima[:, None]
    log_scale += np.log2(maxima)


def regional_cell(prepared, groups, cell, alpha, tilt, pseudocount, window, guard):
    guard()
    rational = dense.proposal_for_cell(groups, cell, tilt, pseudocount)
    probabilities = np.asarray(tuple(map(float, rational)))
    if (np.any(probabilities <= 0) or np.any(~np.isfinite(probabilities))
            or abs(float(probabilities.sum())-1.) > 1e-12):
        raise ArithmeticError("Rational proposal conversion lost positive mass")
    with np.errstate(over="raise", invalid="raise", under="raise"):
        weights = bridge.profile_weights(probabilities)
    if (np.any(weights <= 0) or np.any(~np.isfinite(weights))
            or abs(float(weights.sum())-1.) > 1e-10):
        raise ArithmeticError("Complete positive striped proposal did not normalize")
    count = prepared["count"]
    # Exactly count*4845 products, not count*count*4845. This scratch plane
    # is shared by successive proposals, and bincount coalesces destinations.
    with np.errstate(over="raise", invalid="raise", under="raise"):
        np.multiply(prepared["scale"], weights[None, :], out=prepared["weighted"])
        transfer = np.bincount(prepared["destination"], weights=prepared["weighted"].ravel(),
                               minlength=count*count).reshape(count, count)
    guard()
    final = np.empty((count, 2))
    final[:, 0] = prepared["final_zero"] @ weights
    final[:, 1] = prepared["scale"] @ weights
    rows = np.empty((2, count))
    for a in range(2):
        rows[a] = np.bincount(prepared["initial_destination"][a],
                              weights=weights*prepared["initial_power"][a], minlength=count)
    if (np.any(transfer < 0) or np.any(~np.isfinite(transfer))
            or np.any(final <= 0) or np.any(~np.isfinite(final))):
        raise ArithmeticError("Invalid positive projective transfer")
    log_scale = np.zeros(2)
    normalize_rows(rows, log_scale)
    # The first batch has already been included in rows. Do not round the
    # final batch. Keep each internal matrix sum inside its alpha power.
    for _ in range(window-2):
        guard()
        rows = rows @ transfer
        normalize_rows(rows, log_scale)
    values = rows @ final
    if np.any(values <= 0) or np.any(~np.isfinite(values)):
        raise ArithmeticError("Nonfinite projective window moment")
    grouped = (np.log2(values)+log_scale[:, None]
               +(window-1)*prepared["shift"])
    region = bridge.log_matrix_power(grouped, bridge.BATCHES_PER_REGION//window)
    endpoints = [bridge.histogram_log_probability(bridge.pair_histogram(groups, full),
                                                 probabilities) for full in cell]
    divisor = min(endpoints)
    region -= divisor
    if np.any(~np.isfinite(region)) or not math.isfinite(divisor):
        raise ArithmeticError("Nonfinite complete regional bound")
    guard()
    return region, dict(full_pair_low=cell[0], full_pair_high=cell[1],
                        covered_counts=cell[1]-cell[0]+1, proposal=list(map(str, rational)),
                        endpoint_histogram_log2=endpoints, divisor_log2=divisor,
                        window_iid_log2=grouped.tolist(), regional_log2=region.tolist())


def probe(groups=1024, z=Fraction(1, 2), alpha=Fraction(2, 25), end_groups=None,
          max_cells=128, tail_counts=32, tilt=Fraction(6, 5),
          pseudocount=Fraction(1, 65536), window=4, rho=Fraction(6, 5),
          max_grid=1024, seconds=45., memory_mib=1024, details=False):
    end_groups = groups if end_groups is None else end_groups
    z, alpha, tilt, pseudocount, rho = map(Fraction, (z, alpha, tilt, pseudocount, rho))
    if (not isinstance(groups, int) or not isinstance(end_groups, int)
            or not 1 <= groups < bridge.GROUPS or not groups <= end_groups <= bridge.GROUPS
            or not isinstance(max_cells, int) or not 1 <= max_cells <= 128
            or not isinstance(tail_counts, int) or not 0 <= tail_counts <= 32
            or not 0 < z < 1 or not 0 < alpha <= 1 or tilt <= 0 or pseudocount <= 0
            or window not in (4, 8, 16, 32, 64, 128) or rho not in (Fraction(6, 5), Fraction(9, 8))
            or not isinstance(max_grid, int) or not 2 <= max_grid <= 1024
            or not 0 < float(z) < 1 or not 0 < float(alpha) <= 1
            or bridge.BATCHES_PER_REGION != 128):
        raise ValueError("Require proper anchor, valid rational markers, <=128 bulk+32 tail cells,"
                         "window4/8/16/32/64/128, rho6/5 or9/8, and grid cap2..1024")
    budget = dense.Budget(seconds, memory_mib)
    cells = census.cells_for_group(groups, max_cells, tail_counts)
    regions, reports, result = [], [], {}
    progress = dict(computed_ordered_cell_pairs=0, assigned_run_classes=0,
                    counted_run_classes=0, complete_run_cover=False)
    status, error, complete, completed_sources = "complete-window-route-float", None, False, 0
    try:
        prepared, receipt = prepare(z, alpha, rho, max_grid, budget.guard)
        result["local_kernel"] = receipt
        for index, cell in enumerate(cells):
            region, receipt = regional_cell(prepared, groups, cell, alpha, tilt,
                                            pseudocount, window, budget.guard)
            regions.append(region)
            reports.append(dict(index=index, **receipt))
        result.update(census.cell_census(groups, cells, regions, budget.guard, progress, details))
        source, receipt = dense.source_interval_log(groups, end_groups, z, alpha, budget.guard)
        completed_sources = end_groups-groups+1
        bound = result["raw_anchor_average_operator_log2"]+source
        if not math.isfinite(bound):
            raise ArithmeticError("Nonfinite complete source-interval bound")
        prefix = float(alpha)*(bridge.source_factor(groups, "joint")["log2_prefactor"]
                               - bridge.LIMIT*math.log2(float(z)))
        result.update(log2_uncapped_bound=bound, interval_source_factor_log2=source,
                      source_interval=receipt,
                      anchor_log2_uncapped_bound=result["raw_anchor_census_operator_sum_log2"]+prefix)
        for term in result["top_run_terms"]:
            term["anchor_log2_bound"] = term["log2_counted_route_operator"]+prefix
        selected = sorted({term[key] for term in result["top_run_terms"]
                           for key in ("even_cell_index", "odd_cell_index")})
        result["dominant_cell_reports"] = [reports[index] for index in selected]
        budget.guard()
        complete = True
    except (TimeoutError, MemoryError, ArithmeticError, RuntimeError) as failure:
        status = str(failure) if isinstance(failure, (TimeoutError, MemoryError)) else "arithmetic-error"
        error = str(failure)
    result.update(schema="small-state-window-route-float-v1", status=status, error=error,
                  anchor_groups=groups, end_groups=end_groups, z=str(z), alpha=str(alpha),
                  tilt=str(tilt), pseudocount=str(pseudocount), window_batches=window, rho=str(rho),
                  maximum_grid_states=max_grid, maximum_bulk_cells=max_cells,
                  requested_tail_counts=tail_counts, requested_cells=len(cells), completed_cells=len(regions),
                  unprocessed_cells=len(cells)-len(regions), cell_partition=[list(cell) for cell in cells],
                  requested_run_classes=min(groups, bridge.GROUPS-groups),
                  unprocessed_run_classes=min(groups, bridge.GROUPS-groups)-progress["counted_run_classes"],
                  requested_source_sizes=end_groups-groups+1, completed_source_sizes=completed_sources,
                  unprocessed_source_sizes=end_groups-groups+1-completed_sources,
                  complete_profile_cell_cover=complete, complete_source_interval_cover=complete,
                  census_engine=progress, source_mode="joint",
                  projective_rule="coordinatewise upper rounding; one state per history; final batch unrounded",
                  regional_rule="unconditioned W_window^(128/window) / one full histogram probability",
                  source_accounting="complete averaged raw anchor and actual-r joint source factors",
                  interval_rule="exact grouped raw inclusion monotonicity; fixed markers",
                  empty_input_rule="raw[[1,1],[256^-16,V_empty]]; omega=1",
                  terminal="positive f<=M*ones; no final refresh or reset charge",
                  maximum_quantization_tax_log2=float(alpha)*32768*(1-1/window)*math.log2(float(rho)),
                  directed=False, whole_code_certificate=False, generated_file_bytes=0,
                  output_cap_bytes=OUTPUT_LIMIT, **budget.receipt())
    if details:
        result["cells"] = reports
    paths = {Path(module.__file__).resolve() for module in (census, dense, bridge)}
    paths.add(Path(__file__).resolve())
    provider = sys.modules.get("hd_common_refresh_packets")
    if provider is not None:
        paths.add(Path(provider.__file__).resolve())
    result["dependencies"] = {path.name: hashlib.sha256(path.read_bytes()).hexdigest()
                              for path in sorted(paths)}
    return result


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--groups", type=int, default=1024)
    parser.add_argument("--z", type=Fraction, default=Fraction(1, 2))
    parser.add_argument("--alpha", type=Fraction, default=Fraction(2, 25))
    parser.add_argument("--end-groups", type=int)
    parser.add_argument("--max-cells", type=int, default=128)
    parser.add_argument("--tail-counts", type=int, default=32)
    parser.add_argument("--tilt", type=Fraction, default=Fraction(6, 5))
    parser.add_argument("--pseudocount", type=Fraction, default=Fraction(1, 65536))
    parser.add_argument("--window", type=int, choices=(4, 8, 16, 32, 64, 128), default=4)
    parser.add_argument("--rho", type=Fraction, default=Fraction(6, 5))
    parser.add_argument("--max-grid", type=int, default=1024)
    parser.add_argument("--seconds", type=float, default=45.)
    parser.add_argument("--memory-mib", type=int, default=1024)
    parser.add_argument("--details", action="store_true")
    result = probe(**vars(parser.parse_args()))
    encoded = json.dumps(result, allow_nan=False, sort_keys=True)
    if len(encoded.encode("utf-8")) > OUTPUT_LIMIT:
        result.pop("cells", None)
        result.pop("run_terms", None)
        result["details_omitted_for_output_cap"] = True
        encoded = json.dumps(result, allow_nan=False, sort_keys=True)
    if len(encoded.encode("utf-8")) > OUTPUT_LIMIT:
        raise RuntimeError("Compact window receipt exceeded its output cap")
    print(encoded, flush=True)
    if not result["complete_profile_cell_cover"]:
        raise SystemExit(1)


if __name__ == "__main__":
    main()
