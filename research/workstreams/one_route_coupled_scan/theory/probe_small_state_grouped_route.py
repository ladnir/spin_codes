"""Complete FLOAT source census with ratio-binned two-batch fractionalization.

For positive raw local matrices M,N, an entry of their ordered product is
x+y, where x=M[a,0]*N[0,b] and y=M[a,1]*N[1,b]. For 0<alpha<=1,
(x+y)^alpha=c(x/y)*(x^alpha+y^alpha), with
c(r)=(1+r)^alpha/(1+r^alpha). The ratio factors into a left row ratio of M
and a right column ratio of N. Bin these ratios on the log2 scale. Since
c is symmetric under r->1/r and increases with |log2(r)|, the farther
endpoint of each summed bin interval gives a common correction upper bound.
Positive bin sums factor the IID expectation without enumerating 4845^2
occupancy pairs. All 4845 actual striped occupancies remain represented.

Bin indices are computed EXACTLY from positive rational raw local entries;
the moment powers, correction factors, and sums are ordinary floating point.
This is a diagnostic, not a directed certificate. Widths 1,1/2,1/4 are allowed.

W bounds the UNCONDITIONED IID expectation of (M*N)^entrywise-alpha.
Independence gives W^64 for a 128-batch region. Only then divide by the
4096-pair histogram probability, NOT its alpha power. The endpoint minimum
bounds every full-pair count in a cell. Reuse the complete ordered cell-pair
census, exact cycle-run source counts, and actual-size joint source factors.
No files are written. Root owns all runs; use Python -B and external limits.
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
MAX_CORRECTION_ENTRIES = 2_000_000


def ratio_bin(left, right, subdivisions):
    """Exact floor(subdivisions*log2(left/right)), including dyadic edges."""
    if left <= 0 or right <= 0:
        raise ArithmeticError("A raw ratio requires two strictly positive entries")
    numerator = (left.numerator*right.denominator)**subdivisions
    denominator = (left.denominator*right.numerator)**subdivisions
    exponent = numerator.bit_length()-denominator.bit_length()
    below = (numerator < (denominator << exponent) if exponent >= 0
             else (numerator << -exponent) < denominator)
    return exponent-int(below)


def local_bins(z, alpha, bin_width, guard):
    """Prepare raw moments and proposal-independent left/right ratio bins."""
    guard()
    profiles, _ = bridge.occupancy_profiles()
    if len(profiles) != 4845:
        raise ArithmeticError("Incomplete actual sorted occupancy domain")
    subdivisions = bin_width.denominator
    local = bridge.local_moments(z, q=256, h=4)
    epsilon = Fraction(1, 256**16)
    logs = np.empty((len(profiles), 2, 2))
    left_bins = np.empty((2, len(profiles)), dtype=np.int64)
    right_bins = np.empty_like(left_bins)
    for index, profile in enumerate(profiles):
        if (index & 31) == 0:
            guard()
        quiet, birth, ret, live = bridge.matrix(tuple(map(int, profile)), z,
                                                local, q=256, h=4)
        total = quiet+birth
        entries = (quiet, total, epsilon*quiet+(1-epsilon)*ret,
                   epsilon*total+(1-epsilon)*live)
        if min(entries) <= 0 or max(entries) > 1:
            raise ArithmeticError("A raw local rational entry left (0,1]")
        if index == 0 and entries[:3] != (Fraction(1), Fraction(1), epsilon):
            raise ArithmeticError("Empty input did not retain the raw transitions")
        logs[index] = np.asarray([math.log2(value.numerator)-math.log2(value.denominator)
                                  for value in entries]).reshape(2, 2)
        left_bins[0, index] = ratio_bin(entries[0], entries[1], subdivisions)
        left_bins[1, index] = ratio_bin(entries[2], entries[3], subdivisions)
        right_bins[0, index] = ratio_bin(entries[0], entries[2], subdivisions)
        right_bins[1, index] = ratio_bin(entries[1], entries[3], subdivisions)
    if np.any(~np.isfinite(logs)):
        raise ArithmeticError("Nonfinite positive raw local log moment")
    left = [np.unique(row, return_inverse=True) for row in left_bins]
    right = [np.unique(row, return_inverse=True) for row in right_bins]
    total_entries = sum(len(l[0])*len(r[0]) for l in left for r in right)
    if total_entries > MAX_CORRECTION_ENTRIES:
        raise MemoryError("ratio-bin-pair-cap")
    corrections = [[None, None], [None, None]]
    for a in range(2):
        for b in range(2):
            guard()
            lower = left[a][0][:, None]+right[b][0][None, :]
            # Each exact ratio lies in [2^(k/q),2^((k+1)/q)). Use the
            # closed sum interval [lower/q,(lower+2)/q] conservatively.
            distance = np.maximum(np.abs(lower), np.abs(lower+2))/subdivisions
            if alpha == 1:
                correction = np.zeros(distance.shape)
            else:
                # Stable symmetric formula; no subtraction of large logs.
                correction = (float(alpha)*np.logaddexp2(0., -distance)
                              - np.logaddexp2(0., -float(alpha)*distance))
                correction = np.minimum(correction, 0.)
            if np.any(~np.isfinite(correction)) or np.any(correction > 0):
                raise ArithmeticError("Invalid two-path correction bound")
            corrections[a][b] = correction
    guard()
    receipt = dict(occupancies=len(profiles), raw_positive_entries=4*len(profiles),
                   bin_width=str(bin_width), exact_rational_bin_assignment=True,
                   left_bin_counts=[len(row[0]) for row in left],
                   right_bin_counts=[len(row[0]) for row in right],
                   correction_entries=total_entries,
                   correction_entry_cap=MAX_CORRECTION_ENTRIES,
                   grouping="two ordered batches before entrywise alpha")
    return dict(powered=logs*float(alpha), left=left, right=right,
                corrections=corrections), receipt


def regional_cell(prepared, groups, cell, alpha, tilt, pseudocount, guard):
    """Unconditioned grouped IID moment, then one full-region conditioning."""
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
        raise ArithmeticError("Complete positive striped occupancy law did not normalize")
    terms = np.log2(weights[:, None, None])+prepared["powered"]
    left_sums, right_sums = [], []
    with np.errstate(over="raise", invalid="raise", under="ignore"):
        for a in range(2):
            guard()
            keys, inverse = prepared["left"][a]
            values = np.full((2, len(keys)), -math.inf)
            np.logaddexp2.at(values[0], inverse, terms[:, a, 0])
            np.logaddexp2.at(values[1], inverse, terms[:, a, 1])
            left_sums.append(values)
        for b in range(2):
            guard()
            keys, inverse = prepared["right"][b]
            values = np.full((2, len(keys)), -math.inf)
            np.logaddexp2.at(values[0], inverse, terms[:, 0, b])
            np.logaddexp2.at(values[1], inverse, terms[:, 1, b])
            right_sums.append(values)
        grouped = np.empty((2, 2))
        for a in range(2):
            for b in range(2):
                guard()
                left, right = left_sums[a], right_sums[b]
                paths = np.logaddexp2(left[0, :, None]+right[0, None, :],
                                     left[1, :, None]+right[1, None, :])
                paths += prepared["corrections"][a][b]
                grouped[a, b] = np.logaddexp2.reduce(paths.ravel())
        average = np.logaddexp2.reduce(terms, axis=0)
        old_two_batch = dense.log_matrix_product(average, average)
        if np.any(grouped > old_two_batch+1e-10):
            raise ArithmeticError("Grouped upper exceeded the batchwise fractional upper")
        region = grouped.copy()
        if bridge.BATCHES_PER_REGION != 128:
            raise ArithmeticError("Two-batch grouping requires this 128-batch region")
        for _ in range(6):
            guard()
            region = dense.log_matrix_product(region, region)
    endpoints = [bridge.histogram_log_probability(bridge.pair_histogram(groups, full),
                                                 probabilities) for full in cell]
    divisor = min(endpoints)
    region -= divisor
    if np.any(~np.isfinite(region)) or not math.isfinite(divisor):
        raise ArithmeticError("Nonfinite complete regional majorant")
    guard()
    return region, dict(full_pair_low=cell[0], full_pair_high=cell[1],
                        covered_counts=cell[1]-cell[0]+1, tilt=str(tilt),
                        proposal=list(map(str, rational)),
                        endpoint_histogram_log2=endpoints, divisor_log2=divisor,
                        grouped_iid_log2=grouped.tolist(),
                        two_batch_gain_log2=(old_two_batch-grouped).tolist(),
                        regional_log2=region.tolist())


def probe(groups=2048, z=Fraction(23, 100), alpha=Fraction(1, 10), end_groups=None,
          max_cells=128, tail_counts=0, tilt=Fraction(4, 3),
          pseudocount=Fraction(1, 65536), bin_width=Fraction(1),
          seconds=45., memory_mib=1024, details=False):
    end_groups = groups if end_groups is None else end_groups
    z, alpha, tilt, pseudocount, bin_width = map(Fraction, (z, alpha, tilt, pseudocount, bin_width))
    if (not isinstance(groups, int) or not isinstance(end_groups, int)
            or not 1 <= groups <= end_groups <= bridge.GROUPS or groups == bridge.GROUPS
            or not isinstance(max_cells, int) or not 1 <= max_cells <= 128
            or not isinstance(tail_counts, int) or not 0 <= tail_counts <= 32
            or not 0 < z < 1 or not 0 < alpha <= 1 or tilt <= 0 or pseudocount <= 0
            or bin_width not in (Fraction(1), Fraction(1, 2), Fraction(1, 4))
            or not 0 < float(z) < 1 or not 0 < float(alpha) <= 1):
        raise ValueError("Require proper anchor1..8191, anchor<=end<=8192,1..128 bulk cells,"
                         "0..32 tail singletons, positive rational markers, bin width1,1/2,1/4")
    budget = dense.Budget(seconds, memory_mib)
    cells = census.cells_for_group(groups, max_cells, tail_counts)
    regions, reports, result = [], [], {}
    progress = dict(computed_ordered_cell_pairs=0, assigned_run_classes=0,
                    counted_run_classes=0, complete_run_cover=False)
    status, error, complete, completed_sources = "complete-grouped-route-float", None, False, 0
    try:
        prepared, receipt = local_bins(z, alpha, bin_width, budget.guard)
        result["local_kernel"] = receipt
        for index, cell in enumerate(cells):
            region, receipt = regional_cell(prepared, groups, cell, alpha, tilt,
                                            pseudocount, budget.guard)
            regions.append(region)
            reports.append(dict(index=index, **receipt))
        result.update(census.cell_census(groups, cells, regions, budget.guard, progress, details))
        source, source_receipt = dense.source_interval_log(groups, end_groups, z, alpha, budget.guard)
        completed_sources = end_groups-groups+1
        bound = result["raw_anchor_average_operator_log2"]+source
        if not math.isfinite(bound):
            raise ArithmeticError("Nonfinite source-interval bound")
        prefix = float(alpha)*(bridge.source_factor(groups, "joint")["log2_prefactor"]
                               - bridge.LIMIT*math.log2(float(z)))
        result.update(log2_uncapped_bound=bound, interval_source_factor_log2=source,
                      source_interval=source_receipt,
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
    result.update(schema="small-state-grouped-route-float-v1", status=status, error=error,
                  anchor_groups=groups, end_groups=end_groups, z=str(z), alpha=str(alpha),
                  tilt=str(tilt), pseudocount=str(pseudocount), bin_width=str(bin_width),
                  maximum_bulk_cells=max_cells, requested_tail_counts=tail_counts,
                  requested_cells=len(cells), completed_cells=len(regions),
                  unprocessed_cells=len(cells)-len(regions), cell_partition=[list(cell) for cell in cells],
                  requested_run_classes=min(groups, bridge.GROUPS-groups),
                  unprocessed_run_classes=min(groups, bridge.GROUPS-groups)-progress["counted_run_classes"],
                  requested_source_sizes=end_groups-groups+1, completed_source_sizes=completed_sources,
                  unprocessed_source_sizes=end_groups-groups+1-completed_sources,
                  complete_profile_cell_cover=complete, complete_source_interval_cover=complete,
                  census_engine=progress, source_mode="joint",
                  grouping="two ordered raw batches before entrywise alpha",
                  regional_rule="unconditioned grouped IID W^64 / full4096-pair histogram probability",
                  histogram_cell_rule="minimum endpoint probability; not raised to alpha",
                  profile_accounting="exact run counts times maximum compatible ordered cell-pair bound",
                  source_accounting="averaged anchor times sum of actual-r joint factors and binomials",
                  interval_rule="exact grouped raw inclusion monotonicity; fixed z and alpha",
                  empty_input_rule="raw[[1,1],[256^-16,V_empty]]; omega=1",
                  terminal="positive f<=M*ones; no final reset or refresh charge",
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
    parser.add_argument("--groups", type=int, default=2048)
    parser.add_argument("--z", type=Fraction, default=Fraction(23, 100))
    parser.add_argument("--alpha", type=Fraction, default=Fraction(1, 10))
    parser.add_argument("--end-groups", type=int)
    parser.add_argument("--max-cells", type=int, default=128)
    parser.add_argument("--tail-counts", type=int, default=0)
    parser.add_argument("--tilt", type=Fraction, default=Fraction(4, 3))
    parser.add_argument("--pseudocount", type=Fraction, default=Fraction(1, 65536))
    parser.add_argument("--bin-width", type=Fraction, default=Fraction(1))
    parser.add_argument("--seconds", type=float, default=45.)
    parser.add_argument("--memory-mib", type=int, default=1024)
    parser.add_argument("--details", action="store_true")
    args = parser.parse_args()
    result = probe(**vars(args))
    encoded = json.dumps(result, allow_nan=False, sort_keys=True)
    if len(encoded.encode("utf-8")) > OUTPUT_LIMIT:
        result.pop("cells", None)
        result.pop("run_terms", None)
        result["details_omitted_for_output_cap"] = True
        encoded = json.dumps(result, allow_nan=False, sort_keys=True)
    if len(encoded.encode("utf-8")) > OUTPUT_LIMIT:
        raise RuntimeError("Compact grouped summary exceeded its output cap")
    print(encoded, flush=True)
    if not result["complete_profile_cell_cover"]:
        raise SystemExit(1)


if __name__ == "__main__":
    main()
