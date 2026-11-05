"""Bounded FLOAT raw-operator anchor covering every source matching profile.

For one actual source size a, partition the complete feasible full-pair range
into at most128 cells. Each cell has up to five explicit rational IID proposals. Its
regional matrix is the positive IID matrix power divided by the smaller
histogram probability at the cell's two endpoints. The logarithm of that
probability is concave in the full-pair count, so the same matrix bounds every
count in the cell. Taking the entrywise minimum over these valid regional
majorants preserves the bound. Local matrices use the actual striped pair geometry and
entrywise alpha powers BEFORE the IID average.

Unlike the sparse verifier, this helper retains raw empty-input transitions:
[[1,1],[256^-16,V_empty]]. Exact raw moments decrease when the allowed source
set grows. A complete anchor therefore bounds larger source sizes with the
SAME z and alpha, but their own joint-Jensen density and binomial factors.
The optional interval endpoint does not infer a source factor from the anchor.

All ordered cell pairs are covered, including impossible combinations as a
conservative relaxation. This is not a selected-profile screen. It is still
ordinary floating arithmetic, NOT a certificate or a whole-code claim.
No files are written. Run Python -B; root owns every numerical run.
"""

import argparse
from fractions import Fraction
import hashlib
import json
import math
import os
from pathlib import Path
import shutil
import sys
import time

sys.dont_write_bytecode = True
os.environ.setdefault("OPENBLAS_NUM_THREADS", "1")
os.environ.setdefault("OMP_NUM_THREADS", "1")

import numpy as np
import screen_small_state_bridge as bridge


OUTPUT_LIMIT = 128 * 1024


class Budget:
    def __init__(self, seconds, memory_mib):
        if not 0 < seconds <= 50 or not 64 <= memory_mib <= 1024:
            raise ValueError("Require <=50 seconds and64..1024 MiB")
        self.free = shutil.disk_usage(Path(__file__).resolve().parent).free
        if self.free < 20 * 1024**3:
            raise RuntimeError("At least20 GiB free space is required")
        self.started = time.monotonic()
        self.seconds, self.memory_mib = seconds, memory_mib
        self.next_memory_check = self.started

    def guard(self):
        now = time.monotonic()
        if now - self.started >= self.seconds:
            raise TimeoutError("time-cap")
        if now >= self.next_memory_check:
            if bridge._resident_bytes() > self.memory_mib * 1024**2:
                raise MemoryError("memory-cap")
            self.next_memory_check = now + .05

    def receipt(self):
        return dict(seconds=time.monotonic() - self.started,
                    peak_resident_bytes=bridge._resident_bytes(),
                    time_cap_seconds=self.seconds,
                    memory_cap_bytes=self.memory_mib * 1024**2,
                    free_bytes_at_start=self.free)


def log_matrix_product(left, right):
    """Ordered2x2 product with broadcasting, in base-two logarithms."""
    return np.logaddexp2(left[..., :, 0, None] + right[..., None, 0, :],
                        left[..., :, 1, None] + right[..., None, 1, :])


def cells_for_group(groups, max_cells):
    low, high = max(0, groups - bridge.PAIRS), groups // 2
    width = (high - low + max_cells) // max_cells
    cells = [(start, min(high, start + width - 1))
             for start in range(low, high + 1, width)]
    if (not 1 <= len(cells) <= max_cells or cells[0][0] != low
            or cells[-1][1] != high
            or any(left[1] + 1 != right[0] for left, right in zip(cells, cells[1:]))):
        raise ArithmeticError("Incomplete matching-count cell partition")
    return cells


def proposal_for_cell(groups, cell, tilt, pseudocount):
    """Tilt the midpoint histogram plus a positive rational pseudocount."""
    endpoints = [bridge.pair_histogram(groups, full) for full in cell]
    weights = tuple((Fraction(endpoints[0][i] + endpoints[1][i], 2)
                     + pseudocount) * tilt**i for i in range(3))
    total = sum(weights)
    probabilities = tuple(weight / total for weight in weights)
    if min(probabilities) <= 0 or sum(probabilities) != 1:
        raise ArithmeticError("Rational proposal is not a positive distribution")
    return probabilities


def raw_local_table(z, guard):
    guard()
    table = bridge.local_table(z)["base"].copy()
    guard()
    profiles, _ = bridge.occupancy_profiles()
    if not np.array_equal(profiles[0], (0, 0, 0, 0)):
        raise ArithmeticError("Unexpected location of the empty occupancy")
    # These are the exact missing raw entries, not differences of majorants.
    table[0, 0, 1] = 0.
    table[0, 1, 0] = -128.
    if np.any(~np.isfinite(table)):
        raise ArithmeticError("A positive raw local entry vanished")
    return table


def regional_candidate(table, groups, cell, alpha, tilt, pseudocount, guard):
    guard()
    rational = proposal_for_cell(groups, cell, tilt, pseudocount)
    probabilities = np.asarray(tuple(map(float, rational)), dtype=float)
    if (np.any(probabilities <= 0) or np.any(~np.isfinite(probabilities))
            or abs(float(probabilities.sum()) - 1.) > 1e-12):
        raise ArithmeticError("Rational proposal conversion lost positive mass")
    with np.errstate(over="raise", invalid="raise", under="raise"):
        weights = bridge.profile_weights(probabilities)
    if np.any(weights <= 0) or abs(float(weights.sum()) - 1.) > 1e-10:
        raise ArithmeticError("Positive striped occupancy law did not normalize")
    guard()
    # All occupancies have positive proposal mass because the pseudocount is
    # positive. No dropped tiny occupancy or selected local profile is allowed.
    terms = np.log2(weights[:, None, None]) + float(alpha) * table
    with np.errstate(under="ignore"):
        average = np.logaddexp2.reduce(terms, axis=0)
        region = bridge.log_matrix_power(average, bridge.BATCHES_PER_REGION)
    endpoints = [bridge.histogram_log_probability(
        bridge.pair_histogram(groups, full), probabilities) for full in cell]
    divisor = min(endpoints)
    region -= divisor
    if np.any(~np.isfinite(region)) or not math.isfinite(divisor):
        raise ArithmeticError("Nonfinite regional cell majorant")
    guard()
    return region, dict(full_pair_low=cell[0], full_pair_high=cell[1],
                        covered_counts=cell[1] - cell[0] + 1,
                        tilt=str(tilt),
                        proposal=[str(value) for value in rational],
                        endpoint_histogram_log2=endpoints,
                        divisor_log2=divisor)


def regional_cell(table, groups, cell, alpha, tilt, pseudocount, guard,
                  extra_tilts=()):
    candidates, receipts = [], []
    tilts = tuple(dict.fromkeys((tilt, *extra_tilts)))
    for current in tilts:
        guard()
        region, receipt = regional_candidate(table, groups, cell, alpha, current,
                                             pseudocount, guard)
        candidates.append(region)
        receipts.append(receipt)
    if len(candidates) == 1:
        return candidates[0], receipts[0]
    # Each candidate bounds each of the SAME conditional regional entries
    # throughout the cell. No proposal selection depends on sampled setup.
    stacked = np.asarray(candidates)
    selected = np.argmin(stacked, axis=0)
    result = np.min(stacked, axis=0)
    return result, dict(full_pair_low=cell[0], full_pair_high=cell[1],
                        covered_counts=cell[1] - cell[0] + 1,
                        candidates=receipts,
                        selected_entry_tilts=[
                            [str(tilts[int(selected[row, column])]) for column in range(2)]
                            for row in range(2)],
                        rule="entrywise minimum of complete conditional regional majorants")


def source_interval_log(groups, end_groups, z, alpha, guard):
    """Evaluate every actual-g scalar factor, with no endpoint extrapolation."""
    terms = []
    for actual in range(groups, end_groups + 1):
        if (actual - groups) % 64 == 0:
            guard()
        source = bridge.source_factor(actual, "joint")
        score = (bridge.log2_binomial(bridge.GROUPS, actual)
                 + float(alpha) * (source["log2_prefactor"]
                                   - bridge.LIMIT * math.log2(float(z))))
        terms.append(score)
    guard()
    maximum = max(range(len(terms)), key=terms.__getitem__)
    return bridge.log2_sum(terms), dict(
        first_groups=groups, last_groups=end_groups, source_sizes=len(terms),
        maximum_source_factor_groups=groups + maximum,
        maximum_source_factor_log2=terms[maximum],
        source_accounting="joint Jensen at each actual g; binomial outside alpha")


def probe(groups, z, alpha, end_groups=None, max_cells=128, tilt=Fraction(1),
          pseudocount=Fraction(1, 2), seconds=45., memory_mib=1024, details=False,
          extra_tilts=()):
    end_groups = groups if end_groups is None else end_groups
    z, alpha, tilt, pseudocount = map(Fraction, (z, alpha, tilt, pseudocount))
    extra_tilts = tuple(map(Fraction, extra_tilts))
    if (not isinstance(groups, int) or not isinstance(end_groups, int)
            or not 1 <= groups <= end_groups <= bridge.GROUPS
            or not isinstance(max_cells, int) or not 1 <= max_cells <= 128
            or not 0 < z < 1 or not 0 < alpha <= 1 or tilt <= 0 or pseudocount <= 0
            or len(extra_tilts) > 4 or any(value <= 0 for value in extra_tilts)
            or not 0 < float(z) < 1 or not 0 < float(alpha) <= 1):
        raise ValueError("Require1<=anchor<=end<=8192,1..128 cells, rational0<z<1,"
                         "0<alpha<=1, positive tilt and pseudocount, <=4 extra positive tilts")
    budget, cell_reports, result = Budget(seconds, memory_mib), [], {}
    cells = cells_for_group(groups, max_cells)
    status, error = "complete-dense-anchor-float", None
    completed = False
    try:
        table = raw_local_table(z, budget.guard)
        regions = []
        for cell in cells:
            region, receipt = regional_cell(table, groups, cell, alpha, tilt,
                                           pseudocount, budget.guard, extra_tilts)
            regions.append(region)
            cell_reports.append(receipt)
        budget.guard()
        regions = np.asarray(regions)
        # Cell-pair products keep the even/odd chronology. Including pairs
        # with no realizable source set enlarges the maximum, never reduces it.
        with np.errstate(under="ignore"):
            run = log_matrix_product(regions[:, None], regions[None, :])
            for _ in range(7):
                budget.guard()
                run = log_matrix_product(run, run)
            operators = np.logaddexp2(run[..., 0, 0], run[..., 0, 1])
        if np.any(~np.isfinite(operators)):
            raise ArithmeticError("Nonfinite full cell-pair operator")
        peak = np.unravel_index(int(np.argmax(operators)), operators.shape)
        operator = float(operators[peak])
        source, source_receipt = source_interval_log(groups, end_groups, z, alpha,
                                                     budget.guard)
        result.update(raw_anchor_operator_log2=operator,
                      interval_source_factor_log2=source,
                      log2_uncapped_bound=operator + source,
                      worst_cell_pair=dict(even_index=int(peak[0]), odd_index=int(peak[1]),
                          even=cell_reports[peak[0]], odd=cell_reports[peak[1]]),
                      source_interval=source_receipt)
        completed = True
    except (TimeoutError, MemoryError, ArithmeticError, FloatingPointError) as failure:
        status = str(failure) if isinstance(failure, (TimeoutError, MemoryError)) else "arithmetic-error"
        error = str(failure)
    result.update(schema="small-state-dense-anchor-float-v1", status=status, error=error,
                  anchor_groups=groups, end_groups=end_groups,
                  z=str(z), alpha=str(alpha), tilt=str(tilt), pseudocount=str(pseudocount),
                  extra_tilts=list(map(str, extra_tilts)),
                  complete_profile_cell_cover=completed,
                  requested_cells=len(cells), completed_cells=len(cell_reports),
                  matching_count_range=[cells[0][0], cells[-1][1]],
                  ordered_cell_pairs=len(cells)**2 if completed else None,
                  include_impossible_cell_pairs=True,
                  local_fractional_grouping="entrywise before IID stripe expectation",
                  histogram_cell_rule="concave log probability: minimum at endpoints",
                  empty_input_rule="raw[[1,1],[256^-16,V_empty]]; omega=1",
                  terminal="positive f<=M*ones envelope; no final reset or refresh charge",
                  interval_rule="raw exact-moment inclusion monotonicity; constant markers",
                  directed=False, whole_code_certificate=False, generated_file_bytes=0,
                  **budget.receipt())
    paths = {Path(__file__).resolve(), Path(bridge.__file__).resolve()}
    provider = sys.modules.get("hd_common_refresh_packets")
    if provider is not None:
        paths.add(Path(provider.__file__).resolve())
    result["dependencies"] = {path.name: hashlib.sha256(path.read_bytes()).hexdigest()
                              for path in sorted(paths)}
    if details:
        result["cells"] = cell_reports
    return result


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--groups", type=int, required=True)
    parser.add_argument("--z", type=Fraction, required=True)
    parser.add_argument("--alpha", type=Fraction, required=True)
    parser.add_argument("--end-groups", type=int)
    parser.add_argument("--max-cells", type=int, default=128)
    parser.add_argument("--tilt", type=Fraction, default=Fraction(1))
    parser.add_argument("--extra-tilts", type=Fraction, nargs="*", default=())
    parser.add_argument("--pseudocount", type=Fraction, default=Fraction(1, 2))
    parser.add_argument("--seconds", type=float, default=45.)
    parser.add_argument("--memory-mib", type=int, default=1024)
    parser.add_argument("--details", action="store_true")
    args = parser.parse_args()
    result = probe(args.groups, args.z, args.alpha, args.end_groups, args.max_cells,
                   args.tilt, args.pseudocount, args.seconds, args.memory_mib, args.details,
                   args.extra_tilts)
    encoded = json.dumps(result, allow_nan=False, sort_keys=True)
    if len(encoded.encode("utf-8")) > OUTPUT_LIMIT:
        result.pop("cells", None)
        result["details_omitted_for_output_cap"] = True
        encoded = json.dumps(result, allow_nan=False, sort_keys=True)
    if len(encoded.encode("utf-8")) > OUTPUT_LIMIT:
        raise RuntimeError("Compact output cap exceeded")
    print(encoded, flush=True)
    if not result["complete_profile_cell_cover"]:
        raise SystemExit(1)


if __name__ == "__main__":
    main()
