"""Bounded FLOAT continuation of the exact route gate through groups17..32.

Each completed marker covers every matching profile for its actual source
count. Optimize the maximum complete profile bound at one rational z and
alpha; do not mix partially evaluated markers. First cover all requested
counts, then take a few rational coordinate steps and carry neighboring
markers forward. A time cap can leave either coverage or tuning unfinished.

The coefficient engine and rational local formulas are unchanged. Local
moments are evaluated lazily for the occupancies that engine requests.
Optional cycle-run counts diagnose the same completed marker; optimization
always uses the full-binomial profile maximum. No directed claim follows.

Run with Python -B. Root owns all numerical runs. Output is one bounded JSON
object on stdout; no files or persistent caches are written. The internal
time/RSS guards supplement an external timeout and address-space limit.
"""

import argparse
from collections import OrderedDict
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
import probe_small_state_exact_route as exact


bridge = exact.bridge
DEFAULT_GROUPS = tuple(range(17, 33))
# Root's complete g17 FLOAT screen supplies this continuation seed.
Z17 = Fraction(99, 100)
ALPHA17 = Fraction(1, 2)
OUTPUT_LIMIT = 512 * 1024


class SparseLocalTable:
    """The base entries of bridge.local_table, evaluated on first access."""

    def __init__(self, z, occupancies, guard):
        guard()
        self.z, self.occupancies, self.guard = z, occupancies, guard
        self.local = bridge.local_moments(z, q=bridge.Q0, h=4)
        self.values = {}

    def __getitem__(self, index):
        if index in self.values:
            return self.values[index]
        self.guard()
        occupancy = self.occupancies[index]
        quiet, birth, ret, live = bridge.matrix(
            occupancy, self.z, self.local, q=bridge.Q0, h=4)
        epsilon = Fraction(1, bridge.Q0**16)
        total = quiet + birth
        all_live = epsilon*total + (1-epsilon)*live
        entries = ((quiet, total, epsilon*quiet+(1-epsilon)*ret, all_live)
                   if sum(occupancy) else
                   (Fraction(1), Fraction(0), Fraction(0), all_live))
        if min(entries) < 0 or max(entries) > 1:
            raise ArithmeticError("Exact rational local majorant left [0,1]")
        value = np.asarray([math.log2(float(x)) if x else -math.inf
                            for x in entries], dtype=float).reshape((2, 2))
        self.values[index] = value
        return value


def profiles_for(groups):
    return tuple((even, odd) for even in range(groups//2+1)
                 for odd in range(groups//2+1) if even+odd < groups)


def cycle_counts(groups):
    """Exact counts partition the source subsets by their occupied runs."""
    counts = {}
    for runs in range(1, groups+1):
        count = (Fraction(bridge.GROUPS, runs)*math.comb(groups-1, runs-1)
                 * math.comb(bridge.GROUPS-groups-1, runs-1))
        if count.denominator != 1:
            raise ArithmeticError("Nonintegral source-cycle count")
        counts[runs] = count.numerator
    if sum(counts.values()) != math.comb(bridge.GROUPS, groups):
        raise ArithmeticError("Source-cycle counts do not partition source sets")
    return counts


def tune(groups=DEFAULT_GROUPS, seconds=50., memory_mib=1024,
         evaluations=12, sweeps=2, z=Z17, alpha=ALPHA17,
         source_mode="cutoff", run_census=False, seeds=(), details=False):
    groups = tuple(sorted(groups))
    z, alpha = Fraction(z), Fraction(alpha)
    if (not groups or len(set(groups)) != len(groups)
            or any(not isinstance(g, int) or not 17 <= g <= 32 for g in groups)
            or not 0 < seconds <= 50 or not 64 <= memory_mib <= 1024
            or not isinstance(evaluations, int) or not 1 <= evaluations <= 32
            or not isinstance(sweeps, int) or not 0 <= sweeps <= 4
            or source_mode not in ("cutoff", "joint")
            or not 0 < z < 1 or not 0 < alpha <= 1):
        raise ValueError("Require distinct groups17..32, <=50 seconds, <=1024 MiB, "
                         "1..32 evaluations, 0..4 sweeps, and valid rational markers")
    explicit_seeds = {}
    for group, seed_z, seed_alpha in seeds:
        seed_z, seed_alpha = Fraction(seed_z), Fraction(seed_alpha)
        if (group not in groups or group in explicit_seeds
                or not 0 < seed_z < 1 or not 0 < seed_alpha <= 1):
            raise ValueError("Each seed requires one requested group and valid rational markers")
        explicit_seeds[group] = (seed_z, seed_alpha)

    free = shutil.disk_usage(Path(__file__).resolve().parent).free
    if free < 20*1024**3:
        raise RuntimeError("At least20 GiB free space is required")
    started = time.monotonic()
    next_memory_check = started

    def guard():
        nonlocal next_memory_check
        now = time.monotonic()
        if now-started >= seconds:
            raise TimeoutError("time-cap")
        if now >= next_memory_check:
            if bridge._resident_bytes() > memory_mib*1024**2:
                raise MemoryError("memory-cap")
            next_memory_check = now + .05

    best, tested = {}, {g: {} for g in groups}
    finished_sweeps = {g: 0 for g in groups}
    normalization_checks, attempted = 0, 0
    status, interruption = "complete-bounded-search", None
    local_tables, geometry_rows = OrderedDict(), {}
    active = None

    def evaluate(group, marker_z, marker_alpha, reason):
        nonlocal normalization_checks, attempted, active
        key = (marker_z, marker_alpha)
        if key in tested[group]:
            return tested[group][key]
        if len(tested[group]) >= evaluations:
            return None
        guard()
        active = dict(groups=group, z=str(marker_z), alpha=str(marker_alpha), reason=reason)
        attempted += 1
        evaluation_start = time.monotonic()
        if marker_z not in local_tables:
            # Keep only three temporary tables: the center and two z neighbors.
            if len(local_tables) == 3:
                local_tables.popitem(last=False)
            local_tables[marker_z] = SparseLocalTable(marker_z, occupancies, guard)
        local_tables.move_to_end(marker_z)
        table = local_tables[marker_z]
        regions = {}
        for full in range(group//2+1):
            guard()
            regions[full] = exact.regional(table, float(marker_alpha), group-2*full,
                                          full, indices, occupancy_counts, guard)
            normalization_checks += 1
        source = sources[group]
        source_choices = math.log2(math.comb(bridge.GROUPS, group))
        prefix = source_choices + float(marker_alpha)*(source["log2_prefactor"]
                         - bridge.LIMIT*math.log2(float(marker_z)))
        rows = []
        for even, odd in profiles[group]:
            guard()
            run = bridge.log_matrix_power(
                bridge.log_matrix_product(regions[even], regions[odd]), 128)
            operator = bridge.log2_sum(run[0])
            score = prefix + operator
            if not math.isfinite(score):
                raise ArithmeticError("Nonfinite complete profile score")
            rows.append(dict(even_full=even, odd_full=odd,
                             log2_uncapped_bound=score))
        worst = max(rows, key=lambda row: row["log2_uncapped_bound"])
        result = dict(groups=group, z=str(marker_z), alpha=str(marker_alpha),
                      source_mode=source_mode, reason=reason,
                      checked_profiles=len(rows), complete_profile_cover=True,
                      worst_profile=worst,
                      log2_uncapped_bound=worst["log2_uncapped_bound"],
                      seconds=time.monotonic()-evaluation_start)
        if run_census:
            run_terms = []
            for runs, count in run_counts[group].items():
                candidates = [row for row in rows
                              if group-row["even_full"]-row["odd_full"] == runs]
                if not candidates:
                    raise ArithmeticError("Missing source-cycle run class")
                peak = max(candidates, key=lambda row: row["log2_uncapped_bound"])
                run_terms.append(dict(runs=runs, source_sets=str(count),
                    even_full=peak["even_full"], odd_full=peak["odd_full"],
                    log2_bound=peak["log2_uncapped_bound"]+math.log2(count)-source_choices))
            result["cycle_run_union_log2"] = bridge.log2_sum(
                row["log2_bound"] for row in run_terms)
            result["worst_run_term"] = max(run_terms, key=lambda row: row["log2_bound"])
        if details:
            result["profiles"] = rows
        guard()
        # Publish a witness only after every region and profile has completed.
        tested[group][key] = result
        if group not in best or result["log2_uncapped_bound"] < best[group]["log2_uncapped_bound"]:
            best[group] = result
        active = None
        return result

    try:
        guard()
        raw_counts = exact.load_occupancy_counts()
        occupancies = tuple(tuple(map(int, row)) for row in bridge.occupancy_profiles()[0])
        indices = {row: i for i, row in enumerate(occupancies)}
        profiles = {g: profiles_for(g) for g in groups}
        sources = {g: bridge.source_factor(g, source_mode) for g in groups}
        run_counts = {g: cycle_counts(g) for g in groups} if run_census else {}

        def occupancy_counts(singles, full):
            # At most289 coefficient indices occur for groups<=32.
            key = (singles, full)
            if key not in geometry_rows:
                guard()
                if singles+2*full > 32 or len(geometry_rows) >= 289:
                    raise ArithmeticError("Sparse geometry exceeded its bounded domain")
                geometry_rows[key] = raw_counts(singles, full)
            return geometry_rows[key]

        carry = (z, alpha)
        for group in groups:
            seed = explicit_seeds.get(group, carry)
            evaluate(group, *seed, "explicit-seed" if group in explicit_seeds else "neighbor-seed")
            carry = (Fraction(best[group]["z"]), Fraction(best[group]["alpha"]))

        # Sweep the coordinates, then move to the next source count. This carries
        # each newly improved witness into the neighboring count on the same pass.
        for sweep in range(sweeps):
            for position, group in enumerate(groups):
                if len(tested[group]) >= evaluations:
                    continue
                if position:
                    neighbor = best[groups[position-1]]
                    evaluate(group, Fraction(neighbor["z"]), Fraction(neighbor["alpha"]),
                             "carry-neighbor")
                if len(tested[group]) >= evaluations:
                    continue
                center_z = Fraction(best[group]["z"])
                center_alpha = Fraction(best[group]["alpha"])
                radius_z = (1-center_z)/2**(sweep+3)
                complete_sweep = True
                for direction in (-1, 1):
                    candidate_z = center_z + direction*radius_z
                    if 0 < candidate_z < 1:
                        if evaluate(group, candidate_z, center_alpha, "z-coordinate") is None:
                            complete_sweep = False
                center_z = Fraction(best[group]["z"])
                center_alpha = Fraction(best[group]["alpha"])
                radius_alpha = Fraction(1, 10*2**sweep)
                for direction in (-1, 1):
                    candidate_alpha = center_alpha + direction*radius_alpha
                    if 0 < candidate_alpha <= 1:
                        if evaluate(group, center_z, candidate_alpha, "alpha-coordinate") is None:
                            complete_sweep = False
                if complete_sweep:
                    finished_sweeps[group] += 1
        guard()
    except (TimeoutError, MemoryError, ArithmeticError) as error:
        status = str(error) if isinstance(error, (TimeoutError, MemoryError)) else "arithmetic-error"
        interruption = dict(reason=str(error), incomplete_marker=active)

    reports = []
    for group in groups:
        report = dict(groups=group, completed_markers=len(tested[group]),
                      requested_profiles=len(profiles_for(group)),
                      coordinate_sweeps_completed=finished_sweeps[group],
                      best=best.get(group))
        if details:
            report["marker_summaries"] = [
                {key: value for key, value in row.items() if key != "profiles"}
                for row in tested[group].values()]
        reports.append(report)
    complete = len(best) == len(groups)
    source_path = Path(__file__).resolve()
    return dict(schema="small-state-route-range-float-v1", status=status,
                requested_groups=list(groups), completed_groups=sorted(best),
                unprocessed_groups=[g for g in groups if g not in best],
                complete_source_count_cover=complete,
                complete_requested_range_17_to_32=complete and groups == DEFAULT_GROUPS,
                diagnostic_range_union_log2=bridge.log2_sum(
                    best[g]["log2_uncapped_bound"] for g in groups) if complete else None,
                source_mode=source_mode, reports=reports, interruption=interruption,
                attempted_markers=attempted,
                completed_markers=sum(map(len, tested.values())),
                exact_integer_normalization_checks=normalization_checks,
                sparse_geometry_coefficients=len(geometry_rows),
                resident_local_tables=len(local_tables),
                resident_local_entries=sum(len(table.values) for table in local_tables.values()),
                objective="maximum full-binomial profile score at one common marker per g",
                local_fractional_grouping="entrywise before route averaging",
                terminal="positive f<=M*ones envelope; no final reset or final refresh charge",
                geometry="exact striped-pair coefficient;128 batches;256 alternating regions",
                marker_selection="baseline complete cover, neighbor continuation, rational coordinate steps",
                initial_z=str(z), initial_alpha=str(alpha),
                maximum_evaluations_per_group=evaluations, requested_sweeps=sweeps,
                run_census_diagnostic=run_census, seconds=time.monotonic()-started,
                peak_resident_bytes=bridge._resident_bytes(), free_bytes_at_start=free,
                time_cap_seconds=seconds, memory_cap_bytes=memory_mib*1024**2,
                output_cap_bytes=OUTPUT_LIMIT,
                dependencies={path.name: hashlib.sha256(path.read_bytes()).hexdigest()
                              for path in (source_path, Path(exact.__file__), Path(bridge.__file__))},
                frozen_geometry_sha256=exact.FROZEN_HASH,
                directed=False, whole_code_certificate=False, generated_file_bytes=0)


def parse_seed(value):
    try:
        group, z, alpha = value.split(":")
        return int(group), Fraction(z), Fraction(alpha)
    except (ValueError, ZeroDivisionError) as error:
        raise argparse.ArgumentTypeError("Use GROUP:RATIONAL_Z:RATIONAL_ALPHA") from error


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--groups", type=int, nargs="+", default=DEFAULT_GROUPS)
    parser.add_argument("--seconds", type=float, default=50.)
    parser.add_argument("--memory-mib", type=int, default=1024)
    parser.add_argument("--evaluations", type=int, default=12)
    parser.add_argument("--sweeps", type=int, default=2)
    parser.add_argument("--z", type=Fraction, default=Z17)
    parser.add_argument("--alpha", type=Fraction, default=ALPHA17)
    parser.add_argument("--source-mode", choices=("cutoff", "joint"), default="cutoff")
    parser.add_argument("--run-census", action="store_true")
    parser.add_argument("--seed", type=parse_seed, action="append", default=[])
    parser.add_argument("--details", action="store_true")
    args = parser.parse_args()
    result = tune(args.groups, args.seconds, args.memory_mib, args.evaluations,
                  args.sweeps, args.z, args.alpha, args.source_mode,
                  args.run_census, args.seed, args.details)
    encoded = json.dumps(result, allow_nan=False, sort_keys=True)
    if len(encoded.encode("utf-8")) > OUTPUT_LIMIT:
        for report in result["reports"]:
            report.pop("marker_summaries", None)
            if report["best"] is not None:
                report["best"].pop("profiles", None)
        result["details_omitted_for_output_cap"] = True
        encoded = json.dumps(result, allow_nan=False, sort_keys=True)
    if len(encoded.encode("utf-8")) > OUTPUT_LIMIT:
        raise RuntimeError("Compact summary exceeded its output cap")
    print(encoded, flush=True)
    if not result["complete_source_count_cover"] or result["status"] == "arithmetic-error":
        raise SystemExit(1)


if __name__ == "__main__":
    main()
