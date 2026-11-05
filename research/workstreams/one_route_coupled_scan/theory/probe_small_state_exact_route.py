"""Exact-coefficient route control for a bounded FLOAT two-state bridge.

Use the same local provider and source factors as screen_small_state_bridge.
Only regional averaging changes: extract the exact striped-pair coefficient
from128 ordered batches, instead of bounding it by an IID tilted sum divided
by its histogram probability. Coefficient arithmetic is still binary64;
"exact" describes the combinatorial coefficient, not directed numerical proof.

By default fix the retained g7 z and three alpha values. Apply alpha entrywise to each
local matrix BEFORE route averaging. Sum the two-state paths in their actual
even/odd region order. The terminal uses the same positive f<=M*ones bound.
No final reset or physical final refresh is imposed. Every candidate pair profile
is tested under both cutoff and joint source corrections; no whole-code
certificate is claimed. No files are written. Use Python -B; root owns runs.
"""

import argparse
from fractions import Fraction
import hashlib
import importlib.util
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


GROUPS = 7
PROFILES = tuple((a, b) for a in range(4) for b in range(4))
Z = Fraction(9634151, 9689600)
ANCHOR_ALPHA = Fraction(5470045, 5470051)
ALPHAS = (Fraction(4, 5)*ANCHOR_ALPHA, ANCHOR_ALPHA, Fraction(1))
FROZEN_HASH = "04b62fddce34dcdf60f14d9969573fbf0818a4755eeffb593845e37a955e9053"
FROZEN_PATH = (Path(__file__).resolve().parents[2] / "char2_field_spin" /
               "theory" / "cheap_lane_rank_debit_gate.py")


def load_occupancy_counts():
    path = FROZEN_PATH if FROZEN_PATH.is_file() else Path(__file__).with_name(FROZEN_PATH.name)
    if hashlib.sha256(path.read_bytes()).hexdigest() != FROZEN_HASH:
        raise RuntimeError("frozen striped-geometry hash mismatch")
    spec = importlib.util.spec_from_file_location("small_state_exact_route_geometry", path)
    frozen = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(frozen)
    return frozen.occupancy_counts


def product(left, right, singles, full, guard):
    """Ordered positive 2x2 polynomial product, truncated only above target degrees."""
    guard()
    result = {}
    for (a, b), matrix in left.items():
        for (c, d), other in right.items():
            if a+c <= singles and b+d <= full:
                key = (a+c, b+d)
                value = matrix @ other
                if key in result:
                    result[key] += value
                else:
                    result[key] = value
    return result


def regional(table, alpha, singles, full, indices, occupancy_counts, guard):
    base = {}
    counts = {}
    for a in range(singles+1):
        for b in range(full+1):
            matrix = np.zeros((2, 2))
            count_sum = 0
            for occupancy, count in occupancy_counts(a, b):
                # Frozen geometry sorts descending; the local provider's
                # index order is ascending. Neither adds another multiplicity.
                index = indices[tuple(sorted(occupancy))]
                matrix += count * np.exp2(alpha * table[index])
                count_sum += count
            base[(a, b)] = matrix
            counts[(a, b)] = count_sum
    value = base
    for _ in range(7):
        value = product(value, value, singles, full, guard)
        # Independent exact integer control of the coefficient normalization.
        following = {}
        for (a, b), left in counts.items():
            for (c, d), right in counts.items():
                if a+c <= singles and b+d <= full:
                    key = (a+c, b+d)
                    following[key] = following.get(key, 0) + left*right
        counts = following
    normalizer = (math.comb(4096, singles) * math.comb(4096-singles, full)
                  * 2**singles)
    if counts[(singles, full)] != normalizer:
        raise ArithmeticError("Exact striped-route coefficient normalization failed")
    result = value[(singles, full)] / normalizer
    if np.any(~np.isfinite(result)) or np.any(result < 0):
        raise ArithmeticError("Invalid floating regional coefficient")
    with np.errstate(divide="ignore"):
        return np.log2(result)


def probe(seconds=30., memory_mib=1024, details=False, groups=7, z=Z, alphas=ALPHAS):
    if (not 0 < seconds <= 30 or not 64 <= memory_mib <= 1024
            or not isinstance(groups, int) or not 1 <= groups <= 32):
        raise ValueError("Require <=30 seconds and <=1024 MiB")
    z, alphas = Fraction(z), tuple(Fraction(alpha) for alpha in alphas)
    if not 0 < z < 1 or not alphas or any(not 0 < alpha <= 1 for alpha in alphas):
        raise ValueError("Require rational z in (0,1) and alphas in (0,1]")
    # The union of the two alternating matchings is one cycle on8192 groups.
    # A proper nonempty subset has at most g-1 induced edges. This excludes
    # only the impossible (g/2,g/2) endpoint when g is even in this range.
    profiles = tuple((a, b) for a in range(groups//2+1) for b in range(groups//2+1)
                     if a+b < groups)
    free = shutil.disk_usage(Path(__file__).resolve().parent).free
    if free < 20*1024**3:
        raise RuntimeError("At least20 GiB free space is required")
    started = time.monotonic()

    def guard():
        if time.monotonic()-started >= seconds:
            raise TimeoutError("exact-route control reached its time cap")
        if bridge._resident_bytes() > memory_mib*1024**2:
            raise MemoryError("exact-route control reached its memory cap")

    occupancy_counts = load_occupancy_counts()
    indices = {tuple(map(int, row)): index
               for index, row in enumerate(bridge.occupancy_profiles()[0])}
    sources = {mode: bridge.source_factor(groups, mode) for mode in ("cutoff", "joint")}
    best, records = {}, []
    normalization_checks = 0
    status = "complete-selected-exact-route-control"
    try:
        table = bridge.local_table(z)["base"]
        guard()
        for exact_alpha in alphas:
            alpha = float(exact_alpha)
            regions = {}
            for full in range(groups//2+1):
                regions[full] = regional(table, alpha, groups-2*full, full,
                                         indices, occupancy_counts, guard)
                normalization_checks += 1
            for f0, f1 in profiles:
                guard()
                pair = bridge.log_matrix_product(regions[f0], regions[f1])
                run = bridge.log_matrix_power(pair, 128)
                operator = bridge.log2_sum(run[0])
                for mode, source in sources.items():
                    prefix = (bridge.log2_binomial(bridge.GROUPS, groups)
                              + alpha*(source["log2_prefactor"]
                                       - bridge.LIMIT*math.log2(float(z))))
                    score = prefix + operator
                    if not math.isfinite(score):
                        raise ArithmeticError("Nonfinite exact-route diagnostic score")
                    row = dict(source_mode=mode, even_full=f0, odd_full=f1,
                        z=str(z), alpha=str(exact_alpha), log2_operator=operator,
                        log2_prefactor_and_source_choices=prefix,
                        log2_uncapped_bound=score)
                    records.append(row)
                    key = (mode, f0, f1)
                    if key not in best or score < best[key]["log2_uncapped_bound"]:
                        best[key] = row
    except (TimeoutError, MemoryError) as error:
        status = str(error)
    summaries = []
    run_counts = {runs: Fraction(8192, runs)*math.comb(groups-1, runs-1)
                  * math.comb(8192-groups-1, runs-1) for runs in range(1, groups+1)}
    if (any(count.denominator != 1 for count in run_counts.values())
            or sum(run_counts.values()) != math.comb(8192, groups)):
        raise ArithmeticError("Exact source-cycle run counts do not partition source sets")
    for mode in sources:
        rows = [row for key, row in best.items() if key[0] == mode]
        if not rows:
            continue
        # Two alternating matching edges induce exactly g-f_even-f_odd
        # occupied runs. Count each run class once, maximizing only over
        # profiles within that class. Each row already includes binom(8192,g).
        run_terms = []
        if len(rows) == len(profiles):
            for runs, count in run_counts.items():
                candidates = [row for row in rows
                              if groups-row["even_full"]-row["odd_full"] == runs]
                if not candidates:
                    raise ArithmeticError("Missing source-cycle run class")
                peak = max(candidates, key=lambda row: row["log2_uncapped_bound"])
                run_terms.append(dict(runs=runs, source_sets=str(count.numerator),
                    even_full=peak["even_full"], odd_full=peak["odd_full"],
                    log2_bound=peak["log2_uncapped_bound"]
                               + math.log2(count.numerator)-math.log2(math.comb(8192, groups))))
        summaries.append(dict(source_mode=mode, checked_profiles=len(rows),
            best_profile=min(rows, key=lambda row: row["log2_uncapped_bound"]),
            worst_profile=max(rows, key=lambda row: row["log2_uncapped_bound"]),
            cycle_run_union_log2=bridge.log2_sum(row["log2_bound"] for row in run_terms)
                                 if run_terms else None,
            worst_run_term=max(run_terms, key=lambda row: row["log2_bound"]) if run_terms else None,
            source_factor=sources[mode]))
    result = dict(schema="small-state-bridge-exact-route-float-v2", status=status,
        actual_groups=groups, z=str(z), alphas=list(map(str, alphas)),
        source_modes=list(sources), local_fractional_grouping="entrywise before route averaging",
        origin_log2_penalty=0, terminal="positive f<=M*ones relaxation; no final reset",
        geometry="exact striped-pair coefficient;128 batches;256 alternating regions",
        requested_profiles=len(profiles),
        complete_profile_cover=all((mode, f0, f1) in best for mode in sources for f0, f1 in profiles),
        complete_selected_grid=len(records) == len(alphas)*len(profiles)*len(sources),
        evaluated_profile_markers=len(records),
        exact_integer_normalization_checks=normalization_checks,
        summaries=summaries,
        frozen_dependency=FROZEN_PATH.name, frozen_sha256=FROZEN_HASH,
        seconds=time.monotonic()-started, peak_resident_bytes=bridge._resident_bytes(),
        free_bytes_at_start=free, directed=False, whole_code_certificate=False,
        generated_file_bytes=0)
    if details:
        result["marker_records"] = records
        result["best_profiles"] = [best[key] for key in sorted(best)]
    return result


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--seconds", type=float, default=30.)
    parser.add_argument("--memory-mib", type=int, default=1024)
    parser.add_argument("--details", action="store_true")
    parser.add_argument("--groups", type=int, default=7)
    parser.add_argument("--z", type=Fraction, default=Z)
    parser.add_argument("--alphas", type=Fraction, nargs="+", default=ALPHAS)
    args = parser.parse_args()
    print(json.dumps(probe(args.seconds, args.memory_mib, args.details,
                           args.groups, args.z, args.alphas),
                     allow_nan=False, sort_keys=True), flush=True)


if __name__ == "__main__":
    main()
