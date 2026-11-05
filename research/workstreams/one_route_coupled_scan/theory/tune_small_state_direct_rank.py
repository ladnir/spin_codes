"""Bounded coordinate tuning of the g7 FLOAT source-rank diagnostic.

The objective at each threshold is the maximum over all16 paired profiles.
In log(z),log(x), each profile's positive-polynomial bound is log-convex;
the maximum is convex as well. Fixed golden-section coordinate passes
propose markers. This is not a proof of global optimality or a certificate.

Every proposed z and x is rounded to an explicit rational before evaluation.
The retained schedule initializes all28 thresholds. A deadline preserves
the best complete witness found, even if later coordinate passes stop.
No frozen helper is edited. No files are written. Use Python -B.
"""

import argparse
from collections import OrderedDict
from fractions import Fraction
import json
import math
from pathlib import Path
import shutil
import time

import screen_small_state_direct_rank as gate


PROFILES = tuple((a, b) for a in range(4) for b in range(4))
FIRST = tuple(range(6, 18)) + (26, 27)
ORDER = FIRST + tuple(v for v in range(28) if v not in FIRST)
Z_DENOMINATOR = 1 << 36
X_DENOMINATOR = 1 << 30


def rational_marker(log2_value, denominator, allow_one):
    numerator = round(math.exp2(log2_value) * denominator)
    return Fraction(min(denominator if allow_one else denominator - 1,
                        max(1, numerator)), denominator)


def tune(tau=65280., seconds=50., memory_mib=1024, rounds=2,
         z_steps=10, x_steps=12, thresholds=None, central_z=None):
    if (not 1 < tau <= (1 << 32) - 1 or not 0 < seconds <= 50
            or not 64 <= memory_mib <= 1024 or not 1 <= rounds <= 3
            or not 4 <= z_steps <= 16 or not 4 <= x_steps <= 18):
        raise ValueError("Require bounded markers, <=3 rounds, <=50 seconds and <=1024 MiB")
    thresholds = ORDER if thresholds is None else tuple(thresholds)
    if (not thresholds or len(set(thresholds)) != len(thresholds)
            or any(not isinstance(v, int) or not 0 <= v < 28 for v in thresholds)):
        raise ValueError("Thresholds must be distinct integers0..27")
    central_z = None if central_z is None else Fraction(central_z)
    if central_z is not None and not 0 < central_z < 1:
        raise ValueError("Invalid initial central marker")
    free = shutil.disk_usage(Path(__file__).resolve().parent).free
    if free < 20 * 1024**3:
        raise RuntimeError("At least20 GiB free space is required")
    started = time.monotonic()
    gate.DEADLINE, gate.MEMORY_MIB = started + seconds, memory_mib
    frozen = gate.load_geometry()
    cache = OrderedDict()
    best = {}
    evaluations = 0
    tables = 0
    completed_coordinates = []
    status = "complete-bounded-coordinate-passes"

    def columns(exact_z):
        nonlocal tables
        if exact_z in cache:
            cache.move_to_end(exact_z)
            return cache[exact_z]
        gate.guard()
        z = float(exact_z)
        regions = {f: gate.regional(7 - 2*f, f, tau, z, 32, frozen) for f in range(4)}
        result = {}
        for profile in PROFILES:
            gate.guard()
            half = frozen.column_product(regions[profile[0]], regions[profile[1]])
            result[profile] = frozen.column_product(half, half)
        tables += 1
        cache[exact_z] = result
        while len(cache) > 4:
            cache.popitem(last=False)
        return result

    def evaluate(v, exact_z, exact_x, origin):
        nonlocal evaluations
        gate.guard()
        coefficients = columns(exact_z)
        log2_x = None if v == 27 else math.log2(float(exact_x))
        scores = tuple(gate.threshold_bound(coefficients[p], tau, float(exact_z), v, log2_x)
                       for p in PROFILES)
        if not all(math.isfinite(score) for score in scores):
            raise ArithmeticError("Nonfinite coordinate score")
        evaluations += 1
        worst = max(range(16), key=scores.__getitem__)
        record = dict(v=v, z=str(exact_z), x=None if v == 27 else str(exact_x),
            log2_x=log2_x, max_profile_log2=scores[worst],
            worst_profile=list(PROFILES[worst]), origin=origin, scores=scores)
        if v not in best or record["max_profile_log2"] < best[v]["max_profile_log2"]:
            best[v] = record
        return record["max_profile_log2"]

    def coordinate(v, name, steps, round_index):
        # Each coordinate pass holds the other coordinate fixed. Endpoints
        # are included; a boundary optimum does not trigger range expansion.
        current = best[v]
        exact_z = Fraction(current["z"])
        exact_x = Fraction(current["x"]) if current["x"] is not None else None
        low, high = (-.125, -math.exp2(-20)) if name == "z" else (-8., 0.)
        origin = "coordinate-" + name

        def objective(point):
            if name == "z":
                return evaluate(v, rational_marker(point, Z_DENOMINATOR, False), exact_x, origin)
            return evaluate(v, exact_z, rational_marker(point, X_DENOMINATOR, True), origin)

        objective(low)
        objective(high)
        ratio = (math.sqrt(5.) - 1.) / 2.
        left, right = high - ratio*(high-low), low + ratio*(high-low)
        fl, fr = objective(left), objective(right)
        for _ in range(steps):
            gate.guard()
            if fl <= fr:
                high, right, fr = right, left, fl
                left = high - ratio*(high-low)
                fl = objective(left)
            else:
                low, left, fl = left, right, fr
                right = low + ratio*(high-low)
                fr = objective(right)
        completed_coordinates.append(dict(round=round_index, v=v, coordinate=name))

    try:
        for v, exact_z, log2_x in gate.selected_markers(central_z):
            # Quantize x=1/sqrt2 once; retain the rational z schedule exactly.
            exact_x = (None if log2_x is None else
                       rational_marker(log2_x, X_DENOMINATOR, True))
            evaluate(v, exact_z, exact_x, "retained-initial-marker")
        for round_index in range(1, rounds + 1):
            for v in thresholds:
                coordinate(v, "z", z_steps, round_index)
                if v != 27:
                    coordinate(v, "x", x_steps, round_index)
    except (TimeoutError, MemoryError) as error:
        status = str(error)

    complete = len(best) == 28
    profiles = []
    if complete:
        for index, profile in enumerate(PROFILES):
            scores = [best[v]["scores"][index] for v in range(28)]
            peak = max(scores)
            profiles.append(dict(even_full=profile[0], odd_full=profile[1],
                log2_union_bound=peak+math.log2(sum(math.exp2(score-peak) for score in scores))))
    witnesses = [{key: value for key, value in best[v].items() if key != "scores"}
                 for v in sorted(best)]
    return dict(schema="small-state-direct-rank-coordinate-float-v1", status=status,
        groups=7, coefficient_field_bits=32, tau=tau, distance="7/20",
        objective="maximum of the16 paired-profile bounds at each threshold",
        rounds_requested=rounds, z_steps=z_steps, x_steps=x_steps,
        log2_z_box=[-.125, -math.exp2(-20)], log2_x_box=[-8., 0.],
        tuning_threshold_order=list(thresholds),
        coordinate_passes_completed=len(completed_coordinates),
        completed_coordinates=completed_coordinates,
        rational_denominators=dict(z=Z_DENOMINATOR, x=X_DENOMINATOR),
        central_initial_override=str(central_z) if central_z is not None else None,
        marker_evaluations=evaluations, regional_marker_tables=tables,
        checked_thresholds=len(best), complete_threshold_cover=complete,
        profiles=profiles, thresholds=witnesses,
        max_log2_union_bound=max(row["log2_union_bound"] for row in profiles) if complete else None,
        seconds=time.monotonic()-started, peak_resident_bytes=gate.resident_bytes(),
        free_bytes_at_start=free, directed=False, global_optimum_claim=False,
        whole_code_certificate=False, generated_file_bytes=0)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--tau", type=float, default=65280.)
    parser.add_argument("--seconds", type=float, default=50.)
    parser.add_argument("--memory-mib", type=int, default=1024)
    parser.add_argument("--rounds", type=int, default=2)
    parser.add_argument("--z-steps", type=int, default=10)
    parser.add_argument("--x-steps", type=int, default=12)
    parser.add_argument("--thresholds", type=int, nargs="+")
    parser.add_argument("--central-z", type=Fraction)
    args = parser.parse_args()
    result = tune(args.tau, args.seconds, args.memory_mib, args.rounds,
                  args.z_steps, args.x_steps, args.thresholds, args.central_z)
    print(json.dumps(result, allow_nan=False, sort_keys=True), flush=True)


if __name__ == "__main__":
    main()
