"""Matched exact-rectangle rounding control for actual r776, run count55.

The same coherent current local table supplies both regional calculations.
OLD retains f<=388,s<=188 at markers(341/4096,94/4096); NEW retains
f<=388,s<=110 at markers(360/4096,56/4096). Both contain every target
(f,776-2f) with333<=f<=388. Their unrounded own-normalized coefficients
are identical. Positive256-bit integer rounding can make the reported
upper coefficients differ; no dominance of NEW over OLD is assumed.

For this run class, a+b=721 gives56 ordered profiles and28 shared powers.
Compute both terminal orders and charge the maximum exactly once by the
same indexed run count and actual source factor. Also evaluate entrywise
minima AFTER each candidate's own normalization. Every minimum bounds
the same exact regional entry; it is not an interval intersection.

This completes one run-class majorant only, NOT all r776 source supports
or a whole-code bound. No parameters, local alpha, or source counts change.
Numerical receipts are indexed separately; source alone is not evidence.

Parent-owned serial execution only:180/300s cooperative deadline,512MiB
peak-memory guard,20GiB free,python -B,stdout only. Guards do not interrupt
an in-progress FLINT call. Import performs no numerical work or file writes.
"""

import argparse
from fractions import Fraction as F
import gc
import hashlib
import importlib.util
import json
from math import comb
from pathlib import Path
import shutil
import sys
import time

sys.dont_write_bytecode = True
BASE_PIN = "fec69ce07aa584a63987de825b5cbdde9a507bdd7fdee902317e8d8426e3d115"
WEIGHT, RUNS, POINTS, BITS = 776, 55, 8192, 256
Z, ALPHA = F(984830809, 1600000000), F(119, 1280)
OLD = (294, 388, 188, (F(341, 4096), F(94, 4096)))
NEW = (333, 388, 110, (F(360, 4096), F(56, 4096)))
MEMORY = 512*(1 << 20)


def load_base():
    path = Path(__file__).with_name("cheap_lane_lower_exact_band.py")
    if hashlib.sha256(path.read_bytes()).hexdigest() != BASE_PIN:
        raise RuntimeError("frozen lower-exact helper mismatch")
    spec = importlib.util.spec_from_file_location("tail_frozen_lower_exact", path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def targets():
    degree = WEIGHT-RUNS
    return tuple(range(max(0, degree-WEIGHT//2), min(WEIGHT//2, degree)+1))


def exact_controls(normalizer, join, guard):
    """Own-target cancellation and complete ordered run-class coverage."""
    fulls = targets()
    assert fulls == tuple(range(333, 389))
    ordered = tuple((a, WEIGHT-RUNS-a) for a in fulls)
    shared = tuple((a, b) for a, b in ordered if a <= b)
    assert len(ordered) == 56 and len(shared) == 28
    assert not any(a == b for a, b in ordered)
    assert sorted([pair for a, b in shared for pair in ((a, b), (b, a))]) == sorted(ordered)
    checks = 5
    for f in fulls:
        guard()
        single = WEIGHT-2*f
        mass = comb(4096, f)*comb(4096-f, single)
        coefficient = F(7+2*f+3*single, 31)
        normalized = []
        for lo, hi, cap, marker in (OLD, NEW):
            assert cap == WEIGHT-2*lo and lo <= f <= hi and 0 <= single <= cap
            numerator = coefficient*marker[0]**f*marker[1]**single
            den = normalizer.target_denominator(WEIGHT, f, marker)
            normalized.append(numerator/den)
            assert numerator/den == coefficient/mass
            checks += 2
        assert normalized[0] == normalized[1]
        checks += 1
    counts = join.run_counts(POINTS, WEIGHT)
    formula = F(POINTS, RUNS)*comb(WEIGHT-1, RUNS-1)*comb(POINTS-WEIGHT-1, RUNS-1)
    assert formula.denominator == 1 and counts[RUNS] == formula.numerator
    assert sum(counts) == comb(POINTS, WEIGHT)
    checks += 2
    # Independent exact positive-product check: entrywise minima of two
    # complete upper matrices improve both finite ordered chain bounds.
    old = {}
    new = {}
    minimum = {}
    for f in fulls:
        old[f] = ((F(3+f % 5, 17), F(2, 19), F(5, 23), F(7+f % 3, 29)),
                  (F(11, 31), F(13+f % 4, 37)))
        factors = (F(1, 2), F(3, 2), F(4, 3), F(2, 3), F(5, 4), F(3, 4))
        flat = old[f][0]+old[f][1]
        changed = tuple(x*y for x, y in zip(flat, factors))
        best = tuple(min(x, y) for x, y in zip(flat, changed))
        new[f] = (changed[:4], changed[4:])
        minimum[f] = (best[:4], best[4:])
    for a, b in shared:
        guard()
        oo = join.two_orientations(old[a], old[b], bits=2)
        nn = join.two_orientations(new[a], new[b], bits=2)
        mm = join.two_orientations(minimum[a], minimum[b], bits=2)
        assert all(m <= min(o, n) for m, o, n in zip(mm, oo, nn))
        checks += 2
    return {"tail_fraction_checks": checks, "actual_weight": WEIGHT, "runs": RUNS,
        "ordered_profiles": len(ordered), "shared_powers": len(shared),
        "unrounded_marker_and_rectangle_independence": True,
        "minimum_positive_chain_domination": True, "exact_indexed_run_count_formula": True}


def minimum_regions(v, old, new, guard):
    """Pointwise minimum of two upper majorants, with honest diagnostics."""
    result = {}
    counts = {"old_strictly_smaller": 0, "new_strictly_smaller": 0,
              "overlapping_enclosures": 0, "both_exact_zero": 0,
              "positive_old_zero_new": 0, "zero_old_positive_new": 0}
    ratio_min = ratio_max = None
    checks = 0
    for f in targets():
        guard()
        output = []
        for left, right in zip(old[f], new[f]):
            row = []
            for a, b in zip(left, right):
                value = v.positive(a.min(b))
                if not (value <= a or value.overlaps(a)) or not (value <= b or value.overlaps(b)):
                    raise ArithmeticError("entrywise minimum failed its upper comparisons")
                row.append(value)
                checks += 2
                if a.is_zero() and b.is_zero():
                    counts["both_exact_zero"] += 1
                elif a.is_zero() and b > 0:
                    counts["zero_old_positive_new"] += 1
                elif b.is_zero() and a > 0:
                    counts["positive_old_zero_new"] += 1
                elif a.upper() < b.lower():
                    counts["old_strictly_smaller"] += 1
                elif b.upper() < a.lower():
                    counts["new_strictly_smaller"] += 1
                else:
                    counts["overlapping_enclosures"] += 1
                if a > 0 and b > 0:
                    ratio = (a/b).log()/v.arb(2).log()
                    ratio_min = ratio if ratio_min is None else ratio_min.min(ratio)
                    ratio_max = ratio if ratio_max is None else ratio_max.max(ratio)
            output.append(tuple(row))
        result[f] = tuple(output)
    if sum(counts.values()) != 6*len(targets()):
        raise ArithmeticError("coefficient-comparison coverage mismatch")
    return result, {"compared_entries": 6*len(targets()), "minimum_upper_checks": checks,
        "comparison_counts": counts, "new_dominance_assumed": False,
        "log2_old_over_new_entry_range_intervals": [str(ratio_min), str(ratio_max)]
            if ratio_min is not None else None}


def run_maximum(v, regions, join, guard):
    """One complete run class, before its common source/count factors."""
    if not set(targets()).issubset(regions):
        raise ValueError("missing a marginal needed for this run class")
    degree = WEIGHT-RUNS
    maximum = v.arb(0)
    leader = None
    leader_value = None
    ordered = powers = 0
    for a in targets():
        b = degree-a
        if a > b:
            continue
        guard()
        forward, reverse = join.two_orientations(regions[a], regions[b], diagonal=a == b)
        powers += 1
        for pair, value in (((a, b), forward), ((b, a), reverse)):
            if a == b and pair == (b, a) and ordered % 2:
                raise ArithmeticError("unexpected diagonal in fixed odd degree sum")
            maximum = maximum.max(value)
            if leader_value is None or value.upper() > leader_value.upper():
                leader, leader_value = pair, value
            ordered += 1
    if (ordered, powers) != (56, 28):
        raise ArithmeticError("run-class terminal/profile coverage mismatch")
    return v.positive(maximum), {"ordered_profiles": ordered, "shared_powers": powers,
        "leader_ordered_profile": leader,
        "leader_selection": "largest certified upper endpoint; ties retain first"}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--seconds", type=int, choices=(180, 300), default=300)
    parser.add_argument("--precision", type=int, choices=(256, 384), default=384)
    args = parser.parse_args()
    free = shutil.disk_usage(Path(__file__).resolve().anchor).free
    if free < 20*(1 << 30):
        raise RuntimeError("requires20GiB free on the source volume")
    started = time.monotonic()
    deadline = started+args.seconds
    base = load_base()
    integer = base.load("cheap_lane_integer_regions.py", base.INTEGER_PIN, "tail_integer")
    peak = integer.peak_memory_reader()
    def guard():
        if time.monotonic() > deadline:
            raise TimeoutError("tail rounding deadline; no incomplete comparison")
        if peak() > MEMORY:
            raise MemoryError("tail rounding512MiB peak-memory limit exceeded")
    def progress(phase, **values):
        guard()
        print(json.dumps({"phase": phase, "seconds": time.monotonic()-started,
            "peak_process_bytes": peak(), **values}), flush=True)
    guard()
    anchor = base.load("cheap_lane_cone_cauchy_anchors.py", base.ANCHOR_PIN, "tail_anchor")
    normalizer = base.load("cheap_lane_band_regions.py", base.REGIONS_PIN, "tail_normalizer")
    join = base.load("cheap_lane_directed_run_mass.py", base.JOIN_PIN, "tail_join")
    local = base.load("cheap_lane_ordered_local_directed.py", base.LOCAL_PIN, "tail_local")
    built = local.build_local_pairs(Z, ALPHA, anchor, guard, deadline, args.precision, progress)
    v, core = built["verifier"], built["core"]
    degrees = built["rational_coefficients"]
    before = base.coefficient_hash(degrees)
    if before != built["report"]["rational_coefficient_sha256"]:
        raise ArithmeticError("coherent local table fingerprint mismatch")
    built.pop("coefficients", None)
    tests = {"integer": integer.exact_tests(guard), "normalizer": normalizer.exact_tests(guard),
        "join": join.exact_tests(guard),
        "frozen_wrapper": base.exact_controls(v, anchor, join, normalizer, guard),
        "tail": exact_controls(normalizer, join, guard)}
    candidates, tables = {}, {}
    for label, item in (("old", OLD), ("new", NEW)):
        def emit(phase, **values):
            progress(phase, candidate=label, **values)
        partial, report = base.one_rectangle(v, degrees, item, (WEIGHT,), integer,
                                            normalizer, BITS, guard, emit)
        candidates[label] = {f: partial[WEIGHT][f] for f in targets()}
        tables[label] = report
        del partial
        gc.collect()
        guard()
    if before != base.coefficient_hash(degrees):
        raise ArithmeticError("regional construction changed the shared local table")
    candidates["entrywise_minimum"], comparison = minimum_regions(v, candidates["old"], candidates["new"], guard)
    h, c_h, eta = anchor.source_factor(WEIGHT)
    source_log = v.exact(ALPHA)*(v.exact(eta).log()+WEIGHT*v.exact(core.BETA).log()
                                -core.LIMIT*v.exact(Z).log())
    source = source_log.exp()
    count = join.run_counts(POINTS, WEIGHT)[RUNS]
    logtwo = v.arb(2).log()
    values, reports = {}, {}
    for label, regions in candidates.items():
        maximum, metadata = run_maximum(v, regions, join, guard)
        value = v.positive(source*count*maximum)
        values[label] = value
        reports[label] = {**metadata,
            "run_class_log2_bound_interval": str(value.log()/logtwo) if value > 0 else "zero",
            "conditional_maximum_log2_interval": str(maximum.log()/logtwo) if maximum > 0 else "zero",
            "run_class_below_2_pow_minus_40": bool(value < v.arb(2)**(-40))}
    minimum = values["entrywise_minimum"]
    for label in ("old", "new"):
        if not (minimum <= values[label] or minimum.overlaps(values[label])):
            raise ArithmeticError("refined complete run-class bound failed dominance control")
    gains = {label: str((values["old"]/values[label]).log()/logtwo)
             if values["old"] > 0 and values[label] > 0 else None
             for label in ("new", "entrywise_minimum")}
    guard()
    print(json.dumps({"scope": "one complete actual-r776 run class c55; matched rounding diagnostic",
        "directed": True, "whole_weight_certificate": False, "whole_code_certificate": False,
        "actual_weight": WEIGHT, "runs": RUNS, "degree_sum": WEIGHT-RUNS,
        "required_full_range": [min(targets()), max(targets())],
        "parameters": {"z": str(Z), "alpha": str(ALPHA), "regional_bits": BITS,
            "requested_precision": args.precision, "working_precision": built["report"]["working_precision"],
            "seconds_budget": args.seconds, "memory_budget_bytes": MEMORY, "preflight_free_bytes": free},
        "source_factor": {"h": h, "c_h": str(c_h), "eta_r": str(eta),
            "formula": "(eta_r beta^r z^(-LIMIT))^alpha", "log2_interval": str(source_log/logtwo)},
        "exact_indexed_run_count": str(count),
        "indexed_run_count_sha256": hashlib.sha256(str(count).encode()).hexdigest(),
        "count_outside_alpha": True, "source_applied_once": True,
        "results": reports, "gain_vs_old_bits_intervals": gains,
        "coefficient_comparison": comparison, "tables": tables, "tests": tests,
        "same_unrounded_exact_regional_coefficients": True, "rounding_loss_assumed_zero": False,
        "coefficient_table_unchanged": True,
        "local": {key: built["report"][key] for key in ("rational_coefficient_sha256",
            "raw_table_sha256", "powered_freeze", "tests", "field_tests", "residual_rank_receipt")},
        "nonnegative_intersections": v.intersections, "generated_bytes": 0,
        "elapsed_seconds": time.monotonic()-started, "peak_process_bytes": peak(),
        "resource_guards": "cooperative between operations; no external kill or interruption of FLINT calls",
        "pins": {"lower_exact": BASE_PIN, "integer": base.INTEGER_PIN, "normalizer": base.REGIONS_PIN,
            "join": base.JOIN_PIN, "local": base.LOCAL_PIN, "anchor": base.ANCHOR_PIN},
        "source_sha256": hashlib.sha256(Path(__file__).read_bytes()).hexdigest()}, indent=2), flush=True)


if __name__ == "__main__":
    main()
