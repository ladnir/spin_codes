"""Reuse frozen exact rectangles for up to129 actual source sizes.

This wrapper changes only orchestration. The pinned lower-exact helper builds
the current ordered pair caps and supplies the rectangle schedule, positive
integer products, exact own-target normalizers, and complete per-size joins.
All requested sizes use the same z/alpha and local coefficient table. Each
size r still has its own target(f,r-2f), eta_r, beta^r, indexed run counts,
and both terminal orders. No cross-size matrix envelope is introduced.

The six regional doublings give H^64 and H^63 J. Rectangle scratch depends
on the largest requested size, not the number of extracted sizes. Stored
small marginals and complete joins do grow with band width. Degree extraction
is exact; positive256-bit integer ceilings can still inflate rare targets.
Markers cancel only before rounding. No bound on rounding loss is assumed.

Only a completed sum covers the requested band. A resource failure does not
produce a partial certificate. --tests-only runs small exact controls and
wide-domain coverage checks; it builds no production local or regional table.
Numerical receipts are indexed separately; source alone is not evidence.

Parent-owned serial execution only. Cooperative300/600/1200-second and512MiB
guards do not interrupt an in-progress FLINT call. At least20GiB free is
required. Use python -B; stdout only, no data files or frozen-source edits.
Import performs no numerical work.
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
ORDINARY_PIN = "cb2cc15b9209a95cdfd57b51d98baed99aa205a82099595075c2e72563b826ed"
MEMORY, POINTS, MAX_SIZES, BITS = 512*(1 << 20), 8192, 129, 256


def load_base():
    path = Path(__file__).with_name("cheap_lane_lower_exact_band.py")
    if hashlib.sha256(path.read_bytes()).hexdigest() != BASE_PIN:
        raise RuntimeError("frozen lower-exact source mismatch")
    spec = importlib.util.spec_from_file_location("wide_frozen_lower_exact", path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def validate_band(first, last):
    if (type(first) is not int or type(last) is not int
            or not 225 <= first <= last <= 1024 or last-first >= MAX_SIZES):
        raise ValueError("require225<=first<=last<=1024 and at most129 actual sizes")
    return tuple(range(first, last+1))


def wide_controls(base, anchor, normalizer, join, guard):
    """Pure rational/index controls; no production polynomial tables."""
    checks = targets = 0
    cases = ((225, 232), (225, 353), (384, 512), (549, 613),
             (615, 743), (640, 768), (896, 1024))
    for first, last in cases:
        guard()
        weights = validate_band(first, last)
        schedule = base.rectangle_schedule(last)
        assert len(weights) == last-first+1 <= MAX_SIZES
        assert [f for lo, hi, _, _ in schedule for f in range(lo, hi+1)] == list(range(last//2+1))
        assert len(schedule) <= base.MAX_TABLES
        checks += 3
        for lo, hi, single, marker in schedule:
            assert single == last-2*lo and min(marker) > 0
            assert 0 <= hi <= 1024 and 0 <= single <= 1024
            assert single*(2*hi+1)+hi+1 <= base.PACKED_LIMIT < 300000
            checks += 3
        for r in weights:
            guard()
            seen = []
            for lo, hi, single, marker in schedule:
                for f in range(lo, min(hi, r//2)+1):
                    s = r-2*f
                    assert 0 <= s <= single and f <= hi
                    assert 2*f+s == r and f+s <= 4096
                    seen.append(f)
            assert seen == list(range(r//2+1))
            targets += len(seen)
            counts = join.run_counts(POINTS, r)
            assert len(counts) == r+1 and sum(counts) == comb(POINTS, r)
            assert counts[0] == 0 and counts[1] == POINTS
            checks += 3
    for first, last in ((224, 225), (1024, 1025), (225, 354),
                        (400, 399), (True, 225), (225, F(232))):
        try:
            validate_band(first, last)
        except ValueError:
            checks += 1
        else:
            raise AssertionError("invalid wide domain accepted")
    # These are exact own-target coefficients, not a neighbor's numerator
    # divided by a different normalizer. Check both endpoint parities and
    # every source-factor discontinuity inside the supported domain.
    samples = (225, 255, 256, 257, 511, 512, 513, 767, 768, 769, 1023, 1024)
    for r in samples:
        guard()
        h, c_h, eta = anchor.source_factor(r)
        expected_h = r//256
        expected_c = 1-F(r, 256*(expected_h+1))
        expected_eta = (1+sum((F(comb(r, j), 255) for j in range(1, expected_h+1)), F(0)))/(255*expected_c)
        assert (h, c_h, eta) == (expected_h, expected_c, expected_eta)
        assert c_h > 0 and eta > 0
        checks += 2
        for f in sorted({0, r//4, r//2}):
            s = r-2*f
            coefficient = F(7+2*f+3*s, 31)
            mass = comb(4096, f)*comb(4096-f, s)
            for marker in ((F(2, 7), F(3, 11)), (F(5, 13), F(7, 17))):
                numerator = coefficient*marker[0]**f*marker[1]**s
                denominator = normalizer.target_denominator(r, f, marker)
                assert numerator/denominator == coefficient/mass
                checks += 1
            if f == r//2:
                assert s == r % 2
                checks += 1
    marker = (F(2, 3), F(3, 5))
    coefficients = {(1, 0): F(2), (1, 1): F(7), (1, 2): F(11)}
    for r in (2, 3, 4):
        s = r-2
        numerator = coefficients[1, s]*marker[0]*marker[1]**s
        denominator = normalizer.target_denominator(r, 1, marker, slots=4)
        assert numerator/denominator == coefficients[1, s]/(comb(4, 1)*comb(3, s))
        if r != 2:
            assert numerator != coefficients[1, 0]*marker[0]
        checks += 1
    return {"wide_index_fraction_checks": checks, "covered_target_indices": targets,
        "tested_bands": cases, "maximum_tested_sizes": MAX_SIZES,
        "odd_even_and_source_jumps": True, "fresh_singleton_coefficient_per_size": True,
        "axis_and_packed_limits_unchanged": True, "production_rectangles_built": False}


def compact_table(report):
    """Keep arithmetic provenance, not a duplicate per-weight coverage dump."""
    keys = ("assigned_full_range", "exact_full_limit", "exact_singleton_limit",
            "u", "v", "arithmetic", "regular_exponent", "terminal_exponent",
            "regular_object", "terminal_object", "own_actual_denominators",
            "positive_aliases", "rounded_coefficients_are_exact")
    result = {key: report[key] for key in keys}
    result["extracted_targets"] = sum(number for _, number in report["coverage_weight_target_counts"])
    return result


def main():
    base = load_base()
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--first", type=int, default=225)
    parser.add_argument("--last", type=int, default=232)
    parser.add_argument("--z", type=base.rational, default=F(267, 311))
    parser.add_argument("--alpha", type=base.rational, default=F(119, 1280))
    parser.add_argument("--bits", type=int, choices=(BITS,), default=BITS)
    parser.add_argument("--precision", type=int, choices=(256, 384), default=384)
    parser.add_argument("--seconds", type=int, choices=(300, 600, 1200), default=300)
    parser.add_argument("--tests-only", action="store_true")
    args = parser.parse_args()
    try:
        weights = validate_band(args.first, args.last)
    except ValueError as error:
        parser.error(str(error))
    if not 0 < args.z < 1 or not 0 < args.alpha <= 1:
        parser.error("require0<z<1 and0<alpha<=1")
    free = shutil.disk_usage(Path(__file__).resolve().anchor).free
    if free < 20*(1 << 30):
        raise RuntimeError("requires20GiB free on the source volume")
    started = time.monotonic()
    deadline = started+args.seconds
    integer = base.load("cheap_lane_integer_regions.py", base.INTEGER_PIN, "wide_integer")
    peak = integer.peak_memory_reader()
    def guard():
        if time.monotonic() > deadline:
            raise TimeoutError("wide exact-band deadline; no incomplete band certificate")
        if peak() > MEMORY:
            raise MemoryError("wide exact-band512MiB peak-memory budget exceeded")
    def progress(phase, **values):
        guard()
        print(json.dumps({"phase": phase, "seconds": time.monotonic()-started,
            "peak_process_bytes": peak(), **values}), flush=True)
    guard()
    anchor = base.load("cheap_lane_cone_cauchy_anchors.py", base.ANCHOR_PIN, "wide_anchor")
    normalizer = base.load("cheap_lane_band_regions.py", base.REGIONS_PIN, "wide_normalizer")
    join = base.load("cheap_lane_directed_run_mass.py", base.JOIN_PIN, "wide_join")
    ordinary = base.load("cheap_lane_ordered_ordinary_directed.py", ORDINARY_PIN, "wide_ordinary_tests")
    # The test verifier does not construct local caps or regional polynomials.
    from flint import ctx
    ctx.prec = max(384, args.precision)
    test_core = anchor.load_adapter().load_verifier()
    test_v = test_core.Verifier(2, deadline)
    tests = {"integer": integer.exact_tests(guard),
        "normalizer": normalizer.exact_tests(guard), "join": join.exact_tests(guard),
        "ordinary": ordinary.exact_tests(guard),
        "frozen_wrapper": base.exact_controls(test_v, anchor, join, normalizer, guard),
        "wide_wrapper": wide_controls(base, anchor, normalizer, join, guard)}
    pins = {"lower_exact": BASE_PIN, "ordinary_controls": ORDINARY_PIN,
        "local": base.LOCAL_PIN, "integer": base.INTEGER_PIN,
        "normalizer": base.REGIONS_PIN, "join": base.JOIN_PIN, "anchor": base.ANCHOR_PIN}
    if args.tests_only:
        guard()
        print(json.dumps({"scope": "exact tests only; no source-size certificate",
            "tests": tests, "production_local_table_built": False,
            "production_regional_tables_built": False, "generated_bytes": 0,
            "seconds": time.monotonic()-started, "peak_process_bytes": peak(),
            "budget_seconds": args.seconds, "memory_budget_bytes": MEMORY,
            "pins": pins, "source_sha256": hashlib.sha256(Path(__file__).read_bytes()).hexdigest()},
            indent=2), flush=True)
        return
    del test_v, test_core
    local = base.load("cheap_lane_ordered_local_directed.py", base.LOCAL_PIN, "wide_local")
    schedule = base.rectangle_schedule(args.last)
    expected_marginals = sum(r//2+1 for r in weights)
    progress("wide_exact_preflight", first=args.first, last=args.last, sizes=len(weights),
        marginals=expected_marginals, rectangles=len(schedule), bits=BITS,
        maximum_packed_length=max(s*(2*hi+1)+hi+1 for _, hi, s, _ in schedule))
    built = local.build_local_pairs(args.z, args.alpha, anchor, guard, deadline,
                                   args.precision, progress=progress)
    v, core = built["verifier"], built["core"]
    degrees = built["rational_coefficients"]
    before = base.coefficient_hash(degrees)
    if before != built["report"]["rational_coefficient_sha256"]:
        raise ArithmeticError("local rational coefficient fingerprint mismatch")
    built.pop("coefficients", None)
    regions = {r: [None]*(r//2+1) for r in weights}
    tables = []
    for item in schedule:
        partial, report = base.one_rectangle(v, degrees, item, weights, integer,
                                            normalizer, BITS, guard, progress)
        for r in weights:
            for f, region in partial[r].items():
                if regions[r][f] is not None:
                    raise ArithmeticError("duplicate actual-weight marginal")
                regions[r][f] = region
        tables.append(compact_table(report))
        progress("wide_rectangle_complete", assigned_full_range=report["assigned_full_range"],
            extracted_targets=tables[-1]["extracted_targets"])
        del partial, report
        gc.collect()
        guard()
    if before != base.coefficient_hash(degrees):
        raise ArithmeticError("local table changed during regional extraction")
    if (sum(map(len, regions.values())) != expected_marginals
            or any(x is None for rows in regions.values() for x in rows)):
        raise ArithmeticError("incomplete wide-band actual marginal coverage")
    del degrees
    logtwo = v.arb(2).log()
    beta_log, z_log = v.exact(core.BETA).log(), v.exact(args.z).log()
    alpha = v.exact(args.alpha)
    total = v.arb(0)
    rows, sources = [], []
    exact_total_count = total_profiles = total_marginals = input_intersections = 0
    all_count_checks = True
    for index, r in enumerate(weights):
        guard()
        h, c_h, eta = anchor.source_factor(r)
        source_log = alpha*(v.exact(eta).log()+r*beta_log-core.LIMIT*z_log)
        source = source_log.exp()
        joined = join.directed_join(regions[r], guard, weight=r, points=POINTS)
        contribution = v.positive(source*joined["bound"])
        if not contribution > 0:
            raise ArithmeticError("positive completed per-size contribution required")
        total += contribution
        exact_total_count += comb(POINTS, r)
        expected_profiles = (r//2+1)**2-int(r % 2 == 0)
        if (joined["report"]["weight"] != r
                or joined["report"]["ordered_profiles"] != expected_profiles
                or not joined["report"]["all_run_counts_sum_to_binomial"]):
            raise ArithmeticError("explicit-weight complete join or count control failed")
        total_profiles += expected_profiles
        total_marginals += len(regions[r])
        input_intersections += joined["report"]["input_nonnegative_intersections"]
        all_count_checks &= joined["report"]["all_run_counts_sum_to_binomial"]
        worst = max(range(1, r+1), key=lambda c: joined["run_bounds"][c-1].upper())
        rows.append([r, str(contribution.log()/logtwo), worst])
        sources.append([r, h, str(c_h), str(eta)])
        progress("complete_wide_actual_weight", weight=r, completed_sizes=index+1,
            ordered_profiles=expected_profiles, worst_runs=worst,
            below_2_pow_minus_40=bool(contribution < v.arb(2)**(-40)))
        del joined, regions[r]
        guard()
    if regions or total_marginals != expected_marginals:
        raise ArithmeticError("wide-band completion bookkeeping failed")
    guard()
    print(json.dumps({"scope": "complete actual source-size band; direct per-size exact run envelopes",
        "directed": True, "all_requested_weights_complete": True, "whole_code_certificate": False,
        "parameters": {"first": args.first, "last": args.last, "sizes": len(weights),
            "z": str(args.z), "alpha": str(args.alpha), "regional_bits": BITS,
            "requested_precision": args.precision, "working_precision": built["report"]["working_precision"],
            "seconds_budget": args.seconds, "memory_budget_bytes": MEMORY,
            "packed_limit": base.PACKED_LIMIT, "preflight_free_bytes": free},
        "log2_band_union_bound_interval": str(total.log()/logtwo),
        "band_below_2_pow_minus_40": bool(total < v.arb(2)**(-40)),
        "band_below_2_pow_minus_60": bool(total < v.arb(2)**(-60)),
        "per_weight_columns": ["actual_weight", "log2_complete_bound_interval", "worst_run_count"],
        "per_weight": rows,
        "source_factor_columns": ["actual_weight", "h", "c_h", "eta_r"], "source_factors": sources,
        "source_formula": "(eta_r * beta^r * z^(-LIMIT))^alpha, once per actual size",
        "exact_counts_outside_alpha": True, "all_own_run_counts_sum_to_binomial": all_count_checks,
        "exact_indexed_count_sum_sha256": hashlib.sha256(str(exact_total_count).encode()).hexdigest(),
        "total_marginals": total_marginals, "total_ordered_profiles": total_profiles,
        "own_target_denominators": True, "regional_input_intersections": input_intersections,
        "regional_order": "e0(Ra Rb)^127 Ra Vb; reverse retains Va",
        "regular_region": "H^64", "terminal_region": "H^63 J;127 updates",
        "new_local_alpha_applied": False, "cross_size_envelope_used": False,
        "tables": tables, "tests": tests,
        "local": {key: built["report"][key] for key in ("rational_coefficient_sha256",
            "raw_table_sha256", "powered_freeze", "tests", "field_tests", "residual_rank_receipt")},
        "coefficient_table_unchanged": True, "nonnegative_intersections": v.intersections,
        "elapsed_seconds": time.monotonic()-started, "peak_process_bytes": peak(),
        "generated_bytes": 0, "resource_guards": "cooperative; no external kill or interruption of FLINT calls",
        "pins": pins, "source_sha256": hashlib.sha256(Path(__file__).read_bytes()).hexdigest()},
        indent=2), flush=True)


if __name__ == "__main__":
    main()
