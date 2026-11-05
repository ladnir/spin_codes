"""Directed positive union for a bounded number of cyclic source runs.

Import hash-pinned independent cone arithmetic and exact indexed source
counts. For r=256 and c runs, enumerate ordered defects a+b=c. Matching
profiles are (f0,s0)=(128-a,2a), (f1,s1)=(128-b,2b). Each indexed support
has exactly one ordered profile. Its finite chain is e0(R0 R1)^127 R0 V1;
the exact profile count already accounts for rotations and orientations.

All state paths are covered. Source supports beyond --max-runs are NOT.
This is not a whole-code certificate. No new cap, route relaxation, or
construction change is introduced. Truncation retains all degrees needed
by every selected profile. A failed interval comparison is reported, not
assumed away.48-second cooperative deadline,20GiB free, stdout only.
"""

import argparse
import hashlib
import importlib.util
import json
from pathlib import Path
import shutil
import sys
import time

sys.dont_write_bytecode = True
CORE_PIN = "7ccf0cdab823ef15fdd5edd6156d4b49d9a000f2ae28c37d5d43f4b6ea341f01"
COUNT_PIN = "25e0d9221a62ffdb79d1d625e6e858ffbb8c0b671078d2bde5d2c6fb27e43ed3"
TEST_PIN = "83e00b6f6d148c8058428295b4132d401df08cfccd20cadcda1c1ba82df7b478"


def load_pinned(filename, digest, name):
    path = Path(__file__).with_name(filename)
    if hashlib.sha256(path.read_bytes()).hexdigest() != digest:
        raise RuntimeError("frozen source hash mismatch: "+filename)
    spec = importlib.util.spec_from_file_location(name, path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def ordered_chain(v, r0, r1, terminal1):
    """One ordered matching profile, with no multiplicity or source factor."""
    value = v.matmul(v.matpow(v.matmul(r0, r1), 127), r0)
    return v.positive(v.matmul(value, terminal1)[0][0])


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--max-runs", type=int, choices=(2, 3, 4, 6, 8, 12, 16), default=2)
    parser.add_argument("--precision", type=int, choices=(192, 256, 384, 512), default=192)
    args = parser.parse_args()
    if shutil.disk_usage(Path(__file__).resolve().anchor).free < 20*(1 << 30):
        raise RuntimeError("requires20GiB free on source volume")
    started = time.monotonic()
    exact_tests = load_pinned("cheap_lane_cone_interval_tests.py", TEST_PIN, "run_verify_tests")
    assert exact_tests.PIN == CORE_PIN
    core = exact_tests.load_verifier()
    counts = load_pinned("cheap_lane_run_profile_counts.py", COUNT_PIN, "run_verify_counts")
    counts.DEADLINE = started+48
    from flint import arb, ctx
    ctx.prec = args.precision
    v = core.Verifier(2*args.max_runs, started+48)
    tests = counts.checks()
    tests["toy_positive_convolution"] = v.toy_test()
    records, route_checks = core.route_records(2*args.max_runs)
    v.guard()
    shapes = sorted({key[2] for key in records})
    caps = v.local_caps(shapes)
    regular, terminal = v.regions(*v.polynomials(records, caps))
    rr = [v.extract(regular, 128-a, 2*a) for a in range(args.max_runs+1)]
    tt = [v.extract(terminal, 128-a, 2*a) for a in range(args.max_runs+1)]
    source = (v.exact(core.ETA)*v.exact(core.BETA)**256
              *v.exact(core.Z)**(-core.LIMIT))**v.exact(core.ALPHA)
    logtwo = arb(2).log()
    rows = []
    total = arb(0)
    total_sources = profiles = 0
    worst, worst_info = None, None
    for runs in range(1, args.max_runs+1):
        subtotal = arb(0)
        source_count = 0
        profile_rows = []
        for a in range(runs+1):
            b = runs-a
            number = counts.profile_count(4096, 256, 128-a, 128-b)
            assert number > 0
            value = number*source*ordered_chain(v, rr[a], rr[b], tt[b])
            subtotal += value
            source_count += number
            profiles += 1
            profile_rows.append({"a": a, "b": b, "indexed_sources": number,
                                 "log2_bound_interval": str(value.log()/logtwo)})
            if worst is None or value.upper() > worst.upper():
                worst = value
                worst_info = {"runs": runs, "a": a, "b": b,
                              "indexed_sources": number,
                              "selection": "largest interval upper endpoint",
                              "log2_bound_interval": str(value.log()/logtwo)}
            v.guard()
        assert source_count == counts.run_total(4096, 256, runs)
        total += subtotal
        total_sources += source_count
        # All profiles are retained in the sum; output remains compact by
        # recording their details only for the first two run classes.
        row = {"runs": runs, "profiles": runs+1, "indexed_sources": source_count,
               "log2_bound_interval": str(subtotal.log()/logtwo)}
        if runs <= 2:
            row["profile_details"] = profile_rows
        rows.append(row)
        if runs == 1:
            control = sum(v.chain([rr[0], rr[1]], [tt[0], tt[1]]), arb(0))
            assert subtotal.overlaps(control)
            tests["consecutive_core_control"] = "overlapping directed intervals"
    v.guard()
    passed40 = bool(total < arb(2)**(-40))
    passed100 = bool(total < arb(2)**(-100))
    print(json.dumps({"scope": "all state paths; size256 source supports with bounded cyclic runs",
        "directed": True, "whole_code_certificate": False,
        "source_supports_outside_run_limit_covered": False,
        "parameters": {"max_runs": args.max_runs, "z": str(core.Z), "alpha": str(core.ALPHA),
            "q": core.Q, "distance_threshold": core.LIMIT, "source_weight": 256,
            "singleton_degree": 2*args.max_runs, "full_degree": core.MAX_FULL},
        "precision_bits": args.precision, "run_bounds": rows, "worst_profile": worst_info,
        "total_profiles": profiles, "total_indexed_sources": total_sources,
        "log2_bound_interval": str(total.log()/logtwo),
        "proved_below_2_pow_minus_40": passed40,
        "proved_below_2_pow_minus_100": passed100,
        "tests": {**tests, "exact_route_normalizers": route_checks,
            "run_total_agreements": args.max_runs, "terminal_updates": 32767},
        "route_records": len(records), "occupancy_shapes": len(shapes),
        "nonnegative_intersections": v.intersections,
        "seconds": time.monotonic()-started, "estimated_peak_bytes": 128*(1 << 20),
        "generated_bytes": 0, "core_pin": CORE_PIN, "count_pin": COUNT_PIN, "test_pin": TEST_PIN,
        "source_sha256": hashlib.sha256(Path(__file__).read_bytes()).hexdigest()}, indent=2))
    if not passed40:
        raise SystemExit("directed interval does not prove the fixed40-bit target")


if __name__ == "__main__":
    main()
