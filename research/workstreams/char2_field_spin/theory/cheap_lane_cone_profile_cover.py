"""Directed cover of ALL size256 source supports, using four coefficient bands.

The local cone and carry-safe polynomial arithmetic are hash-pinned.
Geometry, local caps and local alpha powers are computed once. Four
downward-closed rectangles recover all129 marginal matching profiles
(full,single)=(128-a,2a), a=0..128. Exact dyadic scales are independently
removed before profiles from different bands are composed.

For each run count c, the default bound uses exact indexed profile counts.
Also compute the total c-run count times the maximum conditional chain.
This second bound is deliberately conservative; the maximum is NOT a
weighted average and no logarithmic-loss guarantee is assumed. Positive
sums cover all c=1..256. Exact profile totals are checked against each
run total, then against C(8192,256). The common source/distance factor is
powered once. --conservative-counts selects the second bound for the test.

This covers every actual state path at source weight256, not other source
weights and not the whole code.180-second cooperative deadline,20GiB free,
estimated384MiB peak, bounded stdout progress, no generated files.
"""

import argparse
from fractions import Fraction as F
import hashlib
import importlib.util
import json
from math import comb
from pathlib import Path
import shutil
import sys
import time

sys.dont_write_bytecode = True
CORE_PIN = "7ccf0cdab823ef15fdd5edd6156d4b49d9a000f2ae28c37d5d43f4b6ea341f01"
TEST_PIN = "83e00b6f6d148c8058428295b4132d401df08cfccd20cadcda1c1ba82df7b478"
PACK_PIN = "9379de72b355d04e646a52bd560f538d77cf7095e0e8785b7e5bf982a1292293"
COUNT_PIN = "4ce61f179a8bf41e5d769e955a9efeb437e46af81f57b0d784c4cfc57958d8e8"
BANDS = ((0, 32, F(1, 32), F(1, 128)),
         (33, 64, F(5, 256), F(3, 128)),
         (65, 96, F(3, 256), F(5, 128)),
         (97, 128, F(1, 256), F(7, 128)))
REPORT_RUNS = (1, 2, 8, 16, 32, 64, 128, 256)


def load(filename, digest, name):
    path = Path(__file__).with_name(filename)
    if hashlib.sha256(path.read_bytes()).hexdigest() != digest:
        raise RuntimeError("frozen source hash mismatch: "+filename)
    spec = importlib.util.spec_from_file_location(name, path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def local_polynomials(v, core, records, powered):
    """Reuse local alpha powers; only exact variable scales change by band."""
    fmax, smax = min(32, v.packed_full), min(32, v.single)
    values = [[[v.arb(0) for _ in range(fmax+1)] for _ in range(smax+1)]
              for _ in range(6)]
    scales = {(f, s): v.exact(core.USCALE**f*core.VSCALE**s)
              for f in range(fmax+1) for s in range(smax+1)}
    for (f, s, occ), weight in records.items():
        if f <= fmax and s <= smax:
            factor = v.exact(weight)*scales[f, s]
            for j, cap in enumerate(powered[occ]):
                values[j][s][f] += factor*cap
    entries = []
    for item in values:
        entries.append([v.poly(row) for row in item]
                       +[v.poly([]) for _ in range(v.single-smax)])
    return [[entries[0], entries[1]], [entries[2], entries[3]]], [[entries[4]], [entries[5]]]


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--precision", type=int, choices=(256, 384, 512), default=384)
    mode = parser.add_mutually_exclusive_group()
    mode.add_argument("--exact-counts", dest="exact_counts", action="store_true")
    mode.add_argument("--conservative-counts", dest="exact_counts", action="store_false")
    parser.set_defaults(exact_counts=True)
    args = parser.parse_args()
    if shutil.disk_usage(Path(__file__).resolve().anchor).free < 20*(1 << 30):
        raise RuntimeError("requires20GiB free on source volume")
    started = time.monotonic(); deadline = started+180

    def progress(phase, **values):
        print(json.dumps({"phase": phase, "seconds": time.monotonic()-started, **values}), flush=True)

    adapter = load("cheap_lane_cone_interval_tests.py", TEST_PIN, "cover_interval_tests")
    assert adapter.PIN == CORE_PIN
    core = adapter.load_verifier()
    packed_module = load("cheap_lane_cone_packed.py", PACK_PIN, "cover_packed")
    assert packed_module.CORE_PIN == CORE_PIN
    counts = load("cheap_lane_fast_profile_counts.py", COUNT_PIN, "cover_counts")
    counts.DEADLINE = deadline
    from flint import arb, ctx
    ctx.prec = args.precision

    class CoverVerifier(packed_module.packed_class(core)):
        def guard(self):
            if time.monotonic() > self.deadline:
                raise TimeoutError("180-second profile-cover budget exceeded")

    toy_checks = packed_module.tiny_tests(core, CoverVerifier, adapter.exact_point, deadline)
    count_checks = counts.checks()
    setup = CoverVerifier(32, deadline)
    records, route_checks = core.route_records(32)
    setup.guard()
    shapes = sorted({key[2] for key in records})
    caps = setup.local_caps(shapes)
    exponent = setup.exact(core.ALPHA)
    powered = {occ: tuple(setup.positive(x)**exponent for x in cap)
               for occ, cap in caps.items()}
    del caps
    progress("local_caps_ready", records=len(records), occupancy_shapes=len(shapes))
    rr, tt = {}, {}
    band_reports = []
    intersections = 0
    for low, high, uscale, vscale in BANDS:
        setup.guard()
        core.USCALE, core.VSCALE = uscale, vscale
        v = CoverVerifier(2*high, deadline); v.packed_full = 128-low
        regular, terminal = v.regions(*local_polynomials(v, core, records, powered))
        for a in range(low, high+1):
            rr[a] = v.extract(regular, 128-a, 2*a)
            tt[a] = v.extract(terminal, 128-a, 2*a)
        intersections += v.intersections
        band_reports.append({"a_min": low, "a_max": high,
            "full_degree": v.packed_full, "singleton_degree": v.single,
            "full_scale": str(uscale), "singleton_scale": str(vscale)})
        progress("band_ready", a_min=low, a_max=high, retained_profiles=len(rr))
        del regular, terminal
    assert set(rr) == set(range(129)) and set(tt) == set(range(129))
    source = (setup.exact(core.ETA)*setup.exact(core.BETA)**256
              *setup.exact(core.Z)**(-core.LIMIT))**exponent
    logtwo = arb(2).log()
    total_exact = arb(0); total_upper = arb(0)
    source_count_sum = 0; chain_count = 0
    reports = []; worst = None; worst_c = None; balanced = []
    for runs in range(1, 257):
        maximum = arb(0)
        weighted = arb(0); profile_total = 0; profile_numbers = {}
        balanced_chain = None
        for a in range(max(0, runs-128), min(128, runs)+1):
            b = runs-a
            value = packed_module.ordered_chain(setup, rr[a], rr[b], tt[b])
            maximum = maximum.max(value)
            # The exact source count is symmetric; the finite conditional
            # chains are not, so both ordered chains are still evaluated.
            if a <= b:
                number = counts.profile_count(4096, 256, 128-a, 128-b)
                profile_numbers[a] = number
            else:
                number = profile_numbers[b]
            weighted += number*value
            profile_total += number
            if a == b:
                balanced_chain = value
            chain_count += 1
        number = counts.run_total(4096, 256, runs)
        assert profile_total == number
        exact_value, upper_value = source*weighted, source*number*maximum
        total_exact += exact_value; total_upper += upper_value
        value = exact_value if args.exact_counts else upper_value
        source_count_sum += number
        if worst is None or value.upper() > worst.upper():
            worst, worst_c = value, runs
        if runs in REPORT_RUNS:
            reports.append({"runs": runs, "indexed_sources": number,
                            "exact_log2_bound_interval": str(exact_value.log()/logtwo),
                            "conservative_log2_bound_interval": str(upper_value.log()/logtwo)})
        if runs in (32, 64, 128, 256):
            balanced.append({"runs": runs, "source_count_kind": "all c-run supports, upper for balanced subset",
                             "log2_bound_interval": str((source*number*balanced_chain).log()/logtwo)})
        if runs in (64, 128, 192, 256):
            progress("chains_ready", through_runs=runs, ordered_chains=chain_count)
        setup.guard()
    assert chain_count == 129**2-1
    assert source_count_sum == comb(8192, 256)
    total = total_exact if args.exact_counts else total_upper
    passed = bool(total < arb(2)**(-40))
    print(json.dumps({"scope": "all state paths and ALL size256 source supports",
        "directed": True, "whole_code_certificate": False,
        "other_source_weights_covered": False,
        "selected_count_bound": "exact profile counts" if args.exact_counts else "total per c times maximum chain",
        "parameters": {"z": str(core.Z), "alpha": str(core.ALPHA),
            "source_weight": 256, "distance_threshold": core.LIMIT, "bands": band_reports},
        "precision_bits": args.precision, "reported_run_bounds": reports,
        "worst_run": {"runs": worst_c, "selection": "largest upper endpoint",
                      "log2_bound_interval": str(worst.log()/logtwo)},
        "balanced_diagnostics": balanced, "log2_bound_interval": str(total.log()/logtwo),
        "exact_log2_bound_interval": str(total_exact.log()/logtwo),
        "conservative_log2_bound_interval": str(total_upper.log()/logtwo),
        "proved_below_2_pow_minus_40": passed,
        "proved_below_2_pow_minus_100": bool(total < arb(2)**(-100)),
        "tests": {**count_checks, "exact_packed_unpacked_coefficients": toy_checks,
            "exact_route_normalizers": route_checks, "all_marginal_profiles": 129,
            "ordered_profile_chains": chain_count, "profile_to_run_totals": 256,
            "source_counts_sum": "C(8192,256)",
            "terminal_updates": 32767},
        "nonnegative_intersections": intersections+setup.intersections,
        "seconds": time.monotonic()-started, "estimated_peak_bytes": 384*(1 << 20),
        "generated_bytes": 0, "core_pin": CORE_PIN, "test_pin": TEST_PIN,
        "packed_pin": PACK_PIN, "count_pin": COUNT_PIN,
        "source_sha256": hashlib.sha256(Path(__file__).read_bytes()).hexdigest()}, indent=2))
    if not passed:
        raise SystemExit("directed interval does not prove the fixed40-bit target")


if __name__ == "__main__":
    main()
