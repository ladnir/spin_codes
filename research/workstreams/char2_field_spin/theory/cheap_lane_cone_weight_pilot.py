"""One fixed-witness directed pilot at nine predeclared source weights.

Reuse one four-band polynomial pass to extract EVERY marginal profile
(f,s) with2f+s<=256. There are16,641 such profiles, each storing only its
2x2 regular matrix and2x1 terminal. No polynomial is recomputed per source
weight. Local caps, z=21/25 and alpha=119/1280 remain unchanged.

At weight r<=256 use eta_r=256^2/[255^2(257-r)] and beta^r. This is the
exact-support scalar comparison, not a change to the construction. Each
ordered matching profile has s_i=r-2f_i; odd weights are allowed. Sum the
exact indexed profile counts and check C(8192,r). The finite alternating
chain always has32767 updates and one output-only terminal.

Only weights7,8,16,32,64,128,192,255,256 are tested. No interpolation or
whole-code claim. A positive bound is a failed witness, not a verifier
error.180-second cooperative deadline,20GiB free, estimated512MiB peak,
bounded stdout progress, no generated files or frozen-source edits.
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
COVER_PIN = "84ac66e9e3b49b854ce10adf2b2f244d74a530524252736d35d1e72cdfe6e7fb"
WEIGHTS = (7, 8, 16, 32, 64, 128, 192, 255, 256)
CONTROL_CENTER = "-93.5574055068572906550312346546009957059259358406723846743061455591754432859448382230484510871964770748120362074"


def load_cover():
    path = Path(__file__).with_name("cheap_lane_cone_profile_cover.py")
    if hashlib.sha256(path.read_bytes()).hexdigest() != COVER_PIN:
        raise RuntimeError("frozen profile-cover hash mismatch")
    spec = importlib.util.spec_from_file_location("weight_pilot_cover", path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--precision", type=int, choices=(256, 384, 512), default=256)
    args = parser.parse_args()
    if shutil.disk_usage(Path(__file__).resolve().anchor).free < 20*(1 << 30):
        raise RuntimeError("requires20GiB free on source volume")
    started = time.monotonic(); deadline = started+180

    def progress(phase, **values):
        print(json.dumps({"phase": phase, "seconds": time.monotonic()-started, **values}), flush=True)

    cover = load_cover()
    adapter = cover.load("cheap_lane_cone_interval_tests.py", cover.TEST_PIN, "weight_pilot_tests")
    assert adapter.PIN == cover.CORE_PIN
    core = adapter.load_verifier()
    packed = cover.load("cheap_lane_cone_packed.py", cover.PACK_PIN, "weight_pilot_packed")
    assert packed.CORE_PIN == cover.CORE_PIN
    counts = cover.load("cheap_lane_fast_profile_counts.py", cover.COUNT_PIN, "weight_pilot_counts")
    counts.DEADLINE = deadline
    from flint import arb, ctx
    ctx.prec = args.precision

    class PilotVerifier(packed.packed_class(core)):
        def guard(self):
            if time.monotonic() > self.deadline:
                raise TimeoutError("180-second weight-pilot budget exceeded")

    toy_checks = packed.tiny_tests(core, PilotVerifier, adapter.exact_point, deadline)
    count_checks = counts.checks()
    setup = PilotVerifier(32, deadline)
    records, route_checks = core.route_records(32)
    shapes = sorted({key[2] for key in records})
    caps = setup.local_caps(shapes)
    exponent = setup.exact(core.ALPHA)
    powered = {occ: tuple(setup.positive(x)**exponent for x in cap)
               for occ, cap in caps.items()}
    del caps
    progress("local_caps_ready", records=len(records), occupancy_shapes=len(shapes))
    rr, tt = {}, {}
    intersections = 0
    band_reports = []
    for low, high, uscale, vscale in cover.BANDS:
        core.USCALE, core.VSCALE = uscale, vscale
        v = PilotVerifier(2*high, deadline); v.packed_full = 128-low
        regular, terminal = v.regions(*cover.local_polynomials(v, core, records, powered))
        # Assign marginal profiles by full-degree band, not by source weight.
        for f in range(128-high, 128-low+1):
            for s in range(257-2*f):
                assert (f, s) not in rr
                rr[f, s] = v.extract(regular, f, s)
                tt[f, s] = v.extract(terminal, f, s)
        intersections += v.intersections
        band_reports.append({"full_min": 128-high, "full_max": 128-low,
            "singleton_max": v.single, "full_scale": str(uscale), "singleton_scale": str(vscale)})
        progress("band_ready", full_min=128-high, full_max=128-low,
                 retained_marginal_profiles=len(rr))
        del regular, terminal
        setup.guard()
    assert len(rr) == len(tt) == 16641
    assert all((f, s) in rr for f in range(129) for s in range(257-2*f))
    logtwo = arb(2).log()
    chernoff = setup.exact(core.Z)**(-core.LIMIT)
    rows = []; union = arb(0); total_chains = 0
    for weight in WEIGHTS:
        eta = F(256**2, 255**2*(257-weight))
        assert 0 < eta <= core.ETA
        source = (setup.exact(eta)*setup.exact(core.BETA)**weight*chernoff)**exponent
        conditional = arb(0); source_total = 0; chains = 0
        numbers = {}
        by_runs = {}
        maximum = weight//2
        for f0 in range(maximum+1):
            for f1 in range(maximum+1):
                if f0 <= f1:
                    number = counts.profile_count(4096, weight, f0, f1)
                    numbers[f0, f1] = number
                else:
                    number = numbers[f1, f0]
                source_total += number
                runs = weight-f0-f1
                by_runs[runs] = by_runs.get(runs, 0)+number
                if number:
                    s0, s1 = weight-2*f0, weight-2*f1
                    value = packed.ordered_chain(setup, rr[f0, s0], rr[f1, s1], tt[f1, s1])
                    conditional += number*value
                    chains += 1
            setup.guard()
        assert source_total == comb(8192, weight)
        assert all(by_runs.get(c, 0) == counts.run_total(4096, weight, c)
                   for c in range(weight+1))
        value = source*conditional
        union += value; total_chains += chains
        bits = value.log()/logtwo
        if weight == 256:
            assert eta == core.ETA
            assert bits.overlaps(arb(CONTROL_CENTER, "3.25e-110"))
        row = {"source_weight": weight, "eta_exact": str(eta),
               "indexed_sources": source_total, "ordered_profile_chains": chains,
               "log2_bound_interval": str(bits),
               "proved_below_2_pow_minus_40": bool(value < arb(2)**(-40)),
               "proved_below_2_pow_minus_100": bool(value < arb(2)**(-100))}
        rows.append(row)
        progress("weight_ready", source_weight=weight, log2_bound_interval=str(bits),
                 proved_below_2_pow_minus_40=row["proved_below_2_pow_minus_40"])
    setup.guard()
    print(json.dumps({"scope": "all state paths and all supports at nine predeclared source weights",
        "directed": True, "verifier_completed": True, "whole_code_certificate": False,
        "unlisted_weights_covered": False,
        "parameters": {"weights": WEIGHTS, "z": str(core.Z), "alpha": str(core.ALPHA),
            "distance_threshold": core.LIMIT, "bands": band_reports},
        "precision_bits": args.precision, "weight_bounds": rows,
        "selected_weight_union_log2_interval": str(union.log()/logtwo),
        "selected_weight_union_below_2_pow_minus_40": bool(union < arb(2)**(-40)),
        "all_individual_weights_below_2_pow_minus_40": all(row["proved_below_2_pow_minus_40"] for row in rows),
        "tests": {**count_checks, "exact_packed_unpacked_coefficients": toy_checks,
            "exact_route_normalizers": route_checks, "all_triangle_marginals": 16641,
            "all_weight_source_totals": len(WEIGHTS), "per_run_totals": "passed",
            "r256_control": "overlaps frozen384-bit interval", "terminal_updates": 32767},
        "total_ordered_profile_chains": total_chains,
        "nonnegative_intersections": intersections+setup.intersections,
        "seconds": time.monotonic()-started, "estimated_peak_bytes": 512*(1 << 20),
        "generated_bytes": 0, "cover_pin": COVER_PIN, "core_pin": cover.CORE_PIN,
        "test_pin": cover.TEST_PIN, "packed_pin": cover.PACK_PIN, "count_pin": cover.COUNT_PIN,
        "source_sha256": hashlib.sha256(Path(__file__).read_bytes()).hexdigest()}, indent=2))
    # A completed proof calculation may return a positive bound. That is a
    # failed fixed witness, not an arithmetic or source-verification error.


if __name__ == "__main__":
    main()
