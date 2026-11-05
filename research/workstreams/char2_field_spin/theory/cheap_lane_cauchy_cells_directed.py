"""Directed marker-bank replay for all actual source supports of weight 512.

At fixed z and alpha, each local cone cap is powered BEFORE route averaging.
One positive residue128 table per marker cell bounds P^128 and P^127*T.
Every assigned marginal (f,s=512-2f) is extracted and divided by its OWN
M_4096(f,s)*u^f*v^s. Cell assignment selects a valid witness; it does not
interpolate coefficients or change the marker inside a state path.

All 66048 ordered matching profiles are counted exactly. For R0,R1 and
output-only terminals V0,V1, share K=(R0*R1)^127 between the two orientations:
e0*K*R0*V1 and e0*R1*K*V0. Each chain has 32767 updates and one final emission.
The joint-missing source factor is applied once, and source counts are outside
alpha. This covers one weight, NOT the whole code. No FLOAT bound is used.

Only 257 normalized 2x2 matrices and terminals survive between marker cells.
The cooperative budget is 600 or 900 seconds, peak process memory 512MiB,
and source-volume free space at least20GiB. Stdout only; no generated files.
"""

import argparse
from fractions import Fraction as F
import gc
import hashlib
import importlib.util
from itertools import combinations_with_replacement
import json
from math import comb
from pathlib import Path
import shutil
import sys
import time

sys.dont_write_bytecode = True
BASE_PIN = "01356ae20b8b00a3354d76c4c7978bcb3d77d73293b110ea2a9cd4dc9b495d95"
SCREEN_PIN = "6cd1f551da748c461c995f7927144345655cc741156d0550625388fcd6876861"
CELLS_PIN = "0ed129b48ee91cfc456c2e434ea5e279bd6ad777358f836d6152d2d2c9c813b1"
DEFAULT_CENTERS = "16,48,80,112,144,176,208,224,240,252"
MODULUS = 128
CURVE_FULL, CURVE_SINGLE = F(1562175, 990029), F(1140210, 880327)


def load_base():
    path = Path(__file__).with_name("cheap_lane_cauchy_bivariate_modular.py")
    if hashlib.sha256(path.read_bytes()).hexdigest() != BASE_PIN:
        raise RuntimeError("frozen bivariate helper hash mismatch")
    spec = importlib.util.spec_from_file_location("directed_cells_base", path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def orientation_tests():
    """Exact noncommuting matrices and distinct terminals test both orders."""
    def mul(a, b):
        return [[sum((x*b[k][j] for k, x in enumerate(row)), F(0))
                 for j in range(len(b[0]))] for row in a]
    a = [[F(1, 3), F(2, 5)], [F(3, 7), F(5, 11)]]
    b = [[F(7, 13), F(11, 17)], [F(13, 19), F(17, 23)]]
    u, v = [[F(19, 29)], [F(23, 31)]], [[F(29, 37)], [F(31, 41)]]
    eye = [[F(1), F(0)], [F(0), F(1)]]
    ab, ba = eye, eye
    checks = 0
    for _ in range(4):
        assert mul(mul(ab, a), v) == mul(a, mul(ba, v))
        assert mul(mul(b, ab), u) == mul(mul(ba, b), u)
        ab, ba = mul(ab, mul(a, b)), mul(ba, mul(b, a))
        checks += 2
    return checks


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--centers", default=DEFAULT_CENTERS)
    parser.add_argument("--precision", type=int, choices=(256, 384), default=256)
    parser.add_argument("--seconds", type=int, choices=(600, 900), default=600)
    parser.add_argument("--z", default="17447/25000")
    parser.add_argument("--alpha", default="119/1280")
    args = parser.parse_args()
    try:
        centers = tuple(int(x) for x in args.centers.split(","))
        z, alpha = F(args.z), F(args.alpha)
    except (ValueError, ZeroDivisionError):
        parser.error("centers must be integers; z and alpha must be rationals")
    if not (1 <= len(centers) <= 16 and tuple(sorted(set(centers))) == centers
            and 0 <= centers[0] <= centers[-1] <= 256):
        parser.error("require 1..16 strictly increasing centers in0..256")
    if not (0 < z < 1 and 0 < alpha <= 1):
        parser.error("require0<z<1 and0<alpha<=1")
    if max(n.bit_length() for x in (z, alpha)
           for n in (x.numerator, x.denominator)) > 64:
        parser.error("rational numerators/denominators must fit64bits")
    if shutil.disk_usage(Path(__file__).resolve().anchor).free < 20*(1 << 30):
        raise RuntimeError("requires20GiB free on source volume")
    started = time.monotonic(); deadline = started+args.seconds
    base = load_base()
    anchor = base.load("cheap_lane_cone_cauchy_anchors.py", base.ANCHOR_PIN, "cell_anchor")
    counts = base.load("cheap_lane_fast_profile_counts.py", base.COUNT_PIN, "cell_counts")
    screen = base.load("cheap_lane_cauchy_profile_screen.py", SCREEN_PIN, "cell_screen")
    cells = base.load("cheap_lane_cauchy_cells.py", CELLS_PIN, "cell_assignment")
    counts.DEADLINE = deadline
    peak_memory = anchor.peak_memory_reader()

    def guard():
        if time.monotonic() > deadline:
            raise TimeoutError("directed marker-bank time budget exceeded")
        if peak_memory() > 512*(1 << 20):
            raise MemoryError("directed marker-bank512MiB budget exceeded")

    adapter = anchor.load_adapter(); core = adapter.load_verifier()
    core.Z, core.ALPHA = z, alpha
    from flint import arb, ctx
    ctx.prec = args.precision

    class CellVerifier(core.Verifier):
        def guard(self):
            guard()

    v = CellVerifier(2, deadline)
    assignment, assigned = cells.assign_centers(centers)
    assert cells.DEFAULT_CENTERS == DEFAULT_CENTERS
    toy_checks = base.exact_tests(v, adapter)
    scale_checks = cells.exact_scale_tests()
    orientation_checks = orientation_tests()
    count_checks = counts.checks()
    assert screen.marker_curve(16, CURVE_FULL, CURVE_SINGLE) == (
        F(6943, 990029), F(152028, 880327))
    caps = v.local_caps(list(combinations_with_replacement(range(17), 4)))
    degrees, normalizers = base.local_degrees(v, caps, alpha, F(1), F(1))
    del caps
    regular, terminal = [None]*257, [None]*257
    witnesses = []
    dc_checks = 0
    logtwo = arb(2).log()
    print(json.dumps({"phase": "local_coefficients_complete", "seconds": time.monotonic()-started,
                      "marker_cells": len(centers)}), flush=True)
    for center in centers:
        guard()
        u, w = screen.marker_curve(center, CURVE_FULL, CURVE_SINGLE)
        up = [v.exact(u)**k for k in range(33)]
        wp = [v.exact(w)**k for k in range(33)]
        scaled = {(f, s): [up[f]*wp[s]*x for x in row]
                  for (f, s), row in degrees.items()}
        plain = [sum((row[k] for row in scaled.values()), arb(0)) for k in range(6)]
        plain_p = [[plain[0], plain[1]], [plain[2], plain[3]]]
        plain_t = [[plain[4]], [plain[5]]]
        plain_r = v.matpow(plain_p, 128)
        plain_v = v.matmul(v.matpow(plain_p, 127), plain_t)
        packet, last = base.fold(v, scaled, MODULUS, MODULUS)
        rr, tt = base.regions(v, packet, last, MODULUS, MODULUS)
        stride = 2*MODULUS-1
        for cyclic, evaluated in ((rr, plain_r), (tt, plain_v)):
            for i in range(2):
                for j in range(len(cyclic[0])):
                    total = sum((cyclic[i][j][f+stride*s]
                                 for s in range(MODULUS) for f in range(MODULUS)), arb(0))
                    assert total.overlaps(evaluated[i][j])
                    dc_checks += 1
        for f in assigned[center]:
            s = 512-2*f
            index = f % MODULUS+stride*(s % MODULUS)
            denom = (arb(comb(4096, f)*comb(4096-f, s))
                     *v.exact(u)**f*v.exact(w)**s)
            assert denom > 0
            regular[f] = [[v.positive(item[index])/denom for item in row] for row in rr]
            terminal[f] = [[v.positive(item[index])/denom for item in row] for row in tt]
        witnesses.append({"center": center, "u": str(u), "v": str(w),
                          "first_full": min(assigned[center]), "last_full": max(assigned[center])})
        del rr, tt, packet, last, scaled, plain, plain_p, plain_t, plain_r, plain_v
        # Release loop aliases too: otherwise one packed table survives the cell.
        del cyclic, evaluated
        gc.collect(); guard()
        print(json.dumps({"phase": "cell_complete", **witnesses[-1],
            "seconds": time.monotonic()-started, "peak_process_bytes": peak_memory()}), flush=True)
    assert all(x is not None for x in regular+terminal)
    h, c_h, eta = anchor.source_factor(512)
    source = (v.exact(eta)*v.exact(core.BETA)**512*v.exact(z)**(-core.LIMIT))**v.exact(alpha)
    run_bounds = [arb(0) for _ in range(513)]
    run_counts = [0]*513
    ordered = 0
    worst = None
    selected = None
    for f0 in range(257):
        guard()
        for f1 in range(f0, 257):
            number = counts.profile_count(4096, 512, f0, f1)
            if not number:
                assert f0 == f1 == 256
                continue
            runs = 512-f0-f1
            r0, r1 = regular[f0], regular[f1]
            shared = v.matpow(v.matmul(r0, r1), 127)
            forward = v.positive(v.matmul(v.matmul(shared, r0), terminal[f1])[0][0])
            terms = [(f0, f1, forward)]
            if f0 != f1:
                reverse = v.positive(v.matmul(v.matmul(r1, shared), terminal[f0])[0][0])
                terms.append((f1, f0, reverse))
            run_counts[runs] += len(terms)*number
            ordered += len(terms)
            for a, b, value in terms:
                mass = number*value
                run_bounds[runs] += mass
                if worst is None or mass.upper() > worst[2].upper():
                    worst = (a, b, mass)
                if a == b == 223:
                    selected = source*mass
        if f0 % 32 == 0:
            print(json.dumps({"phase": "profile_rows_complete", "through_full0": f0,
                "ordered_profiles": ordered, "seconds": time.monotonic()-started}), flush=True)
    assert ordered == 66048 and selected is not None
    for runs in range(1, 513):
        assert run_counts[runs] == counts.run_total(4096, 512, runs)
    total_count = sum(run_counts)
    assert total_count == comb(8192, 512)
    total = source*sum(run_bounds, arb(0))
    assert total > 0
    worst_run = max(range(1, 513), key=lambda c: run_bounds[c].upper())
    parameters = {"source_weight": 512, "distance_threshold": core.LIMIT,
        "z": str(z), "alpha": str(alpha), "moduli": [MODULUS, MODULUS],
        "curve_full": str(CURVE_FULL), "curve_single": str(CURVE_SINGLE),
        "cells": witnesses, "assignment": assignment}
    witness_hash = hashlib.sha256(json.dumps(parameters, sort_keys=True,
                                              separators=(",", ":")).encode()).hexdigest()
    guard()
    print(json.dumps({"scope": "all actual source supports at weight512, all matching profiles and state paths",
        "directed": True, "all_profiles_at_weight_covered": True, "whole_code_certificate": False,
        "log2_union_bound_interval": str(total.log()/logtwo),
        "weight512_below_2_pow_minus_40": bool(total < arb(2)**(-40)),
        "weight512_below_2_pow_minus_50": bool(total < arb(2)**(-50)),
        "selected_f223_log2_bound_interval": str(selected.log()/logtwo),
        "selected_f223_marker_center": assignment[223],
        "worst_ordered_profile": {"full": list(worst[:2]),
            "log2_bound_interval": str((source*worst[2]).log()/logtwo)},
        "worst_run_count": worst_run,
        "sampled_run_unions": [{"runs": c, "log2_bound_interval": str((source*run_bounds[c]).log()/logtwo)}
            for c in sorted({1, 16, 32, 64, 66, 128, 256, 480, 512, worst_run})],
        "parameters": parameters, "joint_missing": {"h": h, "c_h": str(c_h), "eta": str(eta)},
        "tests": {**count_checks, "exact_packed_toy_checks": toy_checks,
            "exact_marker_scaling_checks": scale_checks, "orientation_checks": orientation_checks,
            "exact_local_route_normalizers": normalizers, "dc_overlap_checks": dc_checks,
            "ordered_profiles": ordered, "all512_run_count_sums_exact": True,
            "total_indexed_supports_equals_binomial8192_512": True,
            "terminal_updates": 32767, "route_packets_per_region": 128},
        "precision_bits": args.precision, "budget_seconds": args.seconds,
        "seconds": time.monotonic()-started, "peak_process_bytes": peak_memory(),
        "nonnegative_intersections": v.intersections, "generated_bytes": 0,
        "witness_sha256": witness_hash, "base_pin": BASE_PIN, "screen_pin": SCREEN_PIN,
        "cells_pin": CELLS_PIN, "anchor_pin": base.ANCHOR_PIN, "count_pin": base.COUNT_PIN,
        "adapter_pin": anchor.ADAPTER_PIN, "core_pin": anchor.CORE_PIN,
        "source_sha256": hashlib.sha256(Path(__file__).read_bytes()).hexdigest()}, indent=2))


if __name__ == "__main__":
    main()
