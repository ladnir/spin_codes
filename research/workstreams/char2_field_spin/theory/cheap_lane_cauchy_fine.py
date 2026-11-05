"""Finer directed coefficient diagnostic for the fixed r512/f223/s66 profile.

Default mode uses the frozen positive two-axis residue64 calculation.
Exact mode retains the downward-closed rectangle f<=223,s<=66 after every
matrix product and extracts its corner from P^128 and P^127*T. Packed stride
447 prevents full-axis carries before truncation. Discarded coefficients
have nonnegative degrees and cannot contribute to the retained rectangle.

Both modes reuse frozen local caps, cap^alpha before route averaging, exact
indexed source count, and the output-only terminal. The exact coefficient
is independent of positive route markers after normalization; residues
retain positive aliases and need not be marker-invariant. This covers ONE
ordered matching profile, not all source supports. No scans or calibration.
120s/512MiB/20GiB, stdout only, no generated files or bytecode.
"""

import argparse
from collections import defaultdict
from fractions import Fraction as F
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
FULL, SINGLE, WEIGHT = 223, 66, 512


def load_base():
    path = Path(__file__).with_name("cheap_lane_cauchy_bivariate_modular.py")
    if hashlib.sha256(path.read_bytes()).hexdigest() != BASE_PIN:
        raise RuntimeError("frozen bivariate source hash mismatch")
    spec = importlib.util.spec_from_file_location("fine_bivariate", path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def rectangle_product(v, left, right, full, single):
    """Carry-safe ordinary product, then positive downset truncation."""
    stride = 2*full+1
    out = []
    for row in left:
        result = []
        for j in range(len(right[0])):
            v.guard()
            raw = v.poly([])
            for k, value in enumerate(row):
                raw += value*right[k][j]
            assert len(raw) <= 2*full+stride*(2*single)+1
            coefficients = [v.arb(0) for _ in range(single*stride+full+1)]
            for s in range(single+1):
                for f in range(full+1):
                    coefficients[f+stride*s] = v.positive(raw[f+stride*s])
            result.append(v.poly(coefficients))
        out.append(result)
    return out


def rectangle_pack(v, degrees, full, single):
    stride = 2*full+1
    entries = [[v.arb(0) for _ in range(single*stride+full+1)] for _ in range(6)]
    for (f, s), values in degrees.items():
        if f <= full and s <= single:
            for k, value in enumerate(values):
                entries[k][f+stride*s] += value
    p = [v.poly(row) for row in entries]
    return [[p[0], p[1]], [p[2], p[3]]], [[p[4]], [p[5]]]


def exact_tests(v, adapter):
    """Unwrapped rational products test truncation and carry boundaries."""
    p = [[{(0, 0): F(1, 3), (2, 0): F(2, 5)}, {(1, 2): F(3, 7)}],
         [{(3, 1): F(5, 11)}, {(0, 1): F(7, 13), (1, 3): F(11, 17)}]]
    t = [[{(4, 0): F(13, 19)}], [{(0, 2): F(17, 23), (2, 1): F(19, 29)}]]

    def exact_product(a, b):
        out = []
        for row in a:
            result = []
            for j in range(len(b[0])):
                item = defaultdict(F)
                for k, left in enumerate(row):
                    for (f, s), mass in left.items():
                        for (ff, ss), other in b[k][j].items():
                            item[f+ff, s+ss] += mass*other
                result.append(dict(item))
            out.append(result)
        return out

    expected = exact_product(exact_product(exact_product(p, p), p), t)
    checks = 0
    for full, single in ((0, 3), (2, 3), (4, 2)):
        stride = 2*full+1

        def pack(table):
            result = []
            for row in table:
                rr = []
                for item in row:
                    coefficients = [F(0)]*(single*stride+full+1)
                    for (f, s), mass in item.items():
                        if f <= full and s <= single:
                            coefficients[f+stride*s] += mass
                    rr.append(v.poly([v.exact(x) for x in coefficients]))
                result.append(rr)
            return result

        pp, tt = pack(p), pack(t)
        actual = rectangle_product(v, rectangle_product(v,
            rectangle_product(v, pp, pp, full, single), pp, full, single), tt, full, single)
        for i in range(2):
            for s in range(single+1):
                for f in range(full+1):
                    adapter.exact_point(actual[i][0][f+stride*s],
                                        expected[i][0].get((f, s), F(0)), v.fmpq)
                    checks += 1
        if full:
            # Overflow full degree must be discarded, never read as y^1.
            a, b = pack([[{(full, 0): F(1)}]]), pack([[{(1, 0): F(1)}]])
            result = rectangle_product(v, a, b, full, single)[0][0]
            for s in range(single+1):
                for f in range(full+1):
                    adapter.exact_point(result[f+stride*s], F(0), v.fmpq)
                    checks += 1
    return checks


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--mode", choices=("residue64", "exact"), default="residue64")
    parser.add_argument("--z", default="17447/25000")
    parser.add_argument("--alpha", default="119/1280")
    parser.add_argument("--u", default="29411/311125")
    parser.add_argument("--v", default="10441/460035")
    parser.add_argument("--precision", type=int, choices=(256, 384, 512), default=256)
    parser.add_argument("--control32", action="store_true")
    args = parser.parse_args()
    z, alpha, u, w = map(F, (args.z, args.alpha, args.u, args.v))
    if not 0 < z < 1 or not 0 < alpha <= 1 or min(u, w) <= 0:
        parser.error("require0<z<1,0<alpha<=1 and positive markers")
    if max(n.bit_length() for x in (z, alpha, u, w)
           for n in (x.numerator, x.denominator)) > 64:
        parser.error("rational numerator/denominator limit is64bits")
    if shutil.disk_usage(Path(__file__).resolve().anchor).free < 20*(1 << 30):
        raise RuntimeError("requires20GiB free on source volume")
    started = time.monotonic(); deadline = started+120
    base = load_base()
    anchor = base.load("cheap_lane_cone_cauchy_anchors.py", base.ANCHOR_PIN, "fine_anchor")
    counts = base.load("cheap_lane_fast_profile_counts.py", base.COUNT_PIN, "fine_counts")
    counts.DEADLINE = deadline
    peak_memory = anchor.peak_memory_reader()

    def guard():
        if time.monotonic() > deadline:
            raise TimeoutError("120-second fine coefficient budget exceeded")
        if peak_memory() > 512*(1 << 20):
            raise MemoryError("512MiB fine coefficient budget exceeded")

    adapter = anchor.load_adapter()
    core = adapter.load_verifier()
    core.Z, core.ALPHA = z, alpha
    from flint import arb, ctx
    ctx.prec = args.precision

    class FineVerifier(core.Verifier):
        def guard(self):
            guard()

    verifier = FineVerifier(2, deadline)
    checks = exact_tests(verifier, adapter)
    cyclic_checks = base.exact_tests(verifier, adapter)
    count_checks = counts.checks()
    caps = verifier.local_caps(list(combinations_with_replacement(range(17), 4)))
    degrees, normalizers = base.local_degrees(verifier, caps, alpha, u, w)
    del caps
    plain = [sum((row[k] for row in degrees.values()), arb(0)) for k in range(6)]
    plain_p = [[plain[0], plain[1]], [plain[2], plain[3]]]
    plain_t = [[plain[4]], [plain[5]]]
    plain_r = verifier.matpow(plain_p, 128)
    plain_v = verifier.matmul(verifier.matpow(plain_p, 127), plain_t)
    denominator = (verifier.exact(u)**FULL*verifier.exact(w)**SINGLE
                   *comb(4096, FULL)*comb(4096-FULL, SINGLE))
    h, c_h, eta = anchor.source_factor(WEIGHT)
    source = (verifier.exact(eta)*verifier.exact(core.BETA)**WEIGHT
              *verifier.exact(z)**(-core.LIMIT))**verifier.exact(alpha)
    number = counts.profile_count(4096, WEIGHT, FULL, FULL)
    assert number > 0 and number <= anchor.run_total(4096, WEIGHT, SINGLE)
    logtwo = arb(2).log()

    def bound_from(r, terminal):
        rr = [[x/denominator for x in row] for row in r]
        tt = [[x/denominator for x in row] for row in terminal]
        return number*source*verifier.matmul(verifier.matpow(rr, 255), tt)[0][0]

    ordinary = bound_from(plain_r, plain_v)
    print(json.dumps({"phase": "local_coefficients_complete", "seconds": time.monotonic()-started,
                      "ordinary_log2_interval": str(ordinary.log()/logtwo)}), flush=True)
    controls = []

    def residue_result(m):
        p, t = base.fold(verifier, degrees, m, m)
        rr, tt = base.regions(verifier, p, t, m, m)
        stride = 2*m-1
        for cyclic, evaluated in ((rr, plain_r), (tt, plain_v)):
            for i in range(2):
                for j in range(len(cyclic[0])):
                    total = sum((cyclic[i][j][f+stride*s]
                                 for s in range(m) for f in range(m)), arb(0))
                    assert total.overlaps(evaluated[i][j])
        index = FULL % m+stride*(SINGLE % m)
        return bound_from([[x[index] for x in row] for row in rr],
                          [[x[index] for x in row] for row in tt])

    if args.control32:
        control = residue_result(32)
        controls.append({"moduli": [32, 32], "log2_interval": str(control.log()/logtwo)})
        print(json.dumps({"phase": "control32_complete", **controls[-1],
                          "seconds": time.monotonic()-started}), flush=True)
    if args.mode == "residue64":
        bound = residue_result(64)
        coefficient_scope = "positive aliases modulo64 in both route degrees"
    else:
        regular, terminal = rectangle_pack(verifier, degrees, FULL, SINGLE)
        for step in range(7):
            terminal = rectangle_product(verifier, regular, terminal, FULL, SINGLE)
            regular = rectangle_product(verifier, regular, regular, FULL, SINGLE)
            print(json.dumps({"phase": "exact_doubling_complete", "route_packets": 1 << (step+1),
                              "seconds": time.monotonic()-started}), flush=True)
        index = FULL+(2*FULL+1)*SINGLE
        r = [[x[index] for x in row] for row in regular]
        terminal = [[x[index] for x in row] for row in terminal]
        bound = bound_from(r, terminal)
        coefficient_scope = "exact selected coefficient of the existing local-cap majorant"
    guard()
    print(json.dumps({"scope": "one ordered actual r512 matching profile only",
        "mode": args.mode, "coefficient_scope": coefficient_scope,
        "directed": True, "all_profiles_at_weight_covered": False, "whole_code_certificate": False,
        "parameters": {"ordered_profile": [[FULL, SINGLE], [FULL, SINGLE]], "source_weight": WEIGHT,
            "z": str(z), "alpha": str(alpha), "u": str(u), "v": str(w),
            "distance_threshold": core.LIMIT, "fresh_marker_calibration": False},
        "indexed_sources_exact": number, "joint_missing": {"h": h, "c_h_exact": str(c_h),
            "eta_exact": str(eta)}, "log2_bound_interval": str(bound.log()/logtwo),
        "ordinary_cauchy_log2_interval": str(ordinary.log()/logtwo),
        "gain_over_ordinary_bits_interval": str((ordinary/bound).log()/logtwo),
        "selected_profile_below_2_pow_minus_40": bool(bound < arb(2)**(-40)),
        "controls": controls, "tests": {**count_checks,
            "exact_rectangle_toy_coefficients": checks, "exact_cyclic_toy_coefficients": cyclic_checks,
            "exact_route_degree_normalizers": normalizers, "terminal_updates": 32767,
            "terminal_route_packets": 128}, "precision_bits": args.precision,
        "seconds": time.monotonic()-started, "peak_process_bytes": peak_memory(),
        "nonnegative_intersections": verifier.intersections, "generated_bytes": 0,
        "base_pin": BASE_PIN, "counts_pin": base.COUNT_PIN, "anchor_pin": base.ANCHOR_PIN,
        "core_pin": anchor.CORE_PIN,
        "source_sha256": hashlib.sha256(Path(__file__).read_bytes()).hexdigest()}, indent=2))


if __name__ == "__main__":
    main()
