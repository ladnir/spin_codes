"""Directed two-axis residue diagnostic for one actual r512 profile.

Local route terms retain their full-pair/singleton degrees (f,s), and the
local cap is powered before route averaging. Positive cyclic multiplication
modulo (x^mF-1,y^mS-1) retains the desired coefficient plus positive aliases.
The two residues are extracted only after128 route packets. Terminal uses
P^127*T, with no final refresh. Original Cauchy denominators and exact indexed
source counts are unchanged. This is NOT an all-profile certificate.

Packed multiplication uses stride2*mF-1, not mF: the product full degree is
at most2*mF-2, so it cannot carry into the singleton axis. Decode and fold
each axis separately. Default fixed witnesses are the calibrated f223/s66
ones; no marker optimization occurs.60s/256MiB/20GiB; stdout only.
"""

import argparse
from collections import defaultdict
from fractions import Fraction as F
import hashlib
from itertools import combinations_with_replacement
import json
from math import comb
from pathlib import Path
import shutil
import sys
import time

sys.dont_write_bytecode = True
ANCHOR_PIN = "319d169a691abe062e56516b96d8e1800b54040375dbb5fc7981984cbb02eba9"
COUNT_PIN = "4ce61f179a8bf41e5d769e955a9efeb437e46af81f57b0d784c4cfc57958d8e8"
MODULAR_PIN = "128ac26691d4e046267e277516c07594094ae7130e8777f3b1d35f99e4e1f014"


def load(name, pin, label):
    import importlib.util
    path = Path(__file__).with_name(name)
    if hashlib.sha256(path.read_bytes()).hexdigest() != pin:
        raise RuntimeError("frozen dependency hash mismatch: "+name)
    spec = importlib.util.spec_from_file_location(label, path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def multiply(v, left, right, mf, ms):
    """Positive packed product, followed by separate two-axis folding."""
    stride = 2*mf-1
    out = []
    for row in left:
        result = []
        for j in range(len(right[0])):
            v.guard()
            raw = v.poly([])
            for k, item in enumerate(row):
                raw += item*right[k][j]
            assert len(raw) <= (2*mf-2)+stride*(2*ms-2)+1
            coefficients = [v.arb(0) for _ in range((ms-1)*stride+mf)]
            for s in range(ms):
                for f in range(mf):
                    value = v.arb(0)
                    for ss in (s, s+ms):
                        if ss > 2*ms-2:
                            continue
                        for ff in (f, f+mf):
                            # ff==stride is NOT a full-axis alias: it would
                            # read the following singleton row accidentally.
                            if ff <= 2*mf-2:
                                value += raw[ff+stride*ss]
                    coefficients[f+stride*s] = v.positive(value)
            result.append(v.poly(coefficients))
        out.append(result)
    return out


def regions(v, regular, terminal, mf, ms):
    for _ in range(7):
        terminal = multiply(v, regular, terminal, mf, ms)
        regular = multiply(v, regular, regular, mf, ms)
    return regular, terminal


def fold(v, degree_entries, mf, ms):
    stride = 2*mf-1
    entries = [[v.arb(0) for _ in range((ms-1)*stride+mf)] for _ in range(6)]
    for (f, s), values in degree_entries.items():
        index = f % mf+stride*(s % ms)
        for j, value in enumerate(values):
            entries[j][index] += value
    polynomials = [v.poly(row) for row in entries]
    return [[polynomials[0], polynomials[1]], [polynomials[2], polynomials[3]]], [
        [polynomials[4]], [polynomials[5]]]


def local_degrees(v, caps, alpha, u, w):
    """Exact dyadic orientation weights; at most33x33x6 Arb coefficients."""
    half = []
    for f in range(17):
        for s in range(17-f):
            for left in range(s+1):
                a, b = f+left, f+s-left
                if a > b:
                    continue
                numerator = (comb(16, f)*comb(16-f, s)*comb(s, left)
                             *(1 if a == b else 2)*(1 << (16-s)))
                half.append((f, s, a, b, numerator))
    powered = {occ: tuple(v.positive(x)**v.exact(alpha) for x in row)
               for occ, row in caps.items()}
    totals = defaultdict(int)
    values = {}
    for i, (f, s, a, b, mass) in enumerate(half):
        v.guard()
        for j in range(i, len(half)):
            ff, ss, c, d, other = half[j]
            degree = (f+ff, s+ss)
            numerator = mass*other*(1 if i == j else 2)
            totals[degree] += numerator
            if degree not in values:
                values[degree] = [v.arb(0) for _ in range(6)]
            weight = v.exact(F(numerator, 1 << 32))
            row = powered[tuple(sorted((a, b, c, d)))]
            for k in range(6):
                values[degree][k] += weight*row[k]
    normalizers = 0
    for f in range(33):
        for s in range(33-f):
            assert totals[f, s] == comb(32, f)*comb(32-f, s)*(1 << 32)
            normalizers += 1
    up = [v.exact(u)**k for k in range(33)]
    wp = [v.exact(w)**k for k in range(33)]
    for (f, s), row in values.items():
        mass = up[f]*wp[s]
        values[f, s] = [mass*x for x in row]
    return values, normalizers


def exact_tests(v, adapter):
    """Independent rational unwrapped matrix products and carry boundaries."""
    p = [[{(0, 0): F(1, 3), (2, 1): F(2, 5)}, {(1, 2): F(3, 7)}],
         [{(3, 0): F(5, 11)}, {(0, 1): F(7, 13), (1, 3): F(11, 17)}]]
    t = [[{(4, 2): F(13, 19)}], [{(0, 3): F(17, 23), (2, 1): F(19, 29)}]]

    def exact_product(a, b):
        result = []
        for row in a:
            output = []
            for j in range(len(b[0])):
                value = defaultdict(F)
                for k, left in enumerate(row):
                    for (f, s), mass in left.items():
                        for (ff, ss), other in b[k][j].items():
                            value[f+ff, s+ss] += mass*other
                output.append(dict(value))
            result.append(output)
        return result

    expected = exact_product(exact_product(exact_product(p, p), p), t)
    checks = 0
    for mf, ms in ((1, 1), (2, 3), (4, 2)):
        stride = 2*mf-1

        def packed(table):
            out = []
            for row in table:
                rr = []
                for item in row:
                    coefficients = [F(0)]*((ms-1)*stride+mf)
                    for (f, s), mass in item.items():
                        coefficients[f % mf+stride*(s % ms)] += mass
                    rr.append(v.poly([v.exact(x) for x in coefficients]))
                out.append(rr)
            return out

        pp, tt = packed(p), packed(t)
        actual = multiply(v, multiply(v, multiply(v, pp, pp, mf, ms), pp, mf, ms), tt, mf, ms)
        for i in range(2):
            for s in range(ms):
                for f in range(mf):
                    value = sum((mass for (ff, ss), mass in expected[i][0].items()
                                 if ff % mf == f and ss % ms == s), F(0))
                    adapter.exact_point(actual[i][0][f+stride*s], value, v.fmpq)
                    checks += 1
        # Both overflow axes, as well as full-only carry, must fold to(0,0).
        for ss in (0, ms-1):
            a = packed([[{(mf-1, ss): F(1)}]])
            b = packed([[{(1, 0 if ss == 0 else 1): F(1)}]])
            result = multiply(v, a, b, mf, ms)[0][0]
            for s in range(ms):
                for f in range(mf):
                    adapter.exact_point(result[f+stride*s], F(int(f == s == 0)), v.fmpq)
                    checks += 1
    return checks


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--full", type=int, default=223)
    parser.add_argument("--u", default="29411/311125")
    parser.add_argument("--v", default="10441/460035")
    parser.add_argument("--z", default="7/10")
    parser.add_argument("--alpha", default="119/1280")
    parser.add_argument("--max-modulus", type=int, choices=(8, 16, 32), default=32)
    parser.add_argument("--precision", type=int, choices=(256, 384, 512), default=256)
    args = parser.parse_args()
    if not 0 <= args.full <= 255:
        parser.error("full must be0..255; repeated f256 has zero actual count")
    f, s = args.full, 512-2*args.full
    u, w, z, alpha = map(F, (args.u, args.v, args.z, args.alpha))
    if min(u, w) <= 0 or not 0 < z < 1 or not 0 < alpha <= 1:
        parser.error("require positive markers,0<z<1,0<alpha<=1")
    if max(n.bit_length() for x in (u, w, z, alpha)
           for n in (x.numerator, x.denominator)) > 64:
        parser.error("rational numerators/denominators must fit64bits")
    if shutil.disk_usage(Path(__file__).resolve().anchor).free < 20*(1 << 30):
        raise RuntimeError("requires20GiB free on source volume")
    started = time.monotonic(); deadline = started+60
    anchor = load("cheap_lane_cone_cauchy_anchors.py", ANCHOR_PIN, "bimod_anchor")
    counts = load("cheap_lane_fast_profile_counts.py", COUNT_PIN, "bimod_counts")
    modular = load("cheap_lane_cauchy_modular.py", MODULAR_PIN, "bimod_oneaxis")
    counts.DEADLINE = deadline
    peak_memory = anchor.peak_memory_reader()

    def guard():
        if time.monotonic() > deadline:
            raise TimeoutError("60-second bivariate diagnostic budget exceeded")
        if peak_memory() > 256*(1 << 20):
            raise MemoryError("256MiB bivariate diagnostic budget exceeded")

    adapter = anchor.load_adapter()
    core = adapter.load_verifier()
    core.Z, core.ALPHA = z, alpha
    from flint import arb, ctx
    ctx.prec = args.precision

    class BivariateVerifier(core.Verifier):
        def guard(self):
            guard()

    verifier = BivariateVerifier(2, deadline)
    tests = exact_tests(verifier, adapter)
    count_checks = counts.checks()
    caps = verifier.local_caps(list(combinations_with_replacement(range(17), 4)))
    degrees, normalizers = local_degrees(verifier, caps, alpha, u, w)
    del caps
    plain = [sum((row[k] for row in degrees.values()), arb(0)) for k in range(6)]
    plain_p = [[plain[0], plain[1]], [plain[2], plain[3]]]
    plain_t = [[plain[4]], [plain[5]]]
    plain_r = verifier.matpow(plain_p, 128)
    plain_v = verifier.matmul(verifier.matpow(plain_p, 127), plain_t)
    denominator = verifier.exact(u)**f*verifier.exact(w)**s*comb(4096, f)*comb(4096-f, s)
    h, c_h, eta = anchor.source_factor(512)
    exponent = verifier.exact(alpha)
    source = (verifier.exact(eta)*verifier.exact(core.BETA)**512
              *verifier.exact(z)**(-core.LIMIT))**exponent
    number = counts.profile_count(4096, 512, f, f)
    assert number > 0 and number <= anchor.run_total(4096, 512, s)
    logtwo = arb(2).log()

    def bound_from(r, terminal):
        rr = [[x/denominator for x in row] for row in r]
        tt = [[x/denominator for x in row] for row in terminal]
        return number*source*verifier.matmul(verifier.matpow(rr, 255), tt)[0][0]

    ordinary = bound_from(plain_r, plain_v)
    reports = []; residue_checks = 0; collapse_checks = 0
    print(json.dumps({"phase": "local_coefficients_complete", "seconds": time.monotonic()-started,
                      "ordinary_log2_interval": str(ordinary.log()/logtwo)}), flush=True)
    for m in (x for x in (1, 8, 16, 32) if x <= args.max_modulus):
        p, t = fold(verifier, degrees, m, m)
        rr, tt = regions(verifier, p, t, m, m)
        stride = 2*m-1
        for cyclic, evaluated in ((rr, plain_r), (tt, plain_v)):
            for i in range(2):
                for j in range(len(cyclic[0])):
                    total = sum((cyclic[i][j][ff+stride*ss]
                                 for ss in range(m) for ff in range(m)), arb(0))
                    assert total.overlaps(evaluated[i][j])
                    residue_checks += 1
        # Independent frozen one-axis arithmetic checks every collapsed residue.
        one = [[arb(0) for _ in range(m)] for _ in range(6)]
        for (ff, ss), row in degrees.items():
            for k in range(6):
                one[k][(2*ff+ss) % m] += row[k]
        entries = [verifier.poly(row) for row in one]
        op = [[entries[0], entries[1]], [entries[2], entries[3]]]
        ot = [[entries[4]], [entries[5]]]
        orr, ott = modular.cyclic_regions(verifier, op, ot, m)
        for two, one_axis in ((rr, orr), (tt, ott)):
            for i in range(2):
                for j in range(len(two[0])):
                    collapsed = [arb(0) for _ in range(m)]
                    for ss in range(m):
                        for ff in range(m):
                            collapsed[(2*ff+ss) % m] += two[i][j][ff+stride*ss]
                    for k in range(m):
                        assert collapsed[k].overlaps(one_axis[i][j][k])
                        collapse_checks += 1
        index = f % m+stride*(s % m)
        r = [[item[index] for item in row] for row in rr]
        terminal = [[item[index] for item in row] for row in tt]
        bound = bound_from(r, terminal)
        if m == 1:
            assert bound.overlaps(ordinary)
        guard()
        report = {"moduli": [m, m], "target_residues": [f % m, s % m],
                  "log2_bound_interval": str(bound.log()/logtwo),
                  "gain_over_ordinary_bits_interval": str((ordinary/bound).log()/logtwo),
                  "selected_profile_below_2_pow_minus_40": bool(bound < arb(2)**(-40)),
                  "seconds": time.monotonic()-started}
        reports.append(report)
        print(json.dumps({"phase": "moduli_complete", **report}), flush=True)
    guard()
    print(json.dumps({"scope": "one ordered actual r512 profile, equal in both matchings",
        "directed": True, "all_profiles_at_weight_covered": False, "whole_code_certificate": False,
        "parameters": {"source_weight": 512, "ordered_profile": [[f, s], [f, s]],
            "z": str(z), "alpha": str(alpha), "u": str(u), "v": str(w),
            "distance_threshold": core.LIMIT, "fresh_marker_calibration": False},
        "indexed_sources_exact": number, "joint_missing": {"h": h, "c_h_exact": str(c_h),
            "eta_exact": str(eta)}, "ordinary_cauchy_log2_interval": str(ordinary.log()/logtwo),
        "bivariate_results": reports, "tests": {**count_checks,
            "exact_packed_toy_coefficients": tests, "exact_route_degree_normalizers": normalizers,
            "sum_residues_replays_evaluated_entries": residue_checks,
            "frozen_total_degree_collapse_checks": collapse_checks,
            "modulus1_replays_ordinary": True, "terminal_updates": 32767,
            "terminal_route_packets": 128}, "precision_bits": args.precision,
        "seconds": time.monotonic()-started, "peak_process_bytes": peak_memory(),
        "nonnegative_intersections": verifier.intersections, "generated_bytes": 0,
        "anchor_pin": ANCHOR_PIN, "counts_pin": COUNT_PIN, "modular_pin": MODULAR_PIN,
        "core_pin": anchor.CORE_PIN,
        "source_sha256": hashlib.sha256(Path(__file__).read_bytes()).hexdigest()}, indent=2))


if __name__ == "__main__":
    main()
