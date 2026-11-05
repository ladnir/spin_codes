"""Directed modular-degree diagnostic for one actual r512 matching profile.

At fixed positive u,v, each local term has total input degree ell=2f+s.
Place its exact route weight*u^f*v^s*cap^alpha in residue ell modulo m.
All coefficients are nonnegative. In the positive cyclic ring t^m=1,
the r512 residue of P^128 or P^127*T contains every desired bivariate
coefficient (f,512-2f), plus aliases. Divide by the unchanged exact
u^f*v^(512-2f)*M_4096(f,512-2f). This refines Cauchy evaluation without
Jensen or a uniform-residue assumption. Terminal route degrees count;
no final update is added. The actual profile is identical in both
matchings, so its finite chain is e0*R^255*V.

Default f223 uses the frozen one-center correction curve. Optional u,v
must be supplied together. Evaluate m=1,4,8, optionally16. These are
diagnostics, not all-profile or whole-code certificates. Exact indexed
profile counts, not the larger complete run count, weight every result.
No coefficient tails are discarded.60s/128MiB/20GiB, stdout only.
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
ANCHOR_PIN = "319d169a691abe062e56516b96d8e1800b54040375dbb5fc7981984cbb02eba9"
COUNT_PIN = "4ce61f179a8bf41e5d769e955a9efeb437e46af81f57b0d784c4cfc57958d8e8"
R = 512


def load(name, pin, label):
    path = Path(__file__).with_name(name)
    if hashlib.sha256(path.read_bytes()).hexdigest() != pin:
        raise RuntimeError("frozen dependency hash mismatch: "+name)
    spec = importlib.util.spec_from_file_location(label, path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def cyclic_product(v, left, right, modulus):
    """Positive wraparound, never a signed remainder or truncation."""
    out = []
    for row in left:
        result = []
        for j in range(len(right[0])):
            v.guard()
            raw = v.poly([])
            for k, item in enumerate(row):
                raw += item*right[k][j]
            # Every operand has degree<m, so the product has degree<=2m-2.
            assert len(raw) <= 2*modulus-1
            result.append(v.poly([v.positive(raw[d]+raw[d+modulus])
                                  for d in range(modulus)]))
        out.append(result)
    return out


def cyclic_regions(v, regular, terminal, modulus):
    for _ in range(7):
        terminal = cyclic_product(v, regular, terminal, modulus)
        regular = cyclic_product(v, regular, regular, modulus)
    return regular, terminal


def exact_tests(v, adapter):
    """Small noncommuting polynomial matrices versus exact unwrapped products."""
    p = [[{0: F(1, 3), 2: F(2, 5)}, {1: F(3, 7)}],
         [{3: F(5, 11)}, {0: F(7, 13), 1: F(11, 17)}]]
    t = [[{4: F(13, 19)}], [{0: F(17, 23), 2: F(19, 29)}]]

    def multiply(a, b):
        result = []
        for row in a:
            output = []
            for j in range(len(b[0])):
                value = defaultdict(F)
                for k, left in enumerate(row):
                    for degree, mass in left.items():
                        for other, weight in b[k][j].items():
                            value[degree+other] += mass*weight
                output.append(dict(value))
            result.append(output)
        return result

    expected = multiply(multiply(multiply(p, p), p), t)
    checks = 0
    for m in (1, 2, 4):
        def fold(table):
            out = []
            for row in table:
                rr = []
                for item in row:
                    coefficients = [F(0)]*m
                    for degree, mass in item.items():
                        coefficients[degree % m] += mass
                    rr.append(v.poly([v.exact(x) for x in coefficients]))
                out.append(rr)
            return out

        pp, tt = fold(p), fold(t)
        actual = cyclic_product(v, cyclic_product(v, cyclic_product(v, pp, pp, m), pp, m), tt, m)
        for i in range(2):
            for residue in range(m):
                value = sum((mass for degree, mass in expected[i][0].items()
                             if degree % m == residue), F(0))
                adapter.exact_point(actual[i][0][residue], value, v.fmpq)
                checks += 1
    return checks


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--full", type=int, default=223)
    parser.add_argument("--u")
    parser.add_argument("--v")
    parser.add_argument("--z", default="7/10")
    parser.add_argument("--alpha", default="119/1280")
    parser.add_argument("--include-16", action="store_true")
    parser.add_argument("--precision", type=int, choices=(256, 384, 512), default=256)
    args = parser.parse_args()
    if not 0 <= args.full <= 255:
        parser.error("full must be0..255; the repeated f256 profile has zero actual source count")
    if (args.u is None) != (args.v is None):
        parser.error("supply both positive rational markers, or neither")
    f, s = args.full, R-2*args.full
    if args.u is None:
        u = F(1562175, 990029)*F(f, 3584+f) if f else F(1, 1 << 30)
        vv = F(1140210, 880327)*F(s, 3584+f)
        marker_source = "one-center correction curve, not a new calibration"
    else:
        u, vv = F(args.u), F(args.v)
        marker_source = "caller-specified fixed rational witnesses"
    z, alpha = F(args.z), F(args.alpha)
    if min(u, vv) <= 0 or not 0 < z < 1 or not 0 < alpha <= 1:
        parser.error("require positive markers,0<z<1,0<alpha<=1")
    if max(n.bit_length() for x in (u, vv, z, alpha)
           for n in (x.numerator, x.denominator)) > 64:
        parser.error("rational parameters must fit64-bit numerators/denominators")
    if shutil.disk_usage(Path(__file__).resolve().anchor).free < 20*(1 << 30):
        raise RuntimeError("requires20GiB free on source volume")
    started = time.monotonic(); deadline = started+60
    anchor = load("cheap_lane_cone_cauchy_anchors.py", ANCHOR_PIN, "modular_anchor")
    counts = load("cheap_lane_fast_profile_counts.py", COUNT_PIN, "modular_counts")
    counts.DEADLINE = deadline
    peak_memory = anchor.peak_memory_reader()

    def guard():
        if time.monotonic() > deadline:
            raise TimeoutError("60-second modular diagnostic budget exceeded")
        if peak_memory() > anchor.MEMORY_LIMIT:
            raise MemoryError("128MiB modular diagnostic budget exceeded")

    adapter = anchor.load_adapter()
    core = adapter.load_verifier()
    core.Z, core.ALPHA = z, alpha
    from flint import arb, ctx
    ctx.prec = args.precision

    class ModularVerifier(core.Verifier):
        def guard(self):
            guard()

    verifier = ModularVerifier(2, deadline)
    checks = exact_tests(verifier, adapter)
    count_checks = counts.checks()
    shapes = list(combinations_with_replacement(range(17), 4))
    weights = anchor.shape_weights(u, vv, guard)
    assert set(weights) == set(shapes)
    caps = verifier.local_caps(shapes)
    exponent = verifier.exact(alpha)
    degrees = [[arb(0) for _ in range(65)] for _ in range(6)]
    plain = [arb(0) for _ in range(6)]
    for occ, weight in weights.items():
        guard()
        mass, degree = verifier.exact(weight), sum(occ)
        for j, cap in enumerate(caps[occ]):
            term = mass*verifier.positive(cap)**exponent
            degrees[j][degree] += term
            plain[j] += term
    del caps, weights
    plain_p = [[plain[0], plain[1]], [plain[2], plain[3]]]
    plain_t = [[plain[4]], [plain[5]]]
    plain_r = verifier.matpow(plain_p, 128)
    plain_v = verifier.matmul(verifier.matpow(plain_p, 127), plain_t)
    denominator = verifier.exact(u)**f*verifier.exact(vv)**s*comb(4096, f)*comb(4096-f, s)
    h, c_h, eta = anchor.source_factor(R)
    source = (verifier.exact(eta)*verifier.exact(core.BETA)**R
              *verifier.exact(z)**(-core.LIMIT))**exponent
    number = counts.profile_count(4096, R, f, f)
    assert number > 0 and number <= anchor.run_total(4096, R, s)
    logtwo = arb(2).log()

    def bound_from(r, terminal):
        rr = [[x/denominator for x in row] for row in r]
        tt = [[x/denominator for x in row] for row in terminal]
        return number*source*verifier.matmul(verifier.matpow(rr, 255), tt)[0][0]

    ordinary = bound_from(plain_r, plain_v)
    reports = []
    residue_checks = 0
    for m in ((1, 4, 8, 16) if args.include_16 else (1, 4, 8)):
        entries = []
        for row in degrees:
            entries.append(verifier.poly([sum((row[d] for d in range(k, 65, m)), arb(0))
                                          for k in range(m)]))
        p = [[entries[0], entries[1]], [entries[2], entries[3]]]
        t = [[entries[4]], [entries[5]]]
        rr, tt = cyclic_regions(verifier, p, t, m)
        for cyclic, evaluated in ((rr, plain_r), (tt, plain_v)):
            for i in range(2):
                for j in range(len(cyclic[0])):
                    summed = sum((cyclic[i][j][k] for k in range(m)), arb(0))
                    assert summed.overlaps(evaluated[i][j])
                    residue_checks += 1
        residue = R % m
        r = [[entry[residue] for entry in row] for row in rr]
        terminal = [[entry[residue] for entry in row] for row in tt]
        bound = bound_from(r, terminal)
        if m == 1:
            assert bound.overlaps(ordinary)
        guard()
        report = {"modulus": m, "target_residue": residue,
                  "log2_bound_interval": str(bound.log()/logtwo),
                  "gain_over_ordinary_bits_interval": str((ordinary/bound).log()/logtwo),
                  "selected_profile_proved_below_2_pow_minus_40": bool(bound < arb(2)**(-40)),
                  "seconds": time.monotonic()-started}
        reports.append(report)
        print(json.dumps({"phase": "modulus_complete", **report}), flush=True)
    guard()
    print(json.dumps({"scope": "one ordered actual r512 profile, identical in both matchings",
        "directed": True, "all_profiles_at_weight_covered": False,
        "whole_code_certificate": False, "fresh_marker_calibration": False,
        "parameters": {"source_weight": R, "ordered_profile": [[f, s], [f, s]],
                       "z": str(z), "alpha": str(alpha), "u": str(u), "v": str(vv),
                       "marker_source": marker_source, "distance_threshold": core.LIMIT},
        "indexed_sources_exact": number, "source_count_log2_interval": str(arb(number).log()/logtwo),
        "joint_missing": {"h": h, "c_h_exact": str(c_h), "eta_exact": str(eta)},
        "ordinary_cauchy_log2_interval": str(ordinary.log()/logtwo), "modular_results": reports,
        "tests": {**count_checks, "exact_cyclic_polynomial_coefficients": checks,
                  "sum_residues_replays_evaluated_entries": residue_checks,
                  "modulus1_replays_ordinary": True,
                  "packet_route_normalizer": "(1+u+v)^32 exactly",
                  "terminal_updates": 32767, "terminal_route_packets": 128},
        "precision_bits": args.precision, "seconds": time.monotonic()-started,
        "peak_process_bytes": peak_memory(), "generated_bytes": 0,
        "nonnegative_intersections": verifier.intersections,
        "anchor_pin": ANCHOR_PIN, "counts_pin": COUNT_PIN, "core_pin": anchor.CORE_PIN,
        "source_sha256": hashlib.sha256(Path(__file__).read_bytes()).hexdigest()}, indent=2))


if __name__ == "__main__":
    main()
