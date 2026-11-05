"""Directed all-r512 cover with exact/cyclic coefficient bands.

Reuse the successful sixteen-marker witness, but refine its coefficients:
cells16/48/80 truncate the full degree exactly at32/64/96 and retain singleton
degree modulo128; cells112/144/176/208 retain both residues modulo128. One
exact rectangle F<=256,S<=78 supplies every remaining target f217..256.
Exact coefficients normalize independently of their marker; marker232 is
used for that rectangle. All truncations remove positive aliases only.

Local cone caps are powered BEFORE route averaging, with no later alpha.
Each target divides by its own M4096(f,s)*u^f*v^s. All66048 ordered source
profiles and their distinct output-only terminals are counted. This proves
only the specified source weight, not a whole-code certificate.

Frozen helpers remain unchanged. Positive Arb arithmetic, exact count/toy
checks,900s/512MiB/20GiB guards, stdout only, and no generated output files.
"""

import argparse
from collections import defaultdict
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
DIRECTED_PIN = "3e4fe22d9caf9419f05d29ec001a7468b3a64a085844e34ac008f5537f439768"
FINE_PIN = "e4fcb2bbb4a73183d31f52ab2b86d5ee9e4906045dd045158fd75966d9aba1a2"
Z, ALPHA = F(17447, 25000), F(119, 1280)
M = 128
SCHEDULE = ((16, 0, 32, "mixed"), (48, 33, 64, "mixed"),
    (80, 65, 96, "mixed"), (112, 97, 128, "cyclic"),
    (144, 129, 160, "cyclic"), (176, 161, 192, "cyclic"),
    (208, 193, 216, "cyclic"), (232, 217, 256, "exact"))


def load_directed():
    path = Path(__file__).with_name("cheap_lane_cauchy_cells_directed.py")
    if hashlib.sha256(path.read_bytes()).hexdigest() != DIRECTED_PIN:
        raise RuntimeError("frozen directed-cell helper mismatch")
    spec = importlib.util.spec_from_file_location("hybrid_directed_base", path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def mixed_pack(v, degrees, full, modulus):
    stride = 2*full+1
    rows = [[v.arb(0) for _ in range((modulus-1)*stride+full+1)] for _ in range(6)]
    for (f, s), values in degrees.items():
        if f <= full:
            for k, value in enumerate(values):
                rows[k][f+stride*(s % modulus)] += value
    p = [v.poly(row) for row in rows]
    return [[p[0], p[1]], [p[2], p[3]]], [[p[4]], [p[5]]]


def mixed_product(v, left, right, full, modulus):
    """Exact full truncation and cyclic singleton folding, without carries."""
    stride = 2*full+1
    out = []
    for row in left:
        result = []
        for j in range(len(right[0])):
            v.guard()
            raw = v.poly([])
            for k, item in enumerate(row):
                raw += item*right[k][j]
            assert len(raw) <= 2*full+stride*(2*modulus-2)+1
            coefficients = [v.arb(0) for _ in range((modulus-1)*stride+full+1)]
            for s in range(modulus):
                for f in range(full+1):
                    value = raw[f+stride*s]
                    if s+modulus <= 2*modulus-2:
                        value += raw[f+stride*(s+modulus)]
                    coefficients[f+stride*s] = v.positive(value)
            result.append(v.poly(coefficients))
        out.append(result)
    return out


def mixed_tests(v, adapter):
    """Exact unwrapped P^3*T, then truncate/fold independently."""
    p = [[{(0, 0): F(1, 3), (1, 2): F(2, 5)}, {(2, 0): F(3, 7)}],
         [{(0, 1): F(5, 11)}, {(1, 0): F(7, 13), (2, 2): F(11, 17)}]]
    t = [[{(0, 1): F(13, 19)}], [{(2, 0): F(17, 23)}]]

    def exact_product(a, b):
        out = []
        for row in a:
            result = []
            for j in range(len(b[0])):
                value = defaultdict(F)
                for k, item in enumerate(row):
                    for (f, s), mass in item.items():
                        for (ff, ss), other in b[k][j].items():
                            value[f+ff, s+ss] += mass*other
                result.append(dict(value))
            out.append(result)
        return out

    expected = exact_product(exact_product(exact_product(p, p), p), t)
    checks = 0
    for full, modulus in ((0, 2), (1, 3), (3, 2)):
        stride = 2*full+1

        def pack(matrix):
            result = []
            for row in matrix:
                rr = []
                for item in row:
                    values = [F(0)]*((modulus-1)*stride+full+1)
                    for (f, s), mass in item.items():
                        if f <= full:
                            values[f+stride*(s % modulus)] += mass
                    rr.append(v.poly([v.exact(x) for x in values]))
                result.append(rr)
            return result

        pp, tt = pack(p), pack(t)
        actual = mixed_product(v, mixed_product(v,
            mixed_product(v, pp, pp, full, modulus), pp, full, modulus), tt, full, modulus)
        for i in range(2):
            for s in range(modulus):
                for f in range(full+1):
                    wanted = sum((mass for (ff, ss), mass in expected[i][0].items()
                                  if ff == f and ss % modulus == s), F(0))
                    adapter.exact_point(actual[i][0][f+stride*s], wanted, v.fmpq)
                    checks += 1
    return checks


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--precision", type=int, choices=(256, 384), default=256)
    args = parser.parse_args()
    if shutil.disk_usage(Path(__file__).resolve().anchor).free < 20*(1 << 30):
        raise RuntimeError("requires20GiB free on source volume")
    started = time.monotonic(); deadline = started+900
    old = load_directed(); base = old.load_base()
    fine = base.load("cheap_lane_cauchy_fine.py", FINE_PIN, "hybrid_exact")
    anchor = base.load("cheap_lane_cone_cauchy_anchors.py", base.ANCHOR_PIN, "hybrid_anchor")
    counts = base.load("cheap_lane_fast_profile_counts.py", base.COUNT_PIN, "hybrid_counts")
    screen = base.load("cheap_lane_cauchy_profile_screen.py", old.SCREEN_PIN, "hybrid_screen")
    counts.DEADLINE = deadline
    peak_memory = anchor.peak_memory_reader()

    def guard():
        if time.monotonic() > deadline:
            raise TimeoutError("900-second hybrid proof budget exceeded")
        if peak_memory() > 512*(1 << 20):
            raise MemoryError("512MiB hybrid proof budget exceeded")

    adapter = anchor.load_adapter(); core = adapter.load_verifier()
    core.Z, core.ALPHA = Z, ALPHA
    from flint import arb, ctx
    ctx.prec = args.precision

    class HybridVerifier(core.Verifier):
        def guard(self):
            guard()

    v = HybridVerifier(2, deadline)
    tests = {"mixed_exact_toy_coefficients": mixed_tests(v, adapter),
        "rectangle_exact_toy_coefficients": fine.exact_tests(v, adapter),
        "cyclic_exact_toy_coefficients": base.exact_tests(v, adapter),
        "orientation_checks": old.orientation_tests(), **counts.checks()}
    assert [f for _, lo, hi, _ in SCHEDULE for f in range(lo, hi+1)] == list(range(257))
    caps = v.local_caps(list(combinations_with_replacement(range(17), 4)))
    degrees, normalizers = base.local_degrees(v, caps, ALPHA, F(1), F(1))
    del caps
    regular, terminal = [None]*257, [None]*257
    witnesses = []; cyclic_dc = truncated_dc = 0
    logtwo = arb(2).log()
    print(json.dumps({"phase": "local_coefficients_complete", "seconds": time.monotonic()-started}), flush=True)
    for center, lo, hi, kind in SCHEDULE:
        guard()
        u, w = screen.marker_curve(center, old.CURVE_FULL, old.CURVE_SINGLE)
        up, wp = [v.exact(u)**k for k in range(33)], [v.exact(w)**k for k in range(33)]
        scaled = {(f, s): [up[f]*wp[s]*x for x in row] for (f, s), row in degrees.items()}
        plain = [sum((row[k] for row in scaled.values()), arb(0)) for k in range(6)]
        plain_p, plain_t = [[plain[0], plain[1]], [plain[2], plain[3]]], [[plain[4]], [plain[5]]]
        plain_r, plain_v = v.matpow(plain_p, 128), v.matmul(v.matpow(plain_p, 127), plain_t)
        if kind == "cyclic":
            full, single, stride = M-1, M-1, 2*M-1
            rr, tt = base.fold(v, scaled, M, M)
            product = lambda a, b: base.multiply(v, a, b, M, M)
        elif kind == "mixed":
            full, single, stride = hi, M-1, 2*hi+1
            rr, tt = mixed_pack(v, scaled, full, M)
            product = lambda a, b: mixed_product(v, a, b, full, M)
        else:
            full, single, stride = 256, 78, 513
            rr, tt = fine.rectangle_pack(v, scaled, full, single)
            product = lambda a, b: fine.rectangle_product(v, a, b, full, single)
        for _ in range(7):
            tt = product(rr, tt)
            rr = product(rr, rr)
        for table, evaluated in ((rr, plain_r), (tt, plain_v)):
            for i in range(2):
                for j in range(len(table[0])):
                    dc = sum((table[i][j][f+stride*s]
                              for s in range(single+1) for f in range(full+1)), arb(0))
                    if kind == "cyclic":
                        assert dc.overlaps(evaluated[i][j]); cyclic_dc += 1
                    else:
                        # Rounding may overlap equality; a disjoint reversed
                        # inequality would contradict positive truncation.
                        assert dc <= evaluated[i][j] or dc.overlaps(evaluated[i][j])
                        truncated_dc += 1
        for f in range(lo, hi+1):
            s = 512-2*f
            ff = f % M if kind == "cyclic" else f
            ss = s if kind == "exact" else s % M
            index = ff+stride*ss
            denominator = arb(comb(4096, f)*comb(4096-f, s))*v.exact(u)**f*v.exact(w)**s
            assert denominator > 0 and regular[f] is None
            regular[f] = [[v.positive(item[index])/denominator for item in row] for row in rr]
            terminal[f] = [[v.positive(item[index])/denominator for item in row] for row in tt]
        witnesses.append({"center": center, "first_full": lo, "last_full": hi, "method": kind,
            "full_limit_or_residue": full+1 if kind == "cyclic" else full,
            "single_limit_or_residue": single if kind == "exact" else single+1,
            "u": str(u), "v": str(w)})
        del rr, tt, table, evaluated, scaled, plain, plain_p, plain_t, plain_r, plain_v, product
        gc.collect(); guard()
        print(json.dumps({"phase": "coefficient_band_complete", **witnesses[-1],
            "seconds": time.monotonic()-started, "peak_process_bytes": peak_memory()}), flush=True)
    assert all(x is not None for x in regular+terminal)
    h, c_h, eta = anchor.source_factor(512)
    source = (v.exact(eta)*v.exact(core.BETA)**512*v.exact(Z)**(-core.LIMIT))**v.exact(ALPHA)
    run_bounds, run_counts = [arb(0) for _ in range(513)], [0]*513
    ordered = 0; worst = selected = None
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
            terms = [(f0, f1, v.positive(v.matmul(v.matmul(shared, r0), terminal[f1])[0][0]))]
            if f0 != f1:
                terms.append((f1, f0, v.positive(v.matmul(v.matmul(r1, shared), terminal[f0])[0][0])))
            run_counts[runs] += len(terms)*number; ordered += len(terms)
            for a, b, value in terms:
                mass = number*value; run_bounds[runs] += mass
                if worst is None or mass.upper() > worst[2].upper():
                    worst = (a, b, mass)
                if a == b == 223:
                    selected = source*mass
        if f0 % 32 == 0:
            print(json.dumps({"phase": "profile_rows_complete", "through_full0": f0,
                "ordered_profiles": ordered, "seconds": time.monotonic()-started}), flush=True)
    assert ordered == 66048 and selected is not None
    assert sum(run_counts) == comb(8192, 512)
    for c in range(1, 513):
        assert run_counts[c] == counts.run_total(4096, 512, c)
    total = source*sum(run_bounds, arb(0)); assert total > 0
    worst_run = max(range(1, 513), key=lambda c: run_bounds[c].upper())
    parameters = {"source_weight": 512, "z": str(Z), "alpha": str(ALPHA),
        "distance_threshold": core.LIMIT, "bands": witnesses}
    witness_hash = hashlib.sha256(json.dumps(parameters, sort_keys=True, separators=(",", ":")).encode()).hexdigest()
    guard()
    print(json.dumps({"scope": "all actual weight512 supports, matching profiles, and state paths",
        "directed": True, "all_profiles_at_weight_covered": True, "whole_code_certificate": False,
        "log2_union_bound_interval": str(total.log()/logtwo),
        "weight512_below_2_pow_minus_40": bool(total < arb(2)**(-40)),
        "weight512_below_2_pow_minus_50": bool(total < arb(2)**(-50)),
        "selected_f223_EXACT_log2_interval": str(selected.log()/logtwo),
        "worst_ordered_profile": {"full": list(worst[:2]),
            "log2_bound_interval": str((source*worst[2]).log()/logtwo)},
        "worst_run_count": worst_run,
        "sampled_run_unions": [{"runs": c, "log2_bound_interval": str((source*run_bounds[c]).log()/logtwo)}
            for c in sorted({1, 16, 32, 64, 66, 128, 256, 480, 512, worst_run})],
        "parameters": parameters, "joint_missing": {"h": h, "c_h": str(c_h), "eta": str(eta)},
        "tests": {**tests, "local_route_normalizers": normalizers,
            "cyclic_dc_overlaps": cyclic_dc, "truncated_dc_inequality_checks": truncated_dc,
            "ordered_profiles": ordered, "all512_run_count_sums_exact": True,
            "total_indexed_supports_equals_binomial8192_512": True, "terminal_updates": 32767},
        "precision_bits": args.precision, "seconds": time.monotonic()-started,
        "peak_process_bytes": peak_memory(), "budget_seconds": 900,
        "nonnegative_intersections": v.intersections, "generated_bytes": 0,
        "witness_sha256": witness_hash, "directed_pin": DIRECTED_PIN, "fine_pin": FINE_PIN,
        "base_pin": old.BASE_PIN, "screen_pin": old.SCREEN_PIN,
        "anchor_pin": base.ANCHOR_PIN, "counts_pin": base.COUNT_PIN,
        "adapter_pin": anchor.ADAPTER_PIN, "core_pin": anchor.CORE_PIN,
        "source_sha256": hashlib.sha256(Path(__file__).read_bytes()).hexdigest()}, indent=2))


if __name__ == "__main__":
    main()
