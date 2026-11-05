"""FLOAT r512 all-profile Cauchy screen with exact indexed source counts.

Keep one rational z,alpha and one deterministic marker curve per run.
For each f=0..256, s=512-2f, use
 u_f=curve_full*f/(3584+f), v_f=curve_single*s/(3584+f).
Replace a zero marker by the positive rational2^-30; do not treat it as
an exact coefficient boundary. Default corrections are the recorded
r512 calibrated ratios. Local cone caps are powered BEFORE route averaging.

Cache only the full33-by33-by6 LOCAL packet polynomial. At each marker
compute the true finite matrices B=P^128 and V=P^127 T. For ordered
profile(f0,f1), the chain is e0(B0 B1)^127 B0 V1 divided by
[M_4096(f0,s0)u0^f0 v0^s0 M_4096(f1,s1)u1^f1 v1^s1]^128.
The source factor (Omega_512 z^-734003)^alpha appears once. Multiply
by exact profile_count(4096,512,f0,f1), not the old run-total upper count.
Both orientations retain their own output-only terminal. Summed source
counts must equal C(8192,512), and every run total is checked exactly.

This is a FLOAT diagnostic, including the final union. A negative result
needs a separate directed replay. There is no Perron approximation, extra
alpha at regional boundaries, automatic parameter grid, or output file.
180-second cooperative guard,256MiB Windows process guard,20GiB free.
"""

import argparse
from fractions import Fraction as F
import hashlib
import heapq
import importlib.util
from itertools import combinations_with_replacement
import json
import math
from math import comb
from pathlib import Path
import shutil
import sys
import time

sys.dont_write_bytecode = True
CALIBRATION_PIN = "f0358e6296357f13e1d0b18b676e06b9672a1ac4ecbb4bbe65872a3391e8054e"
COUNT_PIN = "4ce61f179a8bf41e5d769e955a9efeb437e46af81f57b0d784c4cfc57958d8e8"
WEIGHT, PAIRS, MEMORY_LIMIT = 512, 4096, 256*(1 << 20)
ZERO_MARKER = F(1, 1 << 30)


def load(name, pin, tag):
    path = Path(__file__).with_name(name)
    if hashlib.sha256(path.read_bytes()).hexdigest() != pin:
        raise RuntimeError("frozen helper hash mismatch: "+name)
    spec = importlib.util.spec_from_file_location(tag, path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def normalized_matrix(entries, scale=0.0):
    largest = max(entries)
    if not largest > 0 or any(not math.isfinite(x) or x < 0 for x in entries):
        raise ArithmeticError("invalid FLOAT matrix")
    return tuple(x/largest for x in entries), scale+math.log(largest)


def matrix_product(left, right):
    a, la = left; b, lb = right
    return normalized_matrix((a[0]*b[0]+a[1]*b[2], a[0]*b[1]+a[1]*b[3],
                              a[2]*b[0]+a[3]*b[2], a[2]*b[1]+a[3]*b[3]), la+lb)


def matrix_power(value, exponent):
    result = ((1.0, 0.0, 0.0, 1.0), 0.0)
    while exponent:
        if exponent & 1:
            result = matrix_product(result, value)
        exponent //= 2
        if exponent:
            value = matrix_product(value, value)
    return result


def matrix_vector(matrix, vector):
    a, la = matrix; v, lv = vector
    return normalized_matrix((a[0]*v[0]+a[1]*v[1], a[2]*v[0]+a[3]*v[1]), la+lv)


def first_log(vector):
    values, scale = vector
    if values[0] <= 0:
        return -math.inf
    return scale+math.log(values[0])


def paired_chains(region0, region1):
    """One matrix power supplies both DIFFERENT terminal orientations."""
    b0, v0, d0 = region0[:3]; b1, v1, d1 = region1[:3]
    combined = matrix_power(matrix_product(b0, b1), 127)
    forward = matrix_vector(combined, matrix_vector(b0, v1))
    # (B1 B0)^127 B1 = B1 (B0 B1)^127. Do not swap terminal scales.
    reverse = matrix_vector(b1, matrix_vector(combined, v0))
    denominator = 128*(d0+d1)
    return first_log(forward)-denominator, first_log(reverse)-denominator


def logadd(x, y):
    if x == -math.inf:
        return y
    if y == -math.inf:
        return x
    high, low = (x, y) if x >= y else (y, x)
    return high+math.log1p(math.exp(low-high))


def marker_curve(full, curve_full, curve_single):
    single = WEIGHT-2*full
    u = curve_full*F(full, PAIRS-WEIGHT+full)
    v = curve_single*F(single, PAIRS-WEIGHT+full)
    return (u if u else ZERO_MARKER), (v if v else ZERO_MARKER)


def build_regions(coefficients, curve_full, curve_single, np, guard):
    """FLOAT-only seam: replace this routine with directed evaluation later."""
    indices = [(f, s) for f in range(33) for s in range(33-f)]
    full_degree = np.asarray([x[0] for x in indices], dtype=float)
    single_degree = np.asarray([x[1] for x in indices], dtype=float)
    local = np.asarray([coefficients[f, s] for f, s in indices])
    regions = []
    for full in range(257):
        guard()
        single = WEIGHT-2*full
        u, v = marker_curve(full, curve_full, curve_single)
        lu, lv = math.log(u.numerator)-math.log(u.denominator), math.log(v.numerator)-math.log(v.denominator)
        monomials = np.exp(full_degree*lu+single_degree*lv)
        entries = np.sum(local*monomials[:, None], axis=0)
        p = normalized_matrix(tuple(float(x) for x in entries[:4]))
        t = normalized_matrix(tuple(float(x) for x in entries[4:]))
        b = matrix_power(p, 128)
        terminal = matrix_vector(matrix_power(p, 127), t)
        denominator = (math.log(comb(PAIRS, full)*comb(PAIRS-full, single))
                       +full*lu+single*lv)
        regions.append((b, terminal, denominator, u, v))
    return regions


def chain_tests():
    """Exact noncommuting small powers check forward/reverse terminal order."""
    def mm(a, b):
        return [[sum((a[i][k]*b[k][j] for k in range(len(b))), F(0))
                 for j in range(len(b[0]))] for i in range(len(a))]
    p = [[F(1, 3), F(2, 5)], [F(3, 7), F(4, 11)]]
    q = [[F(5, 13), F(6, 17)], [F(7, 19), F(8, 23)]]
    t = [[F(9, 29)], [F(10, 31)]]
    u = [[F(11, 37)], [F(12, 41)]]
    pf = normalized_matrix(tuple(float(x) for row in p for x in row))
    qf = normalized_matrix(tuple(float(x) for row in q for x in row))
    tf = normalized_matrix(tuple(float(row[0]) for row in t))
    uf = normalized_matrix(tuple(float(row[0]) for row in u))
    checks = 0
    for exponent in (1, 2, 3):
        pq = [[F(int(i == j)) for j in range(2)] for i in range(2)]
        qp = [[F(int(i == j)) for j in range(2)] for i in range(2)]
        for _ in range(exponent):
            pq = mm(pq, mm(p, q)); qp = mm(qp, mm(q, p))
        expected = (mm(mm(pq, p), u)[0][0], mm(mm(qp, q), t)[0][0])
        shared = matrix_power(matrix_product(pf, qf), exponent)
        got = (first_log(matrix_vector(shared, matrix_vector(pf, uf))),
               first_log(matrix_vector(qf, matrix_vector(shared, tf))))
        for actual, exact in zip(got, expected):
            assert abs(actual-math.log(float(exact))) < 1e-12
            checks += 1
    return checks


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--z", default="7/10")
    parser.add_argument("--alpha", default="119/1280")
    parser.add_argument("--curve-full", default="1562175/990029")
    parser.add_argument("--curve-single", default="1140210/880327")
    args = parser.parse_args()
    z, alpha, curve_full, curve_single = map(F, (args.z, args.alpha, args.curve_full, args.curve_single))
    if not 0 < z < 1 or not 0 < alpha <= 1:
        parser.error("require0<z<1 and0<alpha<=1")
    if not all(F(1, 256) <= x <= 256 for x in (curve_full, curve_single)):
        parser.error("curve multipliers must lie in[1/256,256]")
    if any(max(x.numerator.bit_length(), x.denominator.bit_length()) > 64
           for x in (z, alpha, curve_full, curve_single)):
        parser.error("rational witnesses must fit64 bits")
    if shutil.disk_usage(Path(__file__).resolve().anchor).free < 20*(1 << 30):
        raise RuntimeError("requires20GiB free on source volume")
    started = time.monotonic(); deadline = started+180
    calibration = load("cheap_lane_cauchy_marker_calibration.py", CALIBRATION_PIN, "all_profile_calibration")
    gradient = calibration.load_gradient(); anchor = gradient.load_anchor()
    peak_memory = anchor.peak_memory_reader()

    def guard():
        if time.monotonic() > deadline:
            raise TimeoutError("180-second all-profile screen exhausted; no automatic retry")
        if peak_memory() > MEMORY_LIMIT:
            raise MemoryError("256MiB process-memory budget exceeded")

    def progress(phase, **values):
        print(json.dumps({"phase": phase, "seconds": time.monotonic()-started, **values}), flush=True)

    counts = load("cheap_lane_fast_profile_counts.py", COUNT_PIN, "all_profile_exact_counts")
    counts.DEADLINE = deadline
    count_checks = counts.checks()
    checks = chain_tests()
    import numpy as np
    from flint import arb, ctx
    ctx.prec = 256
    adapter = anchor.load_adapter(); core = adapter.load_verifier()

    class ScreenVerifier(core.Verifier):
        def guard(self):
            guard()

    verifier = ScreenVerifier(2, deadline)
    core.Z, core.ALPHA = z, alpha
    shapes = list(combinations_with_replacement(range(17), 4))
    caps = verifier.local_caps(shapes)
    exponent = verifier.exact(alpha)
    powered = {occ: tuple(verifier.positive(x)**exponent for x in cap) for occ, cap in caps.items()}
    del caps
    coefficients = calibration.float_coefficients(powered, np, guard)
    del powered
    regions = build_regions(coefficients, curve_full, curve_single, np, guard)
    h, c, eta = anchor.source_factor(WEIGHT)
    # This is ONLY a source factor: the old selected-profile run count is absent.
    source_ball = (verifier.exact(eta)*verifier.exact(core.BETA)**WEIGHT
                   *verifier.exact(z)**(-core.LIMIT))**exponent
    log_source = float(source_ball.log().mid())
    progress("marginals_ready", marginal_profiles=len(regions), local_coefficient_shape=list(coefficients.shape))
    run_counts = [0]*513; run_logs = [-math.inf]*513
    total_sources = 0; total_log = -math.inf; chains = positive = 0
    top = []

    def add_profile(f0, f1, number, chain):
        nonlocal total_sources, total_log, chains, positive
        runs = WEIGHT-f0-f1
        weighted = log_source+math.log(number)+chain
        total_sources += number; run_counts[runs] += number
        total_log = logadd(total_log, weighted)
        run_logs[runs] = logadd(run_logs[runs], weighted)
        chains += 1; positive += int(weighted >= 0)
        item = (weighted, f0, f1, math.log(number), chain)
        if len(top) < 12:
            heapq.heappush(top, item)
        elif weighted > top[0][0]:
            heapq.heapreplace(top, item)

    for f0 in range(257):
        guard()
        for f1 in range(f0, 257):
            number = counts.profile_count(PAIRS, WEIGHT, f0, f1)
            if not number:
                continue
            forward, reverse = paired_chains(regions[f0], regions[f1])
            add_profile(f0, f1, number, forward)
            if f0 != f1:
                add_profile(f1, f0, number, reverse)
        if f0 % 64 == 0:
            progress("profile_rows", first_full_completed=f0, ordered_chains=chains)
    expected_sources = comb(8192, WEIGHT)
    assert total_sources == expected_sources
    assert chains == 66048
    assert all(run_counts[runs] == counts.run_total(PAIRS, WEIGHT, runs)
               for runs in range(513))
    assert abs(total_log-math.log(sum(math.exp(x-max(run_logs)) for x in run_logs if x != -math.inf))-max(run_logs)) < 1e-8
    logtwo = math.log(2)
    worst_runs = sorted(range(1, 513), key=lambda c: run_logs[c], reverse=True)[:12]
    selected_runs = (1, 2, 4, 8, 16, 32, 64, 128, 256, 384, 480, 512)
    guard()
    print(json.dumps({"scope": "all state paths and all indexed r512 source supports; FLOAT screen only",
        "directed": False, "whole_code_certificate": False, "negative_is_not_a_certificate": True,
        "parameters": {"source_weight": WEIGHT, "z": str(z), "alpha": str(alpha),
                       "curve_full": str(curve_full), "curve_single": str(curve_single),
                       "zero_marker_replacement": str(ZERO_MARKER), "distance_threshold": core.LIMIT},
        "log2_union_bound_FLOAT": total_log/logtwo,
        "powered_source_log2_FLOAT": log_source/logtwo,
        "ordered_profile_chains": chains, "positive_profile_terms": positive,
        "positive_run_unions": sum(x >= 0 for x in run_logs[1:]),
        "indexed_source_count": total_sources,
        "source_count_log2": math.log(total_sources)/logtwo,
        "top_weighted_profiles": [{"full0": a, "full1": b, "single0": WEIGHT-2*a,
            "single1": WEIGHT-2*b, "runs": WEIGHT-a-b,
            "log2_bound_FLOAT": value/logtwo, "log2_source_count": n/logtwo,
            "log2_conditional_chain_FLOAT": chain/logtwo}
            for value, a, b, n, chain in sorted(top, reverse=True)],
        "top_run_unions": [{"runs": c, "log2_bound_FLOAT": run_logs[c]/logtwo} for c in worst_runs],
        "selected_run_unions": [{"runs": c, "log2_bound_FLOAT": run_logs[c]/logtwo} for c in selected_runs],
        "tests": {**count_checks, "noncommuting_terminal_chain_checks": checks,
                  "exact_source_sum": "C(8192,512)", "all_exact_run_totals": True,
                  "both_terminal_orientations_retained": True, "updates": 32767},
        "seconds": time.monotonic()-started, "peak_process_bytes": peak_memory(), "generated_bytes": 0,
        "calibration_pin": CALIBRATION_PIN, "count_pin": COUNT_PIN,
        "source_sha256": hashlib.sha256(Path(__file__).read_bytes()).hexdigest()}, indent=2))


if __name__ == "__main__":
    main()
