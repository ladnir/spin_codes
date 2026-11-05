"""Bounded FLOAT diagnostic for the route-independent positive-vector bound.

Fix rational z and alpha. The frozen cone supplies six nonnegative caps
(Q,A,E,U,T0,U) for every sorted four-block occupancy ell in {0,...,16}^4.
Raise these caps to alpha, without changing the refresh distribution.
For d=sum(ell), t,w>0 and p=(1,w), define the least common constants

 kappa = max_ell max(Q^alpha+A^alpha*w, E^alpha/w+U^alpha)/t^d,
 kappaT = max_ell max(T0^alpha,U^alpha/w)/t^d.

Then M_ell*p <= kappa*t^d*p and T_ell <= kappaT*t^d*p. There are32767
updates and one output-only terminal, with exactly256*r allowed coordinates
at source size r. Hence a band I has the route-independent upper bound

 z^(-alpha*734003) kappa^32767 kappaT
   * sum_(r in I) C(8192,r) Omega_r^alpha t^(256*r).

Each Omega_r is the proved joint-missing source factor, evaluated from
exact integers/rationals. There is no profile count, route average, residual
potential, extra terminal update, or uniform full-field substitution.

The logarithm is convex in(log(t),log(w)). Nested golden-section searches
minimize this objective on an explicit finite box. The output is only a
FLOAT proposal: even a negative result needs independent directed replay;
search convergence is not a proof of the exact/global minimum. A positive
bound is not a counterexample. Uniform conditioning includes rare highly
concentrated routes, which can prevent this method from closing a band.

No SciPy, automatic z/alpha scan, generated data, or output files. One run
has a180-second cooperative deadline,256MiB process guard,20GiB free-space
guard, and at most four configurable inclusive bands. The default bands
start above the analytically excluded1024..4300 range. Optional overlapping
bands retain their harmless overcount when summed; it is reported explicitly.
"""

import argparse
from fractions import Fraction as F
import hashlib
import importlib.util
from itertools import combinations_with_replacement, product
import json
import math
from math import comb
import os
from pathlib import Path
import shutil
import sys
import time

sys.dont_write_bytecode = True
os.environ["OPENBLAS_NUM_THREADS"] = "1"
os.environ["OMP_NUM_THREADS"] = "1"
ANCHOR_PIN = "319d169a691abe062e56516b96d8e1800b54040375dbb5fc7981984cbb02eba9"
CORE_PIN = "7ccf0cdab823ef15fdd5edd6156d4b49d9a000f2ae28c37d5d43f4b6ea341f01"
ADAPTER_PIN = "83e00b6f6d148c8058428295b4132d401df08cfccd20cadcda1c1ba82df7b478"
DEFAULT_BANDS = "4301:6143,6144:7551"
MEMORY_LIMIT, SECONDS = 256*(1 << 20), 180
UPDATES, COORDINATES, GROUPS = 32767, 256, 8192
LN2 = math.log(2)


def load_anchor():
    path = Path(__file__).with_name("cheap_lane_cone_cauchy_anchors.py")
    if hashlib.sha256(path.read_bytes()).hexdigest() != ANCHOR_PIN:
        raise RuntimeError("frozen source-factor helper hash mismatch")
    spec = importlib.util.spec_from_file_location("route_free_anchor", path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    if module.CORE_PIN != CORE_PIN or module.ADAPTER_PIN != ADAPTER_PIN:
        raise RuntimeError("unexpected frozen cone/adapter dependency")
    return module


def log_fraction(value):
    return math.log(value.numerator)-math.log(value.denominator)


def logsum(values, np):
    high = float(np.max(values))
    return high+math.log(float(np.sum(np.exp(values-high))))


def parse_bands(text):
    bands = []
    for item in text.split(","):
        parts = item.split(":")
        if len(parts) != 2:
            raise ValueError("bands must be comma-separated inclusive lo:hi pairs")
        lo, hi = map(int, parts)
        if not 1 <= lo <= hi <= 8191:
            raise ValueError("each band must lie in1..8191")
        bands.append((lo, hi))
    if not 1 <= len(bands) <= 4 or len(set(bands)) != len(bands):
        raise ValueError("require one to four distinct bands")
    return bands


def exact_tests(anchor):
    """Tiny positive products test the theorem independently of the cone."""
    matrices = (
        (0, ((F(1), F(0)), (F(1, 7), F(1, 3)))),
        (1, ((F(1, 3), F(2, 5)), (F(3, 7), F(4, 11)))),
        (2, ((F(2, 9), F(1, 5)), (F(1, 8), F(3, 13)))),
    )
    terminals = ((0, (F(1), F(2, 3))), (1, (F(3, 5), F(1, 2))),
                 (2, (F(1, 4), F(2, 7))))
    checks = 0
    for t, w in ((F(1), F(1)), (F(2, 3), F(5, 7)), (F(7, 5), F(3, 2))):
        p = (F(1), w)
        kappa = max(sum(row[j]*p[j] for j in range(2))/(p[i]*t**d)
                    for d, matrix in matrices for i, row in enumerate(matrix))
        kt = max(value[i]/(p[i]*t**d) for d, value in terminals for i in range(2))
        for length in (0, 1, 2, 3):
            for sequence in product(matrices, repeat=length):
                for td, terminal in terminals:
                    value = terminal
                    for _, matrix in reversed(sequence):
                        value = tuple(sum(matrix[i][j]*value[j] for j in range(2))
                                      for i in range(2))
                    total_degree = td+sum(d for d, _ in sequence)
                    assert value[0] <= kappa**length*kt*t**total_degree
                    checks += 1
    # The source factor changes h at multiples of256; never interpolate eta.
    for r in (1, 255, 256, 257, 511, 512, 7551, 8191):
        h, c, eta = anchor.source_factor(r)
        expected_c = F(256*(h+1)-r, 256*(h+1))
        expected_eta = F(255+sum(comb(r, j) for j in range(1, h+1)), 255**2)/expected_c
        assert h == r//256 and c == expected_c and eta == expected_eta and c > 0
        checks += 1
    for n in (4, 7, 12):
        value = 1; total = value
        for r in range(n):
            value = value*(n-r)//(r+1)
            assert value == comb(n, r+1)
            total += value
            checks += 1
        assert total == 2**n
    assert UPDATES+1 == 32768 and COORDINATES*GROUPS == 32768*64
    return checks


def source_table(bands, anchor, beta, alpha, np, guard):
    requested = {r for lo, hi in bands for r in range(lo, hi+1)}
    first, last = min(requested), max(requested)
    number = comb(GROUPS, first)
    table = {}; total = 0
    logbeta = log_fraction(beta)
    for r in range(first, last+1):
        guard()
        if r in requested:
            h, c, eta = anchor.source_factor(r)
            table[r] = (math.log(number)+alpha*(r*logbeta+log_fraction(eta)),
                        math.log(number), alpha*(r*logbeta+log_fraction(eta)), h, c, eta)
            total += number
        if r in (first, last) or r % 256 == 0:
            assert number == comb(GROUPS, r)
        if r != last:
            numerator = number*(GROUPS-r)
            assert numerator % (r+1) == 0
            number = numerator//(r+1)
    return table, {"unique_weights": len(requested),
                   "exact_unique_support_count_log2_FLOAT": math.log(total)/LN2,
                   "exact_unique_support_count_sha256": hashlib.sha256(str(total).encode()).hexdigest()}


def golden(function, lo, hi, iterations, guard):
    """Bounded convex search; always return an actually evaluated witness."""
    phi = (math.sqrt(5)-1)/2
    a, b = lo, hi
    x1, x2 = b-phi*(b-a), a+phi*(b-a)
    y1, y2 = function(x1), function(x2)
    candidates = [(lo, function(lo)), (hi, function(hi)), (x1, y1), (x2, y2)]
    if lo < 0 < hi:
        candidates.append((0.0, function(0.0)))
    for _ in range(iterations):
        guard()
        if y1 <= y2:
            b, x2, y2 = x2, x1, y1
            x1 = b-phi*(b-a); y1 = function(x1)
            candidates.append((x1, y1))
        else:
            a, x1, y1 = x1, x2, y2
            x2 = a+phi*(b-a); y2 = function(x2)
            candidates.append((x2, y2))
    point, value = min(candidates, key=lambda item: item[1])
    return point, value, (a, b)


class Objective:
    def __init__(self, logcaps, shapes, source, band, alpha, z, limit, np, guard):
        self.np, self.guard = np, guard
        self.logcaps, self.shapes = logcaps, shapes
        self.degree = np.asarray([sum(occ) for occ in shapes], dtype=float)
        lo, hi = band
        self.weights = np.arange(lo, hi+1, dtype=float)
        self.base = np.asarray([source[r][0] for r in range(lo, hi+1)])
        self.threshold = -alpha*limit*log_fraction(z)
        self.evaluations = 0

    def rows(self, x, y):
        np, a = self.np, self.logcaps
        tilt = self.degree*x
        regular = np.column_stack((np.logaddexp(a[:, 0], a[:, 1]+y)-tilt,
                                   np.logaddexp(a[:, 2]-y, a[:, 3])-tilt))
        terminal = np.column_stack((a[:, 4]-tilt, a[:, 5]-y-tilt))
        return regular, terminal

    def local(self, x, y):
        self.guard(); self.evaluations += 1
        regular, terminal = self.rows(x, y)
        return UPDATES*float(self.np.max(regular))+float(self.np.max(terminal))

    def source(self, x):
        return logsum(self.base+COORDINATES*self.weights*x, self.np)

    def value(self, x, y):
        return self.threshold+self.local(x, y)+self.source(x)

    def report(self, x, y):
        regular, terminal = self.rows(x, y)
        kr, kt = float(self.np.max(regular)), float(self.np.max(terminal))
        source_values = self.base+COORDINATES*self.weights*x
        source = logsum(source_values, self.np)
        value = self.threshold+UPDATES*kr+kt+source

        def constraints(array, largest):
            flat = array.reshape(-1)
            order = self.np.argsort(flat)[-12:][::-1]
            return {"within_1e_minus_8_nats": int(self.np.sum(flat >= largest-1e-8)),
                    "largest_rows": [{"shape": self.shapes[int(i)//2], "row": int(i)%2,
                        "degree": int(self.degree[int(i)//2]),
                        "log2_ratio": float(flat[i])/LN2,
                        "log2_slack_to_max": (largest-float(flat[i]))/LN2} for i in order]}

        order = self.np.argsort(source_values)[-5:][::-1]
        return {"log_t": x, "log_w": y, "t_FLOAT": math.exp(x), "w_FLOAT": math.exp(y),
                "log2_bound_FLOAT": value/LN2,
                "components_log2_FLOAT": {"output_threshold": self.threshold/LN2,
                    "updates": UPDATES*kr/LN2, "terminal": kt/LN2,
                    "source_and_coordinate_tilt_sum": source/LN2},
                "regular_constraints": constraints(regular, kr),
                "terminal_constraints": constraints(terminal, kt),
                "largest_source_terms": [{"r": int(self.weights[i]),
                    "log2_term_FLOAT": float(source_values[i])/LN2,
                    "normalized_mass_FLOAT": math.exp(float(source_values[i])-source)} for i in order]}


def optimize(objective, xradius, yradius, iterations, guard):
    witnesses = {}

    def outer(x):
        if x not in witnesses:
            y, local, bracket = golden(lambda yy: objective.local(x, yy),
                                       -yradius, yradius, iterations, guard)
            value = objective.threshold+local+objective.source(x)
            witnesses[x] = (y, value, bracket)
        return witnesses[x][1]

    x, _, bracket = golden(outer, -xradius, xradius, iterations, guard)
    y, _, inner_bracket = witnesses[x]
    best = objective.report(x, y)
    zero = objective.report(0.0, 0.0)
    if zero["log2_bound_FLOAT"] < best["log2_bound_FLOAT"]:
        best = zero
    best["optimizer"] = {"method": "nested convex golden-section on an explicit finite box",
        "iterations_per_search": iterations, "outer_witnesses": len(witnesses),
        "local_evaluations": objective.evaluations,
        "log_t_box": [-xradius, xradius], "log_w_box": [-yradius, yradius],
        "last_outer_bracket": bracket, "selected_inner_bracket": inner_bracket,
        "box_boundary_selected": abs(best["log_t"]) >= xradius-1e-7
                                 or abs(best["log_w"]) >= yradius-1e-7,
        "status": "bounded search complete; neither directed nor a certified global optimum"}
    return best, zero


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--bands", default=DEFAULT_BANDS)
    parser.add_argument("--z", type=F, default=F(17447, 25000))
    parser.add_argument("--alpha", type=F, default=F(119, 1280))
    parser.add_argument("--iterations", type=int, default=48)
    parser.add_argument("--log-t-radius", type=float, default=8.0)
    parser.add_argument("--log-w-radius", type=float, default=64.0)
    args = parser.parse_args()
    bands = parse_bands(args.bands)
    if not 0 < args.z < 1 or not 0 < args.alpha <= 1:
        parser.error("require0<z<1 and0<alpha<=1")
    if not 32 <= args.iterations <= 64:
        parser.error("iterations must be32..64")
    if not 1 <= args.log_t_radius <= 32 or not 1 <= args.log_w_radius <= 128:
        parser.error("finite log radii must be in1..32 and1..128")
    free = shutil.disk_usage(Path(__file__).resolve().anchor).free
    if free < 20*(1 << 30):
        raise RuntimeError("requires20GiB free on source volume")
    started = time.monotonic(); deadline = started+SECONDS
    anchor = load_anchor(); peak_memory = anchor.peak_memory_reader()

    def guard():
        if time.monotonic() > deadline:
            raise TimeoutError("180-second route-free diagnostic expired; no automatic retry")
        if peak_memory() > MEMORY_LIMIT:
            raise MemoryError("256MiB process-memory budget exceeded")

    import numpy as np
    from flint import ctx
    ctx.prec = 192
    tests = exact_tests(anchor)
    adapter = anchor.load_adapter(); core = adapter.load_verifier()
    core.Z, core.ALPHA = args.z, args.alpha

    class LocalVerifier(core.Verifier):
        def guard(self):
            guard()

    verifier = LocalVerifier(2, deadline)
    shapes = list(combinations_with_replacement(range(17), 4))
    assert len(shapes) == 4845
    caps = verifier.local_caps(shapes)
    alpha = float(args.alpha)
    logcaps = np.asarray([[(-math.inf if value == 0 else alpha*float(value.log().mid()))
                           for value in caps[occ]] for occ in shapes])
    if np.any(np.isnan(logcaps)) or np.any(logcaps == math.inf):
        raise ArithmeticError("invalid FLOAT local caps")
    assert np.array_equal(logcaps[:, 3], logcaps[:, 5])
    del caps
    source, source_control = source_table(bands, anchor, core.BETA, alpha, np, guard)
    print(json.dumps({"phase": "local_caps_and_source_ready", "seconds": time.monotonic()-started,
        "shapes": len(shapes), "exact_tests": tests, "source_control": source_control}), flush=True)
    results = []
    for band in bands:
        guard()
        objective = Objective(logcaps, shapes, source, band, alpha, args.z, core.LIMIT, np, guard)
        best, zero = optimize(objective, args.log_t_radius, args.log_w_radius, args.iterations, guard)
        boundary = []
        for r in sorted({band[0], band[1]}):
            _, count_log, omega_log, h, c, eta = source[r]
            boundary.append({"r": r, "h": h, "c_h_exact": str(c),
                "log2_eta_FLOAT": log_fraction(eta)/LN2,
                "log2_source_count_FLOAT": count_log/LN2,
                "log2_powered_Omega_FLOAT": omega_log/LN2})
        report = {"band_inclusive": band, "optimized": best, "zero_parameter_control": zero,
            "source_boundary_controls": boundary,
            "improvement_bits_FLOAT": zero["log2_bound_FLOAT"]-best["log2_bound_FLOAT"],
            "seconds": time.monotonic()-started}
        # Empty/full shapes give a parameter-independent LOWER bound on this
        # majorant, not on the code's failure probability. This is only a
        # FLOAT check here; the accompanying analytic obstruction is separate.
        logg = log_fraction(args.z+(1-args.z)/256)
        lower_terms = [source[r][0]+alpha*(-64*r*math.log(256)+256*r*logg
                       -core.LIMIT*log_fraction(args.z)) for r in range(band[0], band[1]+1)]
        report["empty_full_shape_majorant_lower_log2_FLOAT"] = logsum(np.asarray(lower_terms), np)/LN2
        report["known_obstructed_1024_to_4300_weights"] = max(0, min(band[1], 4300)-max(band[0], 1024)+1)
        results.append(report)
        print(json.dumps({"phase": "band_complete", **report}), flush=True)
    guard()
    summed = logsum(np.asarray([r["optimized"]["log2_bound_FLOAT"]*LN2 for r in results]), np)/LN2
    print(json.dumps({"scope": "route-independent source-band majorants; FLOAT diagnostic only",
        "directed": False, "certificate": False, "whole_code_certificate": False,
        "failure_is_not_counterexample": True, "z": str(args.z), "alpha": str(args.alpha),
        "bands_inclusive": bands, "band_sum_log2_FLOAT": summed,
        "band_endpoint_overcount": sum(hi-lo+1 for lo, hi in bands)-source_control["unique_weights"],
        "conditional_routes": "all routes, including rare concentrated ones",
        "tests": {"exact_small_comparisons": tests, "all_occupancy_shapes": len(shapes),
                  "updates": UPDATES, "terminal_emissions": 1, "total_coordinates": "256*r"},
        "dependencies": {"anchor": ANCHOR_PIN, "adapter": ADAPTER_PIN, "core": CORE_PIN},
        "seconds": time.monotonic()-started, "peak_process_bytes": peak_memory(),
        "free_bytes_preflight": free, "generated_bytes": 0,
        "source_sha256": hashlib.sha256(Path(__file__).read_bytes()).hexdigest()}, indent=2), flush=True)


if __name__ == "__main__":
    main()
