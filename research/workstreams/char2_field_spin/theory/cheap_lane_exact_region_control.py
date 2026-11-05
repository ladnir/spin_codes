"""FLOAT control: exact regional route coefficients, only zero/full states.

Keep the pinned pooling pilot's raw caps, z=.85, alpha=.05, quiet mark,
floor mark1 and consecutive256-group source profiles. Restrict all actual
packet masks to zero or full16 BEFORE raw polynomial multiplication.
This is a positive submajorant, not a bound on the full source class.

At128 packets retain C[f,s], the coefficient of the raw matrix polynomial,
and E[f,s], the coefficient of P**127 T. For each actual regional profile
(128,0) or(127,2), divide by V=C(4096,f)C(4096-f,s) and take alpha. Both
remaining endpoint masks have multiplicity1. This is conditional Jensen
on the exact regional route degree, not a Cauchy evaluation. Assemble the
same256-region finite chains and both orientations; no H/u/w factors.

Compare with the pinned p128 evaluated-polynomial control at the same
parameters. A positive exact control obstructs only this selected upper
majorant; it is not evidence of a low-distance word. All paths, including
paths outside the quiet-count branch, remain with their Chernoff weights.

One48-second guarded invocation,20GiB free,128MiB estimated peak,
stdout only; no output files, caches, scans or parameter optimization.
"""

import hashlib
import importlib.util
import json
import math
import os
from pathlib import Path
import shutil
import sys
import time

sys.dont_write_bytecode = True
os.environ["OPENBLAS_NUM_THREADS"] = "1"
os.environ["OMP_NUM_THREADS"] = "1"
PIN = "2c306223f5b1f876670aede55ca80b905f16fb98364a840a645fea3b19cc32ad"


def main():
    if len(sys.argv) != 1:
        raise ValueError("fixed control accepts no scan arguments")
    if shutil.disk_usage(Path(__file__).resolve().anchor).free < 20*(1 << 30):
        raise RuntimeError("requires20GiB free on source volume")
    started = time.monotonic()
    import numpy as np
    from scipy.special import logsumexp
    path = Path(__file__).with_name("cheap_lane_pooling_scale_pilot.py")
    if hashlib.sha256(path.read_bytes()).hexdigest() != PIN:
        raise RuntimeError("frozen pooling helper hash mismatch")
    spec = importlib.util.spec_from_file_location("frozen_exact_region_pool", path)
    pool = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(pool)
    pool.np = np
    pool.lse = logsumexp
    pool.DEADLINE = started+48
    degree, pilot = pool.load()
    tests = degree.tests(pilot)
    singleton, filtered = pilot.imports()
    shapes, index, degrees, coeff, records, route_checks = pilot.geometry()
    m, t, refinements = pilot.raw_tables(shapes, index, records, degree.Z, singleton, filtered)
    a, b = degree.raw_polynomials(m, t, degrees, coeff, pilot)
    a = a[:, :, [0, 16], :][:, :, :, [0, 16]].copy()
    b = b[:, :, [0, 16], :].copy()
    del m, t
    for _ in range(7):
        b, a = degree.convolution(a, b), degree.convolution(a, a)
    old = pool.score(a, b, 128, "all_packet0_16", degree, pilot)
    assert abs(old["log2_value"]-13187.08) < .02
    regions, tails, diagnostics = [], [], []
    for f, s in ((128, 0), (127, 2)):
        logv = math.log(math.comb(4096, f)*math.comb(4096-f, s))
        reg = degree.ALPHA*(a[s, f]-logv)
        end = degree.ALPHA*(b[s, f]-logv)
        regions.append(reg)
        tails.append(end)
        # The one selected coefficient cannot exceed its Cauchy majorant.
        logu = math.log(degree.OLD[0 if s == 0 else 1])
        logw = -math.inf if s == 0 else math.log(degree.OLD[2])
        logU = math.log1p(math.exp(logu)+(math.exp(logw) if s else 0))
        h = 4096*logU-f*logu-(s*logw if s else 0)-logv
        evaluated = pool.tile_value(degree.normalize(a, 4096, degree.ALPHA, [1, 1]),
                                    4096, logu, logw, logU)+h
        evaluated_t = pool.tile_value(degree.normalize(b, 4096, degree.ALPHA, [1]),
                                      4096, logu, logw, logU)+h
        assert np.all(reg <= evaluated+1e-11)
        assert np.all(end <= evaluated_t+1e-11)
        diagnostics.append({"full": f, "single": s, "log2_route_count": logv/math.log(2),
            "regular_log2": (reg/math.log(2)).tolist(),
            "terminal_log2": (end/math.log(2)).tolist()})
    paths, orientations = pilot.chain(regions, tails)
    source = degree.ALPHA*(2*math.log(256/255)+256*(160*math.log(256)-32*math.log(255)))
    output = -degree.ALPHA*pilot.LIMIT*math.log(degree.Z)
    correction = -degree.ALPHA*24575*math.log(degree.QUIET)
    exact = (source+output+correction+math.log(4096)+paths)/math.log(2)
    assert exact <= old["log2_value"]+1e-8
    pool.guard()
    print(json.dumps({"scope": "FLOAT exact regional coefficient control; all packet masks restricted to0/full16",
        "whole_code_certificate": False, "full_source_class_bound": False, "directed": False,
        "parameters": {"z": degree.Z, "alpha": degree.ALPHA, "quiet_mark": degree.QUIET, "floor_mark": 1},
        "evaluated_region_log2": old["log2_value"], "exact_region_log2": exact,
        "gain_bits": old["log2_value"]-exact,
        "components_bits": {"source": source/math.log(2), "output_Chernoff": output/math.log(2),
            "quiet_correction": correction/math.log(2), "source_sets_per_orientation": 12,
            "finite_two_orientation_sum": paths/math.log(2)},
        "orientation_chain_bits": [x/math.log(2) for x in orientations], "regions": diagnostics,
        "tests": {**tests, "route_normalizers": route_checks, "coefficient_vs_Cauchy_entries": 12,
            "p128_control_reproduction": "passed", "regular_updates": 32767},
        "shared_point_refinements": refinements, "estimated_peak_bytes": 128*(1 << 20),
        "peak_working_set_bytes": pool.peak_working_set(), "generated_bytes": 0,
        "seconds": time.monotonic()-started, "source_sha256": hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
        "pooling_pin": PIN}, indent=2))


if __name__ == "__main__":
    main()
