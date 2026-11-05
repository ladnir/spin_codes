"""One bounded alpha stationarity check, followed by a same-alpha A/B.

Keep z=.85, the frozen two-state return-floor cone, neutral path marks,
and actual consecutive256-group profiles. Reproduce the two endpoints
alpha=.05,.1 and bisect their derivative sign bracket at most six times.
The best sampled SCORE, not an interpolated root, is the selected witness.
At that SAME alpha evaluate the frozen six-state newborn refinement.

The source term is12+alpha*source_rate for EVERY evaluation. Do not call
the frozen cone score wrapper, whose source exponent is fixed at.05.
No z scan, anchor scan, unbounded optimizer, or parameter extrapolation.
This is a FLOAT diagnostic, not a directed root/certificate. One48-second
invocation,20GiB free,128MiB estimated peak, stdout only; no generated data.
"""

import hashlib
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
ALPHA_PIN = "c0b5976135dc91f317209b9618abb5266ad6d626b74feb89bce7a7a929761700"
NEWBORN_PIN = "0d4b1d7f403b9bde6a55aaf1d81c502995f69cfd4eb4ab341607c27e5637641a"
LN2 = math.log(2)


def load_alpha():
    import importlib.util
    path = Path(__file__).with_name("cheap_lane_cone_alpha.py")
    if hashlib.sha256(path.read_bytes()).hexdigest() != ALPHA_PIN:
        raise RuntimeError("frozen alpha helper hash mismatch")
    spec = importlib.util.spec_from_file_location("stationarity_alpha", path)
    module = importlib.util.module_from_spec(spec); spec.loader.exec_module(module)
    return module


def score_at_alpha(m, t, degrees, weights, alpha, sensitivity, cone, degree, pilot):
    """Same finite route-first chain, with the supplied source exponent."""
    a, b = cone.route_polys(m, t, degrees, weights, alpha, pilot)
    for _ in range(7):
        b, a = degree.convolution(a, b), degree.convolution(a, a)
    regions, tails = [], []
    for f, s in ((128, 0), (127, 2)):
        norm = math.log(math.comb(4096, f)*math.comb(4096-f, s))
        regions.append(a[s, f]-norm); tails.append(b[s, f]-norm)
    paths, orientations = pilot.chain(regions, tails)
    source = 12+alpha*sensitivity.source_rate()
    return {"alpha": alpha, "log2_bound": source+paths/LN2,
            "source_bits_per_orientation": source,
            "orientation_chain_bits": [value/LN2 for value in orientations]}


def main():
    if len(sys.argv) != 1:
        raise ValueError("fixed bracket diagnostic takes no scan arguments")
    if shutil.disk_usage(Path(__file__).resolve().anchor).free < 20*(1 << 30):
        raise RuntimeError("requires20GiB free")
    started = time.monotonic(); deadline = started+48
    import numpy as np
    from scipy.special import logsumexp
    sensitivity = load_alpha()
    newborn = sensitivity.load_pinned("cheap_lane_newborn_cone.py", NEWBORN_PIN, "stationarity_newborn")
    floor = sensitivity.load_pinned("cheap_lane_return_floor.py", sensitivity.FLOOR_PIN, "stationarity_floor")
    disjoint = floor.load_disjoint(); cone = disjoint.load_cone()
    dual = sensitivity.load_pinned("cheap_lane_branch_sensitivity.py", sensitivity.DUAL_PIN, "stationarity_dual")
    dual.np, dual.lse, dual.DEADLINE = np, logsumexp, deadline
    degree, pilot = cone.load_degree(np, logsumexp, deadline)
    tests = sensitivity.exact_cost_tests(dual, pilot)
    tests.update(newborn.exact_tests(cone))
    singleton, filtered = pilot.imports()
    shapes, index, degrees, weights, records, route_checks = pilot.geometry()
    m, t, refinements = pilot.raw_tables(shapes, index, records, sensitivity.Z, singleton, filtered)
    p = np.asarray([math.comb(16, j)*min(1., 256/255*255.**(-(16-j))) for j in range(1, 17)])
    coarse, terminal = cone.cone_caps(m, t, p, np)
    old_cap, _, _ = disjoint.return_caps(records, sensitivity.Z, np)
    cap, quiet, _, _, _ = floor.floor_caps(records, sensitivity.Z, disjoint, np)
    coarse[:, 1, 0] = np.minimum(coarse[:, 1, 0], np.minimum(old_cap, cap))
    assert np.allclose(quiet, m[:, 0, 0], rtol=3e-14, atol=0.)
    evaluations = []
    def evaluate(alpha, label):
        value = sensitivity.evaluate_cost(coarse, terminal, degrees, weights,
                                          alpha, dual, cone, degree, pilot)
        value["label"] = label; evaluations.append(value)
        print(json.dumps({"progress": "alpha_evaluated", "alpha": alpha,
            "log2_bound": value["log2_bound"],
            "derivative_bits": value["score_derivative_bits_per_alpha"],
            "elapsed_seconds": time.monotonic()-started}), flush=True)
        return value
    lower = evaluate(.05, "bracket_lower")
    upper = evaluate(.1, "bracket_upper")
    assert abs(lower["log2_bound"]-3262.2132476853913) < 1e-6
    assert abs(upper["log2_bound"]-137.1162226508386) < 1e-6
    assert abs(lower["score_derivative_bits_per_alpha"]+172170.314433482) < 1e-3
    assert abs(upper["score_derivative_bits_per_alpha"]-16628.654825883452) < 1e-3
    assert lower["score_derivative_bits_per_alpha"] < 0 < upper["score_derivative_bits_per_alpha"]
    reserve_stop = False
    for step in range(6):
        # Reserve bounded time for the single six-state comparison.
        if deadline-time.monotonic() < 12:
            reserve_stop = True
            break
        current = evaluate((lower["alpha"]+upper["alpha"])/2, "bisect_"+str(step+1))
        slope = current["score_derivative_bits_per_alpha"]
        assert lower["score_derivative_bits_per_alpha"]-1e-3 <= slope <= upper["score_derivative_bits_per_alpha"]+1e-3
        if slope <= 0:
            lower = current
        else:
            upper = current
    assert len(evaluations) <= 8
    best = min(evaluations, key=lambda value: value["log2_bound"])
    selected = best["alpha"]
    types = np.asarray([sum(x > 0 for x in record[0]) for record in records])
    matrix, final = newborn.newborn_tables(m, t, p, types, (13, 9, 5, 1), coarse, terminal, np)
    for birth in range(1, 5):
        assert np.all(matrix[types+birth <= 4, birth, 0] == 0.)
    refined = score_at_alpha(matrix, final, degrees, weights, selected,
                             sensitivity, cone, degree, pilot)
    assert refined["source_bits_per_orientation"] == best["source_bits_per_orientation"]
    assert refined["log2_bound"] <= best["log2_bound"]+1e-6
    degree.guard(); dual.guard()
    print(json.dumps({"scope": "FLOAT bounded alpha root bracket and same-alpha newborn A/B, consecutive256 supports",
        "whole_code_certificate": False, "directed": False, "all_state_masks_covered": True,
        "all_source_supports_covered": False, "parameters": {"z": sensitivity.Z, "path_marks": 1},
        "derivative_evaluations": evaluations, "bisection_steps": len(evaluations)-2,
        "stopped_to_reserve_final_comparison": reserve_stop,
        "final_float_bracket": {"lower_alpha": lower["alpha"], "upper_alpha": upper["alpha"],
            "lower_derivative": lower["score_derivative_bits_per_alpha"],
            "upper_derivative": upper["score_derivative_bits_per_alpha"]},
        "selected_two_state_witness": best, "same_alpha_newborn": refined,
        "newborn_improvement_bits": best["log2_bound"]-refined["log2_bound"],
        "tests": {**tests, "exact_route_normalizers": route_checks, "frozen_endpoint_replays": 2,
            "monotone_derivative_checks": len(evaluations)-2,
            "same_alpha_source_term": "passed", "short_return_support_exclusions": 4},
        "pair_only_refinements": refinements, "seconds": time.monotonic()-started,
        "estimated_peak_bytes": 128*(1 << 20), "generated_bytes": 0,
        "source_sha256": hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
        "alpha_pin": ALPHA_PIN, "newborn_pin": NEWBORN_PIN}, indent=2))


if __name__ == "__main__":
    main()
