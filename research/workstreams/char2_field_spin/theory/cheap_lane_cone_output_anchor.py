"""Three nearby output marks at the selected alpha, without a joint scan.

The fixed-cone exponent is119/1280=.09296875 from the prior bounded
stationarity test. Rebuild the actual raw caps separately at z=.84,.85,.86.
Keep all path marks1 and the same cone. Compare the two-state floor
operator at those three marks; then apply the newborn refinement only
at the best evaluated mark. Each source term uses its OWN z and alpha.
No numerical differentiation or extrapolated witness is used as a bound.
Only consecutive256-group sources, all masks. FLOAT, not a certificate.
One48-second job,20GiB free,128MiB estimated peak, stdout only.
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
from types import SimpleNamespace

sys.dont_write_bytecode = True
os.environ["OPENBLAS_NUM_THREADS"] = "1"
os.environ["OMP_NUM_THREADS"] = "1"
PIN = "ea5a4670f5135a690d8b8fd06a6bce2a4a4b94b8b6fda28b722a8762a6f27996"
ALPHA = 119/1280


def load_stationarity():
    path = Path(__file__).with_name("cheap_lane_cone_stationarity.py")
    if hashlib.sha256(path.read_bytes()).hexdigest() != PIN:
        raise RuntimeError("frozen stationarity helper hash mismatch")
    spec = importlib.util.spec_from_file_location("output_stationarity", path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def source_for_z(z):
    rate = (math.log2((256/255)**2)+256*(1280-32*math.log2(255))
            -(35*(1 << 21)//100)*math.log2(z))
    return SimpleNamespace(source_rate=lambda: rate)


def main():
    if len(sys.argv) != 1:
        raise ValueError("fixed output-mark check takes no scan arguments")
    if shutil.disk_usage(Path(__file__).resolve().anchor).free < 20*(1 << 30):
        raise RuntimeError("requires20GiB free")
    started = time.monotonic()
    import numpy as np
    from scipy.special import logsumexp
    stationarity = load_stationarity(); sensitivity = stationarity.load_alpha()
    newborn = sensitivity.load_pinned("cheap_lane_newborn_cone.py", stationarity.NEWBORN_PIN, "output_newborn")
    floor = sensitivity.load_pinned("cheap_lane_return_floor.py", sensitivity.FLOOR_PIN, "output_floor")
    disjoint = floor.load_disjoint(); cone = disjoint.load_cone()
    degree, pilot = cone.load_degree(np, logsumexp, started+48)
    singleton, filtered = pilot.imports()
    shapes, index, degrees, weights, records, route_checks = pilot.geometry()
    p = np.asarray([math.comb(16, j)*min(1., 256/255*255.**(-(16-j))) for j in range(1, 17)])
    tests = disjoint.exact_tests()
    rows = []; best = None; retained = None
    for z in (.85, .84, .86):
        m, t, _ = pilot.raw_tables(shapes, index, records, z, singleton, filtered)
        coarse, terminal = cone.cone_caps(m, t, p, np)
        cap, quiet, _, _, _ = floor.floor_caps(records, z, disjoint, np)
        assert np.allclose(quiet, m[:, 0, 0], rtol=3e-14, atol=0.)
        coarse[:, 1, 0] = np.minimum(coarse[:, 1, 0], cap)
        result = stationarity.score_at_alpha(coarse, terminal, degrees, weights, ALPHA,
                                             source_for_z(z), cone, degree, pilot)
        result["z"] = z; rows.append(result)
        if z == .85:
            assert abs(result["log2_bound"]-78.57605118331412) < 1e-6
        if best is None or result["log2_bound"] < best["log2_bound"]:
            best, retained = result, (m, t, coarse, terminal)
    m, t, coarse, terminal = retained
    types = np.asarray([sum(x > 0 for x in record[0]) for record in records])
    matrix, final = newborn.newborn_tables(m, t, p, types, (13, 9, 5, 1), coarse, terminal, np)
    refined = stationarity.score_at_alpha(matrix, final, degrees, weights, ALPHA,
                                          source_for_z(best["z"]), cone, degree, pilot)
    assert refined["source_bits_per_orientation"] == best["source_bits_per_orientation"]
    assert refined["log2_bound"] <= best["log2_bound"]+1e-6
    degree.guard()
    print(json.dumps({"scope": "FLOAT local output-mark comparison, consecutive256 supports only",
        "whole_code_certificate": False, "directed": False,
        "all_state_masks_covered": True, "all_source_supports_covered": False,
        "alpha_rational": "119/1280", "path_marks": 1, "two_state_comparisons": rows,
        "selected_z": best["z"], "same_witness_newborn": refined,
        "tests": {**tests, "exact_route_normalizers": route_checks,
            "old_anchor_reproduced": "passed", "matched_source_terms": "passed"},
        "seconds": time.monotonic()-started, "estimated_peak_bytes": 128*(1 << 20),
        "generated_bytes": 0, "source_sha256": hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
        "stationarity_pin": PIN}, indent=2))


if __name__ == "__main__":
    main()
