"""One positive-dual alpha derivative of the frozen two-state floor cone.

Keep z=.85, the same fixed cone, all path marks1, and the two actual
consecutive256-group profiles. For each occupancy record aggregate
weight*cap**alpha and weight*cap**alpha*(-log(cap)) BEFORE normalizing
each degree coefficient. The cost is nonnegative only because every
raw cap is<=1; fail rather than reuse the positive dual when that fails.
Zero-cap terms have zero mass and zero cost contribution.

The dual propagates total raw log-cost through32767 updates and the
honest output-only terminal. There is no second fractional power.
  d score/d alpha = source_rate_bits - E[total_cost]/log(2).
The4096 source multiplicity per orientation has derivative0.

Evaluate alpha=.05. ONLY if its derivative is negative, compare the one
predeclared alpha=.1 with the same raw caps. No scan or optimization.
FLOAT diagnostic, not directed or a whole-code certificate. One48-second
job,20GiB free,128MiB estimated peak, stdout only; no generated files.
"""

from fractions import Fraction
import hashlib
import importlib.util
import itertools
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
FLOOR_PIN = "2143482d17afe3624bd0f445773c561a47c60cf7ec45ff8bcef57ff3417b3588"
DUAL_PIN = "84b6b757a3fc4e85f28c0c3e9f59f6063ae8daf70cc3704aca6feb9ebe27600b"
Z, ALPHA, FORWARD_ALPHA = .85, .05, .1
LN2 = math.log(2)


def load_pinned(name, pin, label):
    path = Path(__file__).with_name(name)
    if hashlib.sha256(path.read_bytes()).hexdigest() != pin:
        raise RuntimeError("frozen helper hash mismatch: "+name)
    spec = importlib.util.spec_from_file_location(label, path)
    module = importlib.util.module_from_spec(spec); spec.loader.exec_module(module)
    return module


def cost_polys(m, t, degrees, weights, alpha, pilot):
    np = pilot.np
    assert np.isfinite(m).all() and np.isfinite(t).all()
    assert (m >= 0).all() and (t >= 0).all()
    if np.any(m > 1.) or np.any(t > 1.):
        raise RuntimeError("a raw cap exceeds1; nonnegative cost dual is not valid")
    index = (degrees[:, 1].astype(int), degrees[:, 0].astype(int))
    def aggregate(caps):
        shape = (3, 33, caps.shape[1], caps.shape[2])
        mass = np.zeros(shape); cost_mass = np.zeros(shape)
        cost = np.zeros_like(caps)
        np.log(caps, out=cost, where=caps > 0)
        cost = -cost
        local_mass = weights[:, None, None]*caps**alpha
        np.add.at(mass, index, local_mass)
        np.add.at(cost_mass, index, local_mass*cost)
        means = np.zeros(shape+(4,))
        np.divide(cost_mass, mass, out=means[:, :, :, :, 0], where=mass > 0)
        assert (means >= 0).all() and np.isfinite(means).all()
        return pilot.tolog(mass), means
    return aggregate(m), aggregate(t)


def exact_cost_tests(dual, pilot):
    np = pilot.np
    exponents = [[[1, 2], [3, 1]], [[2, 4], [1, 3]]]
    end_exponents = [[2, 3], [1, 4]]
    exact_m = [[[Fraction(1, 2**e) for e in row] for row in table] for table in exponents]
    exact_t = [[[Fraction(1, 2**e)] for e in table] for table in end_exponents]
    weights = (Fraction(1, 3), Fraction(2, 3))
    m, t = np.asarray(exact_m, dtype=float), np.asarray(exact_t, dtype=float)
    a, b = cost_polys(m, t, np.zeros((2, 2)), np.asarray(weights, dtype=float), 1., pilot)
    # Degree0 alone suffices to test record averaging and the positive
    # additive cost; trim irrelevant zero coefficients before convolution.
    a = a[0][:1, :1], a[1][:1, :1]
    b = b[0][:1, :1], b[1][:1, :1]
    product, terminal = dual.convolution(a, a, 0, 0), dual.convolution(a, b, 0, 0)
    checked = 0
    for is_terminal, result in ((False, product), (True, terminal)):
        for i in range(2):
            for j in range(1 if is_terminal else 2):
                total = cost = Fraction(0)
                for r0, r1, mid in itertools.product(range(2), repeat=3):
                    right = exact_t[r1][mid][0] if is_terminal else exact_m[r1][mid][j]
                    exponent = end_exponents[r1][mid] if is_terminal else exponents[r1][mid][j]
                    mass = weights[r0]*weights[r1]*exact_m[r0][i][mid]*right
                    total += mass
                    cost += mass*(exponents[r0][i][mid]+exponent)
                assert abs(result[0][0, 0, i, j]-math.log(total)) < 3e-13
                assert abs(result[1][0, 0, i, j, 0]-float(cost/total)*LN2) < 3e-13
                assert np.all(result[1][0, 0, i, j, 1:] == 0.)
                checked += 1
    zero = np.zeros((1, 2, 2)); end_zero = np.zeros((1, 2, 1))
    za, zb = cost_polys(zero, end_zero, np.zeros((1, 2)), np.ones(1), .05, pilot)
    assert not np.isfinite(za[0]).any() and not np.isfinite(zb[0]).any()
    assert (za[1] == 0).all() and (zb[1] == 0).all()
    failed = False
    try:
        cost_polys(np.full((1, 2, 2), 1.01), end_zero, np.zeros((1, 2)), np.ones(1), .05, pilot)
    except RuntimeError:
        failed = True
    assert failed
    return {"exact_rational_record_and_path_costs": checked,
            "zero_mass_cost_terms": "passed", "cap_above1_rejected": "passed"}


def source_rate():
    limit = 35*(1 << 21)//100
    return math.log2((256/255)**2)+256*(1280-32*math.log2(255))-limit*math.log2(Z)


def evaluate_cost(m, t, degrees, weights, alpha, dual, cone, degree, pilot):
    np = pilot.np
    a, b = cost_polys(m, t, degrees, weights, alpha, pilot)
    plain_a, plain_b = cone.route_polys(m, t, degrees, weights, alpha, pilot)
    assert np.allclose(a[0], plain_a, atol=2e-13)
    assert np.allclose(b[0], plain_b, atol=2e-13)
    for _ in range(7):
        b, a = dual.convolution(a, b), dual.convolution(a, a)
        plain_b, plain_a = degree.convolution(plain_a, plain_b), degree.convolution(plain_a, plain_a)
    regions, tails, plain_regions, plain_tails = [], [], [], []
    for f, s in ((128, 0), (127, 2)):
        norm = math.log(math.comb(4096, f)*math.comb(4096-f, s))
        regions.append((a[0][s, f]-norm, a[1][s, f]))
        tails.append((b[0][s, f]-norm, b[1][s, f]))
        plain_regions.append(plain_a[s, f]-norm); plain_tails.append(plain_b[s, f]-norm)
    paths, means, _ = dual.chain(regions, tails)
    plain_paths, _ = pilot.chain(plain_regions, plain_tails)
    assert abs(paths-plain_paths) < 1e-7
    expected_cost = float(means[0])
    assert np.all(means[1:] == 0.)
    score = 12+alpha*source_rate()+paths/LN2
    slope = source_rate()-expected_cost/LN2
    return {"alpha": alpha, "log2_bound": score,
            "source_bits_per_orientation": 12+alpha*source_rate(),
            "source_derivative_bits": source_rate(),
            "expected_raw_log_cost_nats": expected_cost,
            "path_derivative_bits": -expected_cost/LN2,
            "score_derivative_bits_per_alpha": slope,
            "scalar_dual_agreement": "passed"}


def main():
    if len(sys.argv) != 1:
        raise ValueError("fixed alpha diagnostic takes no scan arguments")
    if shutil.disk_usage(Path(__file__).resolve().anchor).free < 20*(1 << 30):
        raise RuntimeError("requires20GiB free")
    started = time.monotonic(); deadline = started+48
    import numpy as np
    from scipy.special import logsumexp
    floor = load_pinned("cheap_lane_return_floor.py", FLOOR_PIN, "alpha_floor")
    disjoint = floor.load_disjoint(); cone = disjoint.load_cone()
    dual = load_pinned("cheap_lane_branch_sensitivity.py", DUAL_PIN, "alpha_positive_dual")
    dual.np, dual.lse, dual.DEADLINE = np, logsumexp, deadline
    degree, pilot = cone.load_degree(np, logsumexp, deadline)
    tests = exact_cost_tests(dual, pilot); tests.update(dual.tests(pilot))
    singleton, filtered = pilot.imports()
    shapes, index, degrees, weights, records, route_checks = pilot.geometry()
    m, t, refinements = pilot.raw_tables(shapes, index, records, Z, singleton, filtered)
    p = np.asarray([math.comb(16, j)*min(1., 256/255*255.**(-(16-j))) for j in range(1, 17)])
    coarse, terminal = cone.cone_caps(m, t, p, np)
    old_cap, _, _ = disjoint.return_caps(records, Z, np)
    cap, quiet, _, _, _ = floor.floor_caps(records, Z, disjoint, np)
    assert np.allclose(quiet, m[:, 0, 0], rtol=3e-14, atol=0.)
    coarse[:, 1, 0] = np.minimum(coarse[:, 1, 0], np.minimum(old_cap, cap))
    current = evaluate_cost(coarse, terminal, degrees, weights, ALPHA, dual, cone, degree, pilot)
    assert abs(current["log2_bound"]-3262.2132476853913) < 1e-6
    forward = None
    if current["score_derivative_bits_per_alpha"] < 0:
        forward = evaluate_cost(coarse, terminal, degrees, weights, FORWARD_ALPHA, dual, cone, degree, pilot)
        # The positive sum of exponentials has a convex log partition.
        tangent = current["log2_bound"]+(FORWARD_ALPHA-ALPHA)*current["score_derivative_bits_per_alpha"]
        assert forward["log2_bound"] >= tangent-1e-6
    dual.guard(); degree.guard()
    print(json.dumps({"scope": "FLOAT alpha sensitivity of the two-state return-floor cone only",
        "whole_code_certificate": False, "directed": False, "all_state_masks_covered": True,
        "all_source_supports_covered": False, "parameters": {"z": Z, "path_marks": 1},
        "current": current, "predeclared_forward_comparison": forward,
        "forward_rule": "alpha=.1 only if the derivative at alpha=.05 is negative",
        "raw_cap_maxima": {"regular": float(np.max(coarse)), "terminal": float(np.max(terminal))},
        "tests": {**tests, "exact_route_normalizers": route_checks,
            "frozen_floor_control": "passed", "same_degree_mass_before_dual": "passed"},
        "pair_only_refinements": refinements, "seconds": time.monotonic()-started,
        "estimated_peak_bytes": 128*(1 << 20), "generated_bytes": 0,
        "source_sha256": hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
        "floor_pin": FLOOR_PIN, "positive_dual_pin": DUAL_PIN}, indent=2))


if __name__ == "__main__":
    main()
