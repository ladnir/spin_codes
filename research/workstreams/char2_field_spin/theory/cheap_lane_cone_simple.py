"""Test a small independently valid two-state cone at one fixed witness.

Use z=21/25,alpha=119/1280 and the frozen positive cone p. Keep exact
quiet Q. Activation A is the minimum over gamma0,1/2,1 of
 z^gamma*(g(z^(1-gamma))^ell-((1-z^(1-gamma))/256)^ell),
with A=0 at ell0. Let T_h be the minimum gamma output-only cap. Set
 U=sum_h p_h*T_h and L=U by the joint-live bound M_hj<=p_j*T_h.
Use only the positive return-floor cap kappa*q^-16*R for E, dropping
the old cardinality-return minimum. Terminal zero-state output is g^ell.

The simplified operator is independently valid and is compared pointwise
against the pinned full17-state floor-cone construction. Its floating
implementation still reads the output-only final table from that helper;
a directed verifier need only reproduce that smaller output interface.
No alpha or z scan, no path marks, no frozen edits. Consecutive256-group
sources only; not a directed or whole-code certificate. One48-second job,
20GiB free,128MiB estimated peak, stdout only, no generated files.
"""

from fractions import Fraction
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
PIN = "4bac712a2933825de4ade73ad2db3701326894988749d6fc7a79fcfd27c73356"
Z, ALPHA = 21/25, 119/1280


def load_anchor():
    path = Path(__file__).with_name("cheap_lane_cone_output_anchor.py")
    if hashlib.sha256(path.read_bytes()).hexdigest() != PIN:
        raise RuntimeError("frozen output-anchor hash mismatch")
    spec = importlib.util.spec_from_file_location("simple_output_anchor", path)
    module = importlib.util.module_from_spec(spec); spec.loader.exec_module(module)
    return module


def activation_tests():
    # Independent positive subset expansion for the activation expression;
    # remove only the all-requested-zero subset, not the scalarzero atom.
    checked = 0
    for b in (Fraction(0), Fraction(1, 3), Fraction(2, 3), Fraction(1)):
        g = b+(1-b)/256
        for ell in range(9):
            direct = sum(Fraction(math.comb(ell, k))*b**k*((1-b)/256)**(ell-k)
                         for k in range(1, ell+1))
            assert direct == g**ell-((1-b)/256)**ell
            assert 0 <= direct <= g**ell
            checked += 1
    return {"exact_activation_subset_expansions": checked}


def simple_tables(shapes, index, records, p, filtered, floor, disjoint, np, guard):
    output = None
    for gamma in (0., .5, 1.):
        guard()
        _, final = filtered.local_tables(shapes, Z, 1., gamma, 1.)
        current = final[:, 0]+final[:, 2]
        output = current if output is None else np.minimum(output, current)
    output = output[index]
    ell = np.asarray([sum(record[0]) for record in records], dtype=int)
    activation = np.full(len(records), math.inf)
    for gamma in (0., .5, 1.):
        b = Z**(1-gamma); g = b+(1-b)/256
        candidate = Z**gamma*(np.power(g, ell)-np.power((1-b)/256, ell))
        activation = np.minimum(activation, candidate)
    activation[ell == 0] = 0.
    u = output[:, 1:, 0]@p
    returns, quiet, kappa, _, _ = floor.floor_caps(records, Z, disjoint, np)
    matrix = np.empty((len(records), 2, 2)); terminal = np.empty((len(records), 2, 1))
    matrix[:, 0, 0] = quiet; matrix[:, 0, 1] = activation
    matrix[:, 1, 0] = returns; matrix[:, 1, 1] = u
    terminal[:, 0, 0] = np.power(Z+(1-Z)/256, ell)
    terminal[:, 1, 0] = u
    assert np.isfinite(matrix).all() and (matrix >= 0).all()
    assert np.isfinite(terminal).all() and (terminal >= 0).all()
    return matrix, terminal, kappa


def main():
    if len(sys.argv) != 1:
        raise ValueError("fixed simplified-cone diagnostic takes no scan arguments")
    if shutil.disk_usage(Path(__file__).resolve().anchor).free < 20*(1 << 30):
        raise RuntimeError("requires20GiB free")
    started = time.monotonic()
    import numpy as np
    from scipy.special import logsumexp
    anchor = load_anchor(); stationarity = anchor.load_stationarity()
    sensitivity = stationarity.load_alpha()
    floor = sensitivity.load_pinned("cheap_lane_return_floor.py", sensitivity.FLOOR_PIN, "simple_floor")
    disjoint = floor.load_disjoint(); cone = disjoint.load_cone()
    degree, pilot = cone.load_degree(np, logsumexp, started+48)
    tests = activation_tests(); tests.update(cone.exact_tests()); tests.update(floor.exact_tests())
    singleton, filtered = pilot.imports()
    shapes, index, degrees, weights, records, route_checks = pilot.geometry()
    p = np.asarray([math.comb(16, j)*min(1., 256/255*255.**(-(16-j))) for j in range(1, 17)])
    m, t, refinements = pilot.raw_tables(shapes, index, records, Z, singleton, filtered)
    old_matrix, old_terminal = cone.cone_caps(m, t, p, np)
    cap, quiet, _, _, _ = floor.floor_caps(records, Z, disjoint, np)
    old_matrix[:, 1, 0] = np.minimum(old_matrix[:, 1, 0], cap)
    old = stationarity.score_at_alpha(old_matrix, old_terminal, degrees, weights, ALPHA,
                                     anchor.source_for_z(Z), cone, degree, pilot)
    assert abs(old["log2_bound"]+218.6695) < .002
    matrix, terminal, kappa = simple_tables(shapes, index, records, p, filtered,
                                          floor, disjoint, np, degree.guard)
    assert np.allclose(matrix[:, 0, 0], old_matrix[:, 0, 0], rtol=3e-14, atol=0.)
    assert np.all(matrix >= old_matrix*(1-3e-13))
    assert np.all(terminal >= old_terminal*(1-3e-13))
    assert np.allclose(terminal, old_terminal, rtol=3e-14, atol=0.)
    result = stationarity.score_at_alpha(matrix, terminal, degrees, weights, ALPHA,
                                        anchor.source_for_z(Z), cone, degree, pilot)
    assert result["source_bits_per_orientation"] == old["source_bits_per_orientation"]
    assert result["log2_bound"] >= old["log2_bound"]-1e-6
    changes = {}
    for name, i, j in (("quiet", 0, 0), ("activation", 0, 1), ("return", 1, 0), ("live", 1, 1)):
        positive = old_matrix[:, i, j] > 0
        changes[name] = {"strictly_enlarged_records": int(np.count_nonzero(matrix[:, i, j] > old_matrix[:, i, j]*(1+1e-12))),
                         "newly_positive_records": int(np.count_nonzero((~positive) & (matrix[:, i, j] > 0))),
                         "max_ratio_on_old_positive": float(np.max(matrix[positive, i, j]/old_matrix[positive, i, j]))}
    degree.guard()
    print(json.dumps({"scope": "FLOAT simplification of two-state floor cone at one fixed witness",
        "whole_code_certificate": False, "directed": False,
        "all_state_masks_covered": True, "all_source_supports_covered": False,
        "parameters": {"z_rational": "21/25", "alpha_rational": "119/1280", "path_marks": 1},
        "old_full_table_control": old, "simplified_cone": result,
        "simplification_cost_bits": result["log2_bound"]-old["log2_bound"],
        "pointwise_changes": changes, "kappa_exact": str(kappa),
        "tests": {**tests, "exact_route_normalizers": route_checks,
            "pointwise_majorization": "passed", "quiet_and_terminal_equal": "passed",
            "frozen_output_control": "passed", "matched_source_terms": "passed"},
        "pair_only_refinements": refinements, "seconds": time.monotonic()-started,
        "estimated_peak_bytes": 128*(1 << 20), "generated_bytes": 0,
        "source_sha256": hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
        "output_anchor_pin": PIN}, indent=2))


if __name__ == "__main__":
    main()
