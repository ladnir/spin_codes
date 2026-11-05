"""Bounded marker calibration for one r512 Cauchy profile, with Arb replay.

Only the two positive route markers vary. Keep r512, actual profiles
(16,480) in both matchings, z7/10 and alpha119/1280. Cache the frozen
local powered caps. A FLOAT proposal uses positive finite-chain moments,
not a Perron surrogate: H=e0 P^32767 T, with the output-only terminal.
At most12 evaluations apply safeguarded Newton steps in logarithmic
markers; first and second derivatives are count means and covariances.
Failure to reach stationarity does not trigger any extra search.

Round the best proposal to rational markers of denominator at most2^20.
Independently reconstruct its exact rational route weights and evaluate
the full positive bound in directed Arb arithmetic. The FLOAT objective
and derivatives are never used as a certificate. Source Ω, conservative
run-total source count, alpha placement and terminal remain unchanged.
The final bound covers only this selected actual profile, not all r512
supports or the whole code.60s/128MiB/20GiB, stdout only, no output files.
"""

import argparse
from fractions import Fraction as F
import hashlib
import importlib.util
from itertools import combinations_with_replacement, product
import json
import math
from math import comb
from pathlib import Path
import shutil
import sys
import time

sys.dont_write_bytecode = True
GRADIENT_PIN = "2bdf98517687e3fa549f412bab29db45752e3b4e6b5e8984f50838c9916f5631"
MAX_EVALUATIONS = 12


def load_gradient():
    path = Path(__file__).with_name("cheap_lane_cauchy_marker_gradient.py")
    if hashlib.sha256(path.read_bytes()).hexdigest() != GRADIENT_PIN:
        raise RuntimeError("frozen marker-gradient hash mismatch")
    spec = importlib.util.spec_from_file_location("marker_calibration_gradient", path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def float_coefficients(powered, np, guard):
    """Proposal-only full local polynomial, not a regional coefficient array."""
    half = []
    for f in range(17):
        for s in range(17-f):
            for left in range(s+1):
                a, b = f+left, f+s-left
                if a <= b:
                    weight = (comb(16, f)*comb(16-f, s)*comb(s, left)
                              *(1 if a == b else 2)/2**s)
                    half.append((f, s, a, b, weight))
    coefficients = np.zeros((33, 33, 6), dtype=float)
    masses = np.zeros((33, 33), dtype=float)
    approximate = {occ: np.asarray([float(x.mid()) for x in values])
                   for occ, values in powered.items()}
    for i, (f, s, a, b, weight) in enumerate(half):
        guard()
        for j in range(i, len(half)):
            ff, ss, c, d, other = half[j]
            mass = weight*other*(1 if i == j else 2)
            coefficients[f+ff, s+ss] += mass*approximate[tuple(sorted((a, b, c, d)))]
            masses[f+ff, s+ss] += mass
    for f in range(33):
        for s in range(33-f):
            exact = comb(32, f)*comb(32-f, s)
            assert abs(masses[f, s]-exact) <= 2e-13*max(1, exact)
    return coefficients


def normalized_product(left, right, np):
    """Six positive moments: mass,F,S,F²,FS,S²; scales are not differentiated."""
    a, la = left; b, lb = right
    c = np.empty((6, a.shape[1], b.shape[2]), dtype=float)
    c[0] = a[0]@b[0]
    c[1] = a[1]@b[0]+a[0]@b[1]
    c[2] = a[2]@b[0]+a[0]@b[2]
    c[3] = a[3]@b[0]+2*(a[1]@b[1])+a[0]@b[3]
    c[4] = a[4]@b[0]+a[1]@b[2]+a[2]@b[1]+a[0]@b[4]
    c[5] = a[5]@b[0]+2*(a[2]@b[2])+a[0]@b[5]
    scale = float(np.max(c[0]))
    if not scale > 0 or not np.all(np.isfinite(c)):
        raise ArithmeticError("nonfinite FLOAT proposal moment")
    return c/scale, la+lb+math.log(scale)


def normalized_power(value, exponent, np):
    identity = np.zeros((6, 2, 2), dtype=float)
    identity[0] = np.eye(2)
    result = identity, 0.0
    while exponent:
        if exponent & 1:
            result = normalized_product(result, value, np)
        exponent //= 2
        if exponent:
            value = normalized_product(value, value, np)
    return result


def proposal_tests(np):
    """Exact short path enumeration checks normalized six-moment products."""
    entries = [[[(0, 0, F(1, 3)), (1, 0, F(2, 5))], [(0, 1, F(3, 7))]],
               [[(1, 0, F(5, 11))], [(0, 0, F(7, 13)), (0, 1, F(11, 17))]]]
    terminal = [[(1, 0, F(13, 19))], [(0, 1, F(17, 23)), (0, 0, F(19, 29))]]
    u, v = F(2, 7), F(3, 11)

    def moment(terms):
        result = [F(0)]*6
        for f, s, weight in terms:
            mass = weight*u**f*v**s
            for d, count in enumerate((1, f, s, f*f, f*s, s*s)):
                result[d] += count*mass
        return result

    p = np.zeros((6, 2, 2)); t = np.zeros((6, 2, 1))
    for i in range(2):
        t[:, i, 0] = [float(x) for x in moment(terminal[i])]
        for j in range(2):
            p[:, i, j] = [float(x) for x in moment(entries[i][j])]
    actual, scale = normalized_product(normalized_power((p, 0.0), 3, np), (t, 0.0), np)
    checks = 0
    for start in range(2):
        terms = []
        for states in product((0, 1), repeat=3):
            route = (start,)+states
            choices = [entries[route[k]][route[k+1]] for k in range(3)]
            choices.append(terminal[route[-1]])
            for selected in product(*choices):
                f = sum(x[0] for x in selected); s = sum(x[1] for x in selected)
                weight = F(1)
                for x in selected:
                    weight *= x[2]
                terms.append((f, s, weight))
        expected = moment(terms)
        for d in range(6):
            assert abs(actual[d, start, 0]*math.exp(scale)-float(expected[d])) <= 1e-12*max(1, float(expected[d]))
            checks += 1
    return checks


def calibrate(coefficients, constant, targets, start, np, guard):
    """At most12 proposals; no evaluation follows a budget stop."""
    full = np.arange(33, dtype=float)[:, None]
    single = np.arange(33, dtype=float)[None, :]
    counts = [np.ones((33, 33)), full+np.zeros((33, 33)),
              single+np.zeros((33, 33)), full*full+np.zeros((33, 33)),
              full*single, single*single+np.zeros((33, 33))]
    history = []; best = None

    def evaluate(x):
        nonlocal best
        guard()
        if len(history) >= MAX_EVALUATIONS:
            raise RuntimeError("internal proposal-budget overrun")
        monomials = np.exp(full*x[0]+single*x[1])
        values = np.asarray([np.sum(coefficients*(monomials*count)[:, :, None], axis=(0, 1))
                             for count in counts])
        p = values[:, (0, 1, 2, 3)].reshape((6, 2, 2))
        t = values[:, (4, 5)].reshape((6, 2, 1))
        result, scale = normalized_product(normalized_power((p, 0.0), 32767, np), (t, 0.0), np)
        moments = result[:, 0, 0]/result[0, 0, 0]
        means = moments[1:3]
        covariance = np.asarray([[moments[3]-means[0]**2, moments[4]-means[0]*means[1]],
                                 [moments[4]-means[0]*means[1], moments[5]-means[1]**2]])
        objective = scale+math.log(result[0, 0, 0])-float(targets@x)+constant
        point = {"x": x.copy(), "objective": objective, "gradient": means-targets,
                 "covariance": covariance, "means": means}
        history.append({"evaluation": len(history)+1, "log2_proposal_bound": objective/math.log(2),
                        "markers": [math.exp(float(t)) for t in x],
                        "regional_means": [float(t/256) for t in means]})
        if best is None or objective < best["objective"]:
            best = point
        return point

    current = evaluate(np.log(np.asarray(start, dtype=float)))
    status = "evaluation budget exhausted"
    while len(history) < MAX_EVALUATIONS:
        if np.max(np.abs(current["gradient"]/256)) < 1e-5:
            status = "FLOAT first-order tolerance reached"
            break
        covariance = current["covariance"]
        if np.all(np.isfinite(covariance)) and np.linalg.eigvalsh(covariance)[0] > 1e-8:
            direction = -np.linalg.solve(covariance, current["gradient"])
        else:
            direction = -current["gradient"]/np.maximum(targets, 1)
        if float(direction@current["gradient"]) >= 0:
            direction = -current["gradient"]/np.maximum(targets, 1)
        direction /= max(1.0, float(np.max(np.abs(direction)))/2)
        accepted = False
        step = 1.0
        while len(history) < MAX_EVALUATIONS:
            trial = np.clip(current["x"]+step*direction, -20.0, 4.0)
            displacement = trial-current["x"]
            if float(np.max(np.abs(displacement))) < 1e-12:
                status = "proposal box or negligible step reached"
                break
            candidate = evaluate(trial)
            if candidate["objective"] <= current["objective"]+1e-4*float(current["gradient"]@displacement):
                current = candidate; accepted = True
                break
            step *= 0.5
        if not accepted:
            break
    assert best is not None and len(history) <= MAX_EVALUATIONS
    return best, history, status


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--precision", type=int, choices=(256, 384, 512), default=256)
    args = parser.parse_args()
    if shutil.disk_usage(Path(__file__).resolve().anchor).free < 20*(1 << 30):
        raise RuntimeError("requires20GiB free on source volume")
    started = time.monotonic(); deadline = started+60
    gradient = load_gradient(); anchor = gradient.load_anchor()
    peak_memory = anchor.peak_memory_reader()

    def guard():
        if time.monotonic() > deadline:
            raise TimeoutError("60-second calibration budget exceeded; no automatic retry")
        if peak_memory() > anchor.MEMORY_LIMIT:
            raise MemoryError("128MiB process-memory budget exceeded")

    import numpy as np
    from flint import arb, ctx
    ctx.prec = args.precision
    checks = proposal_tests(np)+gradient.exact_tests()
    adapter = anchor.load_adapter(); core = adapter.load_verifier()

    class CalibrationVerifier(core.Verifier):
        def guard(self):
            guard()

    verifier = CalibrationVerifier(2, deadline)
    r, f, s, old_u, old_v, z = anchor.ANCHORS[0]
    core.Z, core.ALPHA = z, anchor.ALPHA
    shapes = list(combinations_with_replacement(range(17), 4))
    caps = verifier.local_caps(shapes)
    alpha = verifier.exact(anchor.ALPHA)
    powered = {occ: tuple(verifier.positive(x)**alpha for x in values) for occ, values in caps.items()}
    del caps
    h, c, eta = anchor.source_factor(r)
    source = (verifier.exact(eta)*verifier.exact(core.BETA)**r*verifier.exact(z)**(-core.LIMIT))**alpha
    number = anchor.run_total(4096, r, s)
    normalizer = comb(4096, f)*comb(4096-f, s)
    constant_ball = (number*source).log()-256*arb(normalizer).log()
    constant = float(constant_ball.mid())
    coefficients = float_coefficients(powered, np, guard)
    best, history, status = calibrate(coefficients, constant, np.asarray([256*f, 256*s], dtype=float),
                                     [float(old_u), float(old_v)], np, guard)
    chosen = tuple(F(math.exp(float(x))).limit_denominator(1 << 20) for x in best["x"])
    # Rational rounding is not part of the proof until this independent replay.
    if any(x <= 0 for x in chosen):
        raise ArithmeticError("nonpositive rounded marker")

    def directed_receipt(u, v):
        moments = gradient.shape_moments(u, v, guard)
        entries = [[arb(0), arb(0), arb(0)] for _ in range(6)]
        for occ, weights in moments.items():
            guard()
            masses = [verifier.exact(x) for x in weights]
            for j, cap in enumerate(powered[occ]):
                for axis in range(3):
                    entries[j][axis] += masses[axis]*cap
        p = [[tuple(entries[0]), tuple(entries[1])], [tuple(entries[2]), tuple(entries[3])]]
        t = [[tuple(entries[4])], [tuple(entries[5])]]
        result = gradient.dual_matmul(gradient.dual_matpow(p, 32767, arb(0), arb(1)), t, arb(0))[0][0]
        denominator = verifier.exact(u)**f*verifier.exact(v)**s*normalizer
        bound = number*source*result[0]/denominator**256
        logtwo = arb(2).log()
        return bound, {"markers": [str(u), str(v)], "log2_bound_interval": str(bound.log()/logtwo),
            "regional_mean_intervals": [str(result[1]/result[0]/256), str(result[2]/result[0]/256)],
            "log_marker_gradient_intervals": [str((result[1]/result[0]-256*f)/logtwo),
                                               str((result[2]/result[0]-256*s)/logtwo)],
            "selected_profile_proved_below_2_pow_minus_40": bool(bound < arb(2)**(-40))}

    old_bound, old_receipt = directed_receipt(old_u, old_v)
    new_bound, new_receipt = directed_receipt(*chosen)
    guard()
    print(json.dumps({"scope": "one selected actual r512 matching profile only",
        "directed_final_replay": True, "whole_code_certificate": False,
        "all_profiles_at_weight_covered": False,
        "parameters": {"source_weight": r, "ordered_profile": [[f, s], [f, s]],
                       "z": str(z), "alpha": str(anchor.ALPHA), "updates": 32767},
        "pilot": {"arithmetic": "FLOAT proposals only", "status": status,
                  "evaluations": len(history), "maximum_evaluations": MAX_EVALUATIONS,
                  "log_marker_box": [-20, 4], "history": history},
        "directed_initial": old_receipt, "directed_selected": new_receipt,
        "directed_log2_gain_interval": str((old_bound/new_bound).log()/arb(2).log()),
        "tests": {"toy_moment_checks": checks, "full_packet_route_mass": "checked",
                  "directed_marker_weights": "exact rational reconstruction"},
        "precision_bits": args.precision, "seconds": time.monotonic()-started,
        "peak_process_bytes": peak_memory(), "generated_bytes": 0,
        "gradient_pin": GRADIENT_PIN,
        "source_sha256": hashlib.sha256(Path(__file__).read_bytes()).hexdigest()}, indent=2))


if __name__ == "__main__":
    main()
