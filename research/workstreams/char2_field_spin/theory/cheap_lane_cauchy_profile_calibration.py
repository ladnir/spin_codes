"""One equal-matching r512 profile: bounded marker proposals and Arb replay.

Default full count223 selects actual profiles(223,66) in both matchings.
Keep z7/10, alpha119/1280 and the EXACT indexed profile count fixed while
changing only the two positive route markers. The starting curve is
 (1562175/990029)*f/(3584+f), (1140210/880327)*s/(3584+f).
Zero markers are replaced by2^-30, not by an exact boundary extraction.

Reuse the frozen calibrator's cached local powered caps, positive FLOAT
finite-chain moments and at most12 guided proposals. The chosen markers
are rounded to rationals and independently replayed with exact rational
route weights and directed Arb arithmetic. The final chain P^32767 T
retains the output-only terminal, with denominator[M(f,s)u^f v^s]^256.
Source Omega, state paths, alpha placement and source count do not change.
The optional run-count control changes only the displayed count multiplier.

An empty matching profile is reported as such without a numerical search.
At a zero target degree, keep that marker fixed and collapse its FLOAT
proposal axis. The final directed evaluation still includes every local
route term. No FLOAT value supplies a certificate. Selected profile only;
no parameter grid, no whole-code claim,60s/128MiB/20GiB/stdout-only.
"""

import argparse
from fractions import Fraction as F
import hashlib
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
CURVE_FULL, CURVE_SINGLE = F(1562175, 990029), F(1140210, 880327)
ZERO_MARKER = F(1, 1 << 30)
WEIGHT = 512


def load(name, pin, tag):
    path = Path(__file__).with_name(name)
    if hashlib.sha256(path.read_bytes()).hexdigest() != pin:
        raise RuntimeError("frozen helper hash mismatch: "+name)
    spec = importlib.util.spec_from_file_location(tag, path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def initial_markers(full):
    single = WEIGHT-2*full
    u = CURVE_FULL*F(full, 3584+full)
    v = CURVE_SINGLE*F(single, 3584+full)
    return (u if u else ZERO_MARKER), (v if v else ZERO_MARKER)


def freeze_zero_axes(coefficients, targets, markers, np):
    """Proposal-only evaluation of a marker whose target degree is zero."""
    out = coefficients.copy()
    if targets[0] == 0:
        weights = np.power(float(markers[0]), np.arange(33))
        evaluated = np.sum(out*weights[:, None, None], axis=0)
        out[:] = 0; out[0] = evaluated
    if targets[1] == 0:
        weights = np.power(float(markers[1]), np.arange(33))
        evaluated = np.sum(out*weights[None, :, None], axis=1)
        out[:] = 0; out[:, 0] = evaluated
    return out


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--full", type=int, default=223)
    parser.add_argument("--z", default="7/10")
    parser.add_argument("--alpha", default="119/1280")
    parser.add_argument("--precision", type=int, choices=(256, 384, 512), default=256)
    parser.add_argument("--run-count-control", action="store_true")
    args = parser.parse_args()
    z, alpha = F(args.z), F(args.alpha)
    if not 0 <= args.full <= 256:
        parser.error("full count must be0..256")
    if not 0 < z < 1 or not 0 < alpha <= 1:
        parser.error("require0<z<1 and0<alpha<=1")
    if any(max(x.numerator.bit_length(), x.denominator.bit_length()) > 64 for x in (z, alpha)):
        parser.error("rational witnesses must fit64 bits")
    if shutil.disk_usage(Path(__file__).resolve().anchor).free < 20*(1 << 30):
        raise RuntimeError("requires20GiB free on source volume")
    started = time.monotonic(); deadline = started+60
    calibration = load("cheap_lane_cauchy_marker_calibration.py", CALIBRATION_PIN, "selected_profile_calibration")
    gradient = calibration.load_gradient(); anchor = gradient.load_anchor()
    peak_memory = anchor.peak_memory_reader()

    def guard():
        if time.monotonic() > deadline:
            raise TimeoutError("60-second profile calibration exhausted; no automatic retry")
        if peak_memory() > anchor.MEMORY_LIMIT:
            raise MemoryError("128MiB process-memory budget exceeded")

    counts = load("cheap_lane_fast_profile_counts.py", COUNT_PIN, "selected_profile_counts")
    counts.DEADLINE = deadline
    full, single = args.full, WEIGHT-2*args.full
    number = counts.profile_count(4096, WEIGHT, full, full)
    if number == 0:
        guard()
        print(json.dumps({"scope": "empty actual equal-matching profile",
            "source_weight": WEIGHT, "ordered_profile": [[full, single], [full, single]],
            "exact_source_count": 0, "numerical_search_run": False,
            "whole_code_certificate": False, "generated_bytes": 0,
            "count_pin": COUNT_PIN, "seconds": time.monotonic()-started,
            "source_sha256": hashlib.sha256(Path(__file__).read_bytes()).hexdigest()}, indent=2))
        return
    import numpy as np
    from flint import arb, ctx
    ctx.prec = args.precision
    checks = calibration.proposal_tests(np)+gradient.exact_tests()
    adapter = anchor.load_adapter(); core = adapter.load_verifier()

    class ProfileVerifier(core.Verifier):
        def guard(self):
            guard()

    verifier = ProfileVerifier(2, deadline)
    core.Z, core.ALPHA = z, alpha
    shapes = list(combinations_with_replacement(range(17), 4))
    caps = verifier.local_caps(shapes)
    exponent = verifier.exact(alpha)
    powered = {occ: tuple(verifier.positive(x)**exponent for x in values) for occ, values in caps.items()}
    del caps
    h, c, eta = anchor.source_factor(WEIGHT)
    source = (verifier.exact(eta)*verifier.exact(core.BETA)**WEIGHT
              *verifier.exact(z)**(-core.LIMIT))**exponent
    normalizer = comb(4096, full)*comb(4096-full, single)
    constant = float(((number*source).log()-256*arb(normalizer).log()).mid())
    coefficients = calibration.float_coefficients(powered, np, guard)
    old_markers = initial_markers(full)
    targets = np.asarray([256*full, 256*single], dtype=float)
    proposal_coefficients = freeze_zero_axes(coefficients, targets, old_markers, np)
    # A frozen axis is already evaluated into the FLOAT coefficients. Its
    # optimizer coordinate is immaterial and is discarded after selection.
    start = [float(old_markers[i]) if targets[i] else 1.0 for i in range(2)]
    best, history, status = calibration.calibrate(proposal_coefficients, constant, targets,
                                                 start, np, guard)
    chosen = tuple(max(ZERO_MARKER, F(math.exp(float(best["x"][i]))).limit_denominator(1 << 20))
                   if targets[i] else old_markers[i] for i in range(2))
    if any(x <= 0 for x in chosen):
        raise ArithmeticError("nonpositive rounded marker")
    run_number = counts.run_total(4096, WEIGHT, WEIGHT-2*full)
    assert number <= run_number

    def directed_receipt(markers):
        u, v = markers
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
        value = gradient.dual_matmul(gradient.dual_matpow(p, 32767, arb(0), arb(1)), t, arb(0))[0][0]
        denominator = verifier.exact(u)**full*verifier.exact(v)**single*normalizer
        bound = number*source*value[0]/denominator**256
        logtwo = arb(2).log()
        receipt = {"markers": [str(u), str(v)], "log2_bound_interval": str(bound.log()/logtwo),
            "regional_mean_intervals": [str(value[1]/value[0]/256), str(value[2]/value[0]/256)],
            "log_marker_gradient_intervals": [str((value[1]/value[0]-256*full)/logtwo),
                                               str((value[2]/value[0]-256*single)/logtwo)],
            "selected_profile_proved_below_2_pow_minus_40": bool(bound < arb(2)**(-40))}
        if args.run_count_control:
            control = bound*run_number/number
            receipt["run_count_control_log2_interval"] = str(control.log()/logtwo)
            receipt["run_count_control_scope"] = "same selected profile with a larger count, not all run-class profiles"
        return bound, receipt

    old_bound, old_receipt = directed_receipt(old_markers)
    new_bound, new_receipt = directed_receipt(chosen)
    guard()
    print(json.dumps({"scope": "one selected actual equal-matching r512 profile only",
        "directed_final_replay": True, "whole_code_certificate": False,
        "all_profiles_at_weight_covered": False,
        "parameters": {"source_weight": WEIGHT, "ordered_profile": [[full, single], [full, single]],
                       "runs": WEIGHT-2*full, "z": str(z), "alpha": str(alpha), "updates": 32767},
        "exact_profile_source_count": number,
        "profile_source_count_log2_interval": str(arb(number).log()/arb(2).log()),
        "source_count_policy": "same exact indexed profile count in initial and selected receipts",
        "pilot": {"arithmetic": "FLOAT proposals only", "status": status,
                  "evaluations": len(history), "maximum_evaluations": calibration.MAX_EVALUATIONS,
                  "zero_target_markers_frozen": [bool(x == 0) for x in targets], "history": history},
        "directed_initial": old_receipt, "directed_selected": new_receipt,
        "directed_log2_gain_interval": str((old_bound/new_bound).log()/arb(2).log()),
        "tests": {"toy_moment_checks": checks,
                  "packet_route_mass": "checked", "directed_marker_weights": "exact rational reconstruction"},
        "precision_bits": args.precision, "seconds": time.monotonic()-started,
        "peak_process_bytes": peak_memory(), "generated_bytes": 0,
        "calibration_pin": CALIBRATION_PIN, "count_pin": COUNT_PIN,
        "source_sha256": hashlib.sha256(Path(__file__).read_bytes()).hexdigest()}, indent=2))


if __name__ == "__main__":
    main()
