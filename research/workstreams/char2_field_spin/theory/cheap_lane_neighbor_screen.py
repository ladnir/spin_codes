"""FLOAT shared-table neighbor screen and run-envelope diagnostic.

One fixed residue128 marker bank serves up to three selected weights in
508..516. Markers remain on the reference-r512 curve; for each actual(r,f)
the target is(f,r-2f), with its own exact route denominator and initial
source factor. Every selected weight retains exact ordered profile counts,
both terminal orientations, and all run/source count checks.

Alongside the exact profile-weighted sum, compute the conservative sum over
runs c of N(r,c) times the largest CONDITIONAL chain in that run. This
measures counting-envelope loss without confusing rare profiles with the
whole class. Unselected intermediate weights are not covered.

The explicit default15-center bank omits reference center256: its nearzero
singleton marker is unsuitable for neighboring positive-singleton targets.
No unresolved coefficient is clipped, skipped, or silently retried.
FLOAT only; negative output requires directed replay.180s/256MiB/20GiB,
stdout only, bytecode disabled, no generated data or automatic scans.
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
CELLS_PIN = "0ed129b48ee91cfc456c2e434ea5e279bd6ad777358f836d6152d2d2c9c813b1"
DEFAULT_CENTERS = "16,48,80,112,144,176,208,224,232,236,240,244,248,252,254"
MODULUS = 128


def load_cells():
    path = Path(__file__).with_name("cheap_lane_cauchy_cells.py")
    if hashlib.sha256(path.read_bytes()).hexdigest() != CELLS_PIN:
        raise RuntimeError("frozen cell screen hash mismatch")
    spec = importlib.util.spec_from_file_location("neighbor_cells", path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def build_marginals(coefficients, weights, centers, cf, cs, cells, fft, screen, np, guard, progress):
    max_full = max(weights)//2
    assignment = [min(centers, key=lambda c: (abs(f-c), c)) for f in range(max_full+1)]
    reference, _ = cells.assign_centers(centers)
    assert assignment[:min(257, len(assignment))] == reference[:len(assignment)]
    indices = [(f, s) for f in range(33) for s in range(33-f)]
    fd = np.asarray([f for f, _ in indices], dtype=int)
    sd = np.asarray([s for _, s in indices], dtype=int)
    local = np.asarray([coefficients[f, s] for f, s in indices])
    marginals = {r: [None]*(r//2+1) for r in weights}
    reports = []
    stats = {"max_imaginary_over_DC": 0.0, "max_negative_over_DC": 0.0,
             "min_selected_over_DC": 1.0}
    for center in centers:
        guard()
        targets = [(r, f) for r in weights for f in range(r//2+1) if assignment[f] == center]
        if not targets:
            reports.append({"center": center, "assigned_targets": 0})
            continue
        u, w = screen.marker_curve(center, cf, cs)
        lu = math.log(u.numerator)-math.log(u.denominator)
        lw = math.log(w.numerator)-math.log(w.denominator)
        mass = local*np.exp(fd*lu+sd*lw)[:, None]
        degree = np.zeros((MODULUS, MODULUS, 6))
        np.add.at(degree, (fd % MODULUS, sd % MODULUS), mass)
        p = degree[:, :, :4].reshape(MODULUS, MODULUS, 2, 2)
        t = degree[:, :, 4:].reshape(MODULUS, MODULUS, 2, 1)
        entries = np.sum(mass, axis=0)
        assert np.allclose(np.sum(p, axis=(0, 1)).ravel(), entries[:4], rtol=2e-13, atol=0)
        assert np.allclose(np.sum(t, axis=(0, 1)).ravel(), entries[4:], rtol=2e-13, atol=0)
        ordinary_p = screen.normalized_matrix(tuple(float(x) for x in entries[:4]))
        ordinary_t = screen.normalized_matrix(tuple(float(x) for x in entries[4:]))
        ordinary_b = screen.matrix_power(ordinary_p, 128)
        ordinary_v = screen.matrix_vector(screen.matrix_power(ordinary_p, 127), ordinary_t)
        spectra = fft.regional_spectra(p, t, np, guard)
        tables = []
        for kind, spectrum, expected in zip(("regular", "terminal"), spectra, (ordinary_b, ordinary_v)):
            raw, scale = spectrum
            try:
                fft.scaled_check(tuple(float(x.real) for x in raw[0, 0].ravel()), scale, expected)
                residue, _, _, diagnostic = fft.recover(spectrum, np)
                thresholds = cells.selection_thresholds(fft, residue, raw, np)
            except ArithmeticError as error:
                raise ArithmeticError(f"center={center},modulus128,{kind}: {error}") from error
            tables.append((residue, scale, kind, thresholds))
            for key in ("max_imaginary_over_DC", "max_negative_over_DC"):
                stats[key] = max(stats[key], diagnostic[key])
        for r, f in targets:
            guard()
            s = r-2*f
            target = (f % MODULUS, s % MODULUS)
            projected = []
            for residue, scale, kind, thresholds in tables:
                try:
                    values, minimum = cells.select_entries(residue, target, thresholds)
                except ArithmeticError as error:
                    raise ArithmeticError(f"weight={r},center={center},full={f},single={s},"
                                          f"modulus128,{kind}: {error}") from error
                projected.append(screen.normalized_matrix(values, scale))
                stats["min_selected_over_DC"] = min(stats["min_selected_over_DC"], minimum)
            denominator = math.log(comb(4096, f)*comb(4096-f, s))+f*lu+s*lw
            assert marginals[r][f] is None
            marginals[r][f] = (projected[0], projected[1], denominator, u, w)
        report = {"center": center, "u": str(u), "v": str(w), "assigned_targets": len(targets),
            "first_full": min(f for _, f in targets), "last_full": max(f for _, f in targets)}
        reports.append(report)
        progress("shared_marker_cell_complete", **report)
        del tables, spectra, degree, p, t
    assert all(all(x is not None for x in row) for row in marginals.values())
    return marginals, assignment, reports, stats


def join_weight(r, regions, counts, screen, log_source, guard, progress):
    """Exact indexed profile sum and its same-chain run-count envelope."""
    run_counts = [0]*(r+1)
    run_logs = [-math.inf]*(r+1)
    max_chains = [-math.inf]*(r+1)
    max_profiles = [None]*(r+1)
    source_total = ordered = positive = 0
    total_log = -math.inf
    top = []
    selected = None

    def add(a, b, number, chain):
        nonlocal source_total, ordered, positive, total_log, selected
        c = r-a-b
        weighted = log_source+math.log(number)+chain
        source_total += number; ordered += 1; positive += int(weighted >= 0)
        run_counts[c] += number
        run_logs[c] = screen.logadd(run_logs[c], weighted)
        total_log = screen.logadd(total_log, weighted)
        if chain > max_chains[c]:
            max_chains[c], max_profiles[c] = chain, (a, b)
        entry = (weighted, a, b, math.log(number), chain)
        if len(top) < 6:
            heapq.heappush(top, entry)
        elif weighted > top[0][0]:
            heapq.heapreplace(top, entry)
        if a == b == 223:
            selected = {"full": [a, b], "single": [r-2*a, r-2*b],
                "u": str(regions[a][3]), "v": str(regions[a][4]),
                "log2_weighted_profile_FLOAT": weighted/math.log(2),
                "log2_conditional_chain_FLOAT": chain/math.log(2)}

    for f0 in range(r//2+1):
        guard()
        for f1 in range(f0, r//2+1):
            number = counts.profile_count(4096, r, f0, f1)
            if not number:
                continue
            forward, reverse = screen.paired_chains(regions[f0], regions[f1])
            add(f0, f1, number, forward)
            if f0 != f1:
                add(f1, f0, number, reverse)
        if f0 % 64 == 0:
            progress("weight_profile_rows", weight=r, through_full0=f0, ordered_profiles=ordered)
    assert source_total == comb(8192, r)
    assert ordered == (r//2+1)**2-(1 if r % 2 == 0 else 0)
    envelope_log = -math.inf
    envelopes = [-math.inf]*(r+1)
    for c in range(r+1):
        number = counts.run_total(4096, r, c)
        assert run_counts[c] == number
        if number:
            assert max_chains[c] != -math.inf
            envelopes[c] = log_source+math.log(number)+max_chains[c]
            assert envelopes[c]+1e-8 >= run_logs[c]
            envelope_log = screen.logadd(envelope_log, envelopes[c])
    assert envelope_log+1e-8 >= total_log
    replay = -math.inf
    for value in run_logs:
        replay = screen.logadd(replay, value)
    assert abs(replay-total_log) < 1e-8
    logtwo = math.log(2)
    largest_runs = sorted(range(1, r+1), key=lambda c: envelopes[c], reverse=True)[:6]
    return {"weight": r, "ordered_profiles": ordered, "positive_profile_terms": positive,
        "log2_exact_profile_sum_FLOAT": total_log/logtwo,
        "log2_run_envelope_FLOAT": envelope_log/logtwo,
        "run_envelope_loss_bits_FLOAT": (envelope_log-total_log)/logtwo,
        "powered_source_log2_FLOAT": log_source/logtwo,
        "top_weighted_profiles": [{"full": [a, b], "single": [r-2*a, r-2*b],
            "runs": r-a-b, "log2_bound_FLOAT": value/logtwo,
            "log2_exact_source_count": count/logtwo, "log2_conditional_chain_FLOAT": chain/logtwo}
            for value, a, b, count, chain in sorted(top, reverse=True)],
        "top_run_envelopes": [{"runs": c, "maximizing_ordered_profile": max_profiles[c],
            "log2_envelope_FLOAT": envelopes[c]/logtwo,
            "log2_exact_profile_sum_FLOAT": run_logs[c]/logtwo,
            "loss_bits_FLOAT": (envelopes[c]-run_logs[c])/logtwo} for c in largest_runs],
        "selected_f223": selected, "exact_count_checks": {
            "all_run_totals": True, "source_sum_equals_binomial8192_weight": True,
            "source_count_log2": math.log(source_total)/logtwo,
            "both_terminal_orientations": True, "updates": 32767}}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--weights", default="512")
    parser.add_argument("--centers", default=DEFAULT_CENTERS)
    parser.add_argument("--z", default="17447/25000")
    parser.add_argument("--alpha", default="119/1280")
    parser.add_argument("--curve-full", default="1562175/990029")
    parser.add_argument("--curve-single", default="1140210/880327")
    args = parser.parse_args()
    try:
        weights = tuple(sorted(int(x) for x in args.weights.split(",")))
        centers = tuple(sorted(int(x) for x in args.centers.split(",")))
        z, alpha, cf, cs = map(F, (args.z, args.alpha, args.curve_full, args.curve_single))
    except (ValueError, ZeroDivisionError):
        parser.error("invalid integer list or rational parameter")
    if not 1 <= len(weights) <= 3 or len(set(weights)) != len(weights) or not 508 <= weights[0] <= weights[-1] <= 516:
        parser.error("choose1..3 distinct weights in508..516")
    if not 1 <= len(centers) <= 16 or len(set(centers)) != len(centers) or not 0 <= centers[0] <= centers[-1] <= 256:
        parser.error("choose1..16 distinct reference-r512 centers in0..256")
    if not 0 < z < 1 or not 0 < alpha <= 1 or not all(F(1, 256) <= x <= 256 for x in (cf, cs)):
        parser.error("require0<z<1,0<alpha<=1,curve multipliers in[1/256,256]")
    if max(n.bit_length() for x in (z, alpha, cf, cs)
           for n in (x.numerator, x.denominator)) > 64:
        parser.error("rational numerator/denominator limit64bits")
    if shutil.disk_usage(Path(__file__).resolve().anchor).free < 20*(1 << 30):
        raise RuntimeError("requires20GiB free on source volume")
    started = time.monotonic(); deadline = started+180
    cells = load_cells(); fft = cells.load_fft(); screen = fft.load_screen()
    calibration = screen.load("cheap_lane_cauchy_marker_calibration.py", screen.CALIBRATION_PIN, "neighbor_calibration")
    gradient = calibration.load_gradient(); anchor = gradient.load_anchor()
    peak_memory = anchor.peak_memory_reader()

    def guard():
        if time.monotonic() > deadline:
            raise TimeoutError("180-second neighbor screen exhausted; no automatic retry")
        if peak_memory() > 256*(1 << 20):
            raise MemoryError("256MiB neighbor screen budget exceeded")

    def progress(phase, **values):
        print(json.dumps({"phase": phase, "seconds": time.monotonic()-started, **values}), flush=True)

    counts = screen.load("cheap_lane_fast_profile_counts.py", screen.COUNT_PIN, "neighbor_counts")
    counts.DEADLINE = deadline
    tests = {**counts.checks(), "noncommuting_terminal_checks": screen.chain_tests(),
             "exact_scale_tests": cells.exact_scale_tests()}
    import numpy as np
    from flint import ctx
    ctx.prec = 256
    tests["exact_FFT_toy_coefficients"] = fft.exact_toy(np, guard)
    adapter = anchor.load_adapter(); core = adapter.load_verifier()
    core.Z, core.ALPHA = z, alpha

    class NeighborVerifier(core.Verifier):
        def guard(self):
            guard()

    verifier = NeighborVerifier(2, deadline)
    caps = verifier.local_caps(list(combinations_with_replacement(range(17), 4)))
    exponent = verifier.exact(alpha)
    powered = {occ: tuple(verifier.positive(x)**exponent for x in row) for occ, row in caps.items()}
    del caps
    coefficients = calibration.float_coefficients(powered, np, guard)
    del powered
    marginals, assignment, bank, fft_stats = build_marginals(
        coefficients, weights, centers, cf, cs, cells, fft, screen, np, guard, progress)
    results = []; joint_exact = joint_envelope = -math.inf
    for r in weights:
        guard()
        h, c_h, eta = anchor.source_factor(r)
        source_ball = (verifier.exact(eta)*verifier.exact(core.BETA)**r
                       *verifier.exact(z)**(-core.LIMIT))**exponent
        log_source = float(source_ball.log().mid())
        result = join_weight(r, marginals[r], counts, screen, log_source, guard, progress)
        result["joint_missing"] = {"h": h, "c_h_exact": str(c_h), "eta_exact": str(eta)}
        results.append(result)
        joint_exact = screen.logadd(joint_exact, result["log2_exact_profile_sum_FLOAT"]*math.log(2))
        joint_envelope = screen.logadd(joint_envelope, result["log2_run_envelope_FLOAT"]*math.log(2))
        progress("weight_complete", **result)
    guard()
    print(json.dumps({"scope": "complete source supports only at the explicitly selected weights",
        "directed": False, "negative_is_not_a_certificate": True, "whole_code_certificate": False,
        "unselected_weights_covered": False, "parameters": {"weights": weights, "centers": centers,
            "marker_reference_weight": 512, "moduli": [MODULUS, MODULUS],
            "z": str(z), "alpha": str(alpha), "curve_full": str(cf), "curve_single": str(cs),
            "distance_threshold": core.LIMIT, "fallback_or_scanning": False},
        "selected_weight_joint_log2_exact_sum_FLOAT": joint_exact/math.log(2),
        "selected_weight_joint_log2_run_envelope_FLOAT": joint_envelope/math.log(2),
        "results": results, "bank": bank, "assignment_by_full": assignment,
        "FFT_controls": fft_stats, "tests": tests,
        "seconds": time.monotonic()-started, "peak_process_bytes": peak_memory(), "generated_bytes": 0,
        "cells_pin": CELLS_PIN, "FFT_pin": cells.FFT_PIN, "screen_pin": fft.SCREEN_PIN,
        "calibration_pin": screen.CALIBRATION_PIN, "count_pin": screen.COUNT_PIN,
        "core_pin": anchor.CORE_PIN,
        "source_sha256": hashlib.sha256(Path(__file__).read_bytes()).hexdigest()}, indent=2))


if __name__ == "__main__":
    main()
