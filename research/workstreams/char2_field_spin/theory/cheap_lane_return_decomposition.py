"""Selected-profile FLOAT decomposition of the frozen route-aware majorant.

Keep the local caps, alpha, own-f marker curve, source count, route
normalizers, residue modulus and ordered terminal fixed. Multiply only the
POWERED return channel E (live->quiet) by nu=1,2^-16,2^-32 before taking
regional powers. Thus a path with k returns receives nu^k, not nu^(alpha*k).
No global correction is applied. The nu<1 values are positive submajorant
diagnostics, NOT upper bounds for the code or for a support class.

Also extract the pure-quiet path through scalar Q^128 and Q^127*T0
polynomials. This is a positive component of the SAME residue-projected
majorant. Independent scalar scaling prevents the live channel from
obscuring its FFT coefficients. No subtraction is used to construct a
nonquiet bound; quiet/total log ratios are descriptive only.

Defaults r1024, ordered profiles(384,500),(384,384), residue512,
z304397809/625000000, alpha119/1280. Fixed own-f markers; no calibration or
parameter scan. Strict frozen FFT noise checks remain FLOAT diagnostics,
not interval certificates.120s/256MiB/20GiB, stdout only, no generated data.
"""

import argparse
from collections import defaultdict
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
CALIBRATION_PIN = "e9f90bdda919247df80e71e5d2555f24ef7ca12323d95b1df622b22556327c56"
SECONDS, MEMORY_LIMIT, LN2 = 120, 256*(1 << 20), math.log(2)
RETURN_MARKS = (F(1), F(1, 1 << 16), F(1, 1 << 32))


def load_calibration():
    path = Path(__file__).with_name("cheap_lane_broad_profile_calibration.py")
    if hashlib.sha256(path.read_bytes()).hexdigest() != CALIBRATION_PIN:
        raise RuntimeError("frozen broad-profile calibration hash mismatch")
    spec = importlib.util.spec_from_file_location("return_decomposition_calibration", path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def scalar_spectra(q, terminal, np, fft, guard, doublings=7):
    """Positive scalar convolutions; one scale for the entire FFT array."""
    def normalize(values, scale):
        dc = values[0, 0]
        mass = float(dc.real)
        if not mass > 0 or not np.isfinite(values).all():
            raise ArithmeticError("invalid quiet scalar DC mass")
        if abs(float(dc.imag)) > fft.RELATIVE_TOLERANCE*mass:
            raise ArithmeticError("imaginary quiet scalar DC mass")
        return values/mass, scale+math.log(mass)

    qh, qscale = normalize(np.fft.fft2(q), 0.0)
    th, tscale = normalize(np.fft.fft2(terminal), 0.0)
    for _ in range(doublings):
        guard()
        th, tscale = normalize(qh*th, qscale+tscale)
        qh, qscale = normalize(qh*qh, 2*qscale)
    return (qh, qscale), (th, tscale)


def scalar_selected(spectrum, target, expected_log, np, fft):
    raw, scale = spectrum
    dc = float(raw[0, 0].real)
    if abs(math.log(dc)+scale-expected_log) > 5e-9:
        raise ArithmeticError("quiet scalar DC/ordinary disagreement")
    residue = np.fft.ifft2(raw)
    imaginary = float(np.max(np.abs(residue.imag)))/dc
    negative = max(0.0, -float(np.min(residue.real)))/dc
    if max(imaginary, negative) > fft.RELATIVE_TOLERANCE:
        raise ArithmeticError("material imaginary/negative quiet FFT residue")
    if abs(float(np.sum(residue.real))-dc) > fft.RELATIVE_TOLERANCE*dc:
        raise ArithmeticError("quiet scalar residue sum/DC disagreement")
    value = float(residue[target].real)
    relative = value/dc
    noise = max(imaginary, negative, 32*np.finfo(float).eps)
    if relative <= max(fft.SELECTED_FLOOR, fft.NOISE_MULTIPLIER*noise):
        raise ArithmeticError("quiet selected coefficient numerically unresolved")
    return scale+math.log(value), {"max_imaginary_over_DC": imaginary,
        "max_negative_over_DC": negative, "selected_over_DC": relative}


class QuietReplay:
    def __init__(self, coefficients, r, modulus, calibration, fft, np, guard):
        self.r, self.modulus = r, modulus
        self.calibration, self.fft, self.np, self.guard = calibration, fft, np, guard
        indices = [(f, s) for f in range(33) for s in range(33-f)]
        self.fd = np.asarray([f for f, s in indices])
        self.sd = np.asarray([s for f, s in indices])
        local = np.asarray([[coefficients[f, s, j] for j in (0, 4)] for f, s in indices])
        self.positive = local > 0
        self.loglocal = np.full(local.shape, -math.inf)
        self.loglocal[self.positive] = np.log(local[self.positive])
        self.cache = {}

    def marginal(self, full, markers):
        key = (full, markers)
        if key in self.cache: return self.cache[key]
        np, fft, m = self.np, self.fft, self.modulus
        self.guard()
        lu, lv = map(self.calibration.logq, markers)
        single = self.r-2*full
        logs = self.loglocal+(self.fd*lu+self.sd*lv)[:, None]
        shifts = np.max(logs, axis=0)
        mass = np.exp(logs-shifts[None, :])
        if np.any(self.positive & (mass == 0)):
            raise ArithmeticError("positive quiet local coefficient underflow")
        degree = np.zeros((m, m, 2))
        np.add.at(degree, (self.fd % m, self.sd % m), mass)
        packet_logs = shifts+np.log(np.sum(mass, axis=0))
        denominator = math.log(comb(4096, full)*comb(4096-full, single))+full*lu+single*lv
        try:
            spectra = scalar_spectra(degree[:, :, 0], degree[:, :, 1], np, fft, self.guard)
            selected, diagnostics = [], []
            scale_additions = (128*shifts[0], 127*shifts[0]+shifts[1])
            expected_logs = (128*packet_logs[0], 127*packet_logs[0]+packet_logs[1])
            for (raw, scale), extra, expected in zip(spectra, scale_additions, expected_logs):
                value, diagnostic = scalar_selected((raw, scale+float(extra)),
                    (full % m, single % m), float(expected), np, fft)
                selected.append(value); diagnostics.append(diagnostic)
            answer = {"status": "resolved", "regional_logs": tuple(selected),
                "denominator": denominator, "FFT_controls": diagnostics}
        except ArithmeticError as error:
            answer = {"status": "unresolved", "reason": str(error)}
        self.cache[key] = answer
        return answer

    def profile(self, fulls, markers, log_source_count):
        values = [self.marginal(f, m) for f, m in zip(fulls, markers)]
        if any(x["status"] != "resolved" for x in values):
            return {"status": "unresolved", "reasons": [x.get("reason") for x in values]}
        (q0, t0), (q1, t1) = [x["regional_logs"] for x in values]
        conditional = 128*q0+127*q1+t1-128*sum(x["denominator"] for x in values)
        return {"status": "resolved", "log2_component_FLOAT": (log_source_count+conditional)/LN2,
            "log2_conditional_chain_FLOAT": conditional/LN2,
            "FFT_controls": [x["FFT_controls"] for x in values]}


def exact_tests(np, fft, guard):
    """Exact return exponents/terminal and scalar polynomial FFT controls."""
    base = ((F(2, 7), F(3, 11)), (F(5, 13), F(7, 17)))
    terminal = (F(11, 19), F(13, 23))
    checks, totals = 0, []
    for nu in (F(0),)+tuple(reversed(RETURN_MARKS)):
        matrix = ((base[0][0], base[0][1]), (nu*base[1][0], base[1][1]))
        row = (F(1), F(0))
        for _ in range(4):
            row = tuple(sum(row[i]*matrix[i][j] for i in (0, 1)) for j in (0, 1))
        actual = sum(row[j]*terminal[j] for j in (0, 1))
        expected = F(0); quiet = F(0)
        for states in product((0, 1), repeat=4):
            path = (0,)+states
            returns = sum(path[i] == 1 and path[i+1] == 0 for i in range(4))
            mass = terminal[states[-1]]
            for i in range(4): mass *= base[path[i]][path[i+1]]
            expected += mass*nu**returns
            if not any(states): quiet += mass
        assert actual == expected and quiet == base[0][0]**4*terminal[0]
        assert actual >= quiet
        totals.append(actual); checks += 3
    assert totals == sorted(totals); checks += 1

    q = {(0, 0): F(1, 3), (1, 0): F(2, 5)}
    terminal_poly = {(0, 1): F(3, 7), (1, 1): F(5, 11)}
    def multiply(a, b):
        out = defaultdict(F)
        for (f, s), x in a.items():
            for (ff, ss), y in b.items(): out[f+ff, s+ss] += x*y
        return dict(out)
    q2 = multiply(q, q)
    truth = (multiply(q2, q2), multiply(multiply(q2, q), terminal_poly))
    arrays = []
    for poly in (q, terminal_poly):
        array = np.zeros((4, 4))
        for (f, s), value in poly.items(): array[f % 4, s % 4] += float(value)
        arrays.append(array)
    for (raw, scale), expected in zip(scalar_spectra(*arrays, np, fft, guard, doublings=2), truth):
        recovered = np.fft.ifft2(raw)*math.exp(scale)
        exact = np.zeros((4, 4))
        for (f, s), value in expected.items(): exact[f % 4, s % 4] += float(value)
        assert np.max(np.abs(recovered-exact)) < 2e-12
        checks += exact.size
    return checks


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--weight", type=int, default=1024)
    parser.add_argument("--profiles", default="384:500,384:384")
    parser.add_argument("--modulus", type=int, choices=(128, 256, 512), default=512)
    parser.add_argument("--z", type=F, default=F(304397809, 625000000))
    parser.add_argument("--alpha", type=F, default=F(119, 1280))
    args = parser.parse_args(); r = args.weight
    try: profiles = [tuple(map(int, item.split(":"))) for item in args.profiles.split(",")]
    except ValueError: parser.error("profiles must be comma-separated full0:full1 pairs")
    if not 7 <= r <= 4096 or not 1 <= len(profiles) <= 2 or len(set(profiles)) != len(profiles):
        parser.error("require r7..4096 and one or two distinct profiles")
    if any(len(pair) != 2 or not all(0 <= f <= r//2 for f in pair) for pair in profiles):
        parser.error("full counts must lie in0..floor(r/2)")
    if not 0 < args.z < 1 or not 0 < args.alpha <= 1:
        parser.error("require0<z<1 and0<alpha<=1")
    if max(x.bit_length() for v in (args.z, args.alpha) for x in (v.numerator, v.denominator)) > 64:
        parser.error("rational witnesses limited to64 bits")
    if shutil.disk_usage(Path(__file__).resolve().anchor).free < 20*(1 << 30):
        raise RuntimeError("requires20GiB free on source volume")
    started = time.monotonic(); deadline = started+SECONDS
    calibration = load_calibration(); calibration.MODULI = (args.modulus,)
    broad = calibration.load_broad()
    neighbor = broad.load("cheap_lane_neighbor_screen.py", broad.NEIGHBOR_PIN, "decomposition_neighbor")
    cells = neighbor.load_cells(); fft = cells.load_fft(); screen = fft.load_screen()
    old = screen.load("cheap_lane_cauchy_marker_calibration.py", screen.CALIBRATION_PIN, "decomposition_old")
    gradient = old.load_gradient(); anchor = gradient.load_anchor(); peak = anchor.peak_memory_reader()

    def guard():
        if time.monotonic() > deadline: raise TimeoutError("120-second decomposition deadline; no retry")
        if peak() > MEMORY_LIMIT: raise MemoryError("256MiB decomposition budget exceeded")

    def progress(phase, **values):
        print(json.dumps({"phase": phase, "seconds": time.monotonic()-started, **values}), flush=True)

    import numpy as np
    from flint import ctx
    ctx.prec = 256
    tests = exact_tests(np, fft, guard)+fft.exact_toy(np, guard)+cells.exact_scale_tests()
    counts = screen.load("cheap_lane_fast_profile_counts.py", screen.COUNT_PIN, "decomposition_counts")
    counts.DEADLINE = deadline
    adapter = anchor.load_adapter(); core = adapter.load_verifier()
    core.Z, core.ALPHA = args.z, args.alpha
    class LocalVerifier(core.Verifier):
        def guard(self): guard()
    verifier = LocalVerifier(2, deadline)
    caps = verifier.local_caps(list(combinations_with_replacement(range(17), 4)))
    exponent = verifier.exact(args.alpha)
    powered = {occ: tuple(verifier.positive(value)**exponent for value in row) for occ, row in caps.items()}
    del caps
    for row in powered.values():
        for value in row:
            x = float(value.mid())
            if not math.isfinite(x) or x < 0 or (value > 0 and x == 0):
                raise ArithmeticError("positive powered cap not representable in FLOAT")
    coefficients = old.float_coefficients(powered, np, guard); del powered
    h, c_h, eta = anchor.source_factor(r)
    source = exponent*(verifier.exact(eta).log()+r*verifier.exact(core.BETA).log()
                       -core.LIMIT*verifier.exact(args.z).log())
    logsource = float(source.mid())
    cf, cs = F(1562175, 990029), F(1140210, 880327)
    quiet = QuietReplay(coefficients, r, args.modulus, calibration, fft, np, guard)
    replays = []
    for nu in RETURN_MARKS:
        modified = coefficients.copy()
        modified[:, :, 2] *= float(nu)
        replays.append(calibration.Replay(modified, r, cells, fft, screen, np, guard))
    progress("cached_local_caps_ready", tests=tests)
    for fulls in profiles:
        guard(); number = counts.profile_count(4096, r, *fulls)
        if not number:
            progress("empty_actual_profile", fulls=fulls); continue
        markers = tuple(broad.marker_curve(r, f, cf, cs) for f in fulls)
        log_source_count = logsource+math.log(number)
        pure_quiet = quiet.profile(fulls, markers, log_source_count)
        results = []
        baseline = None
        for nu, replay in zip(RETURN_MARKS, replays):
            replayed = replay.profile(fulls, markers, log_source_count)
            result = {"return_mark_after_alpha": str(nu), "results": replayed,
                "code_bound_claimed": False}
            selected = replayed["residue"+str(args.modulus)]
            if selected["status"] == "resolved":
                value = selected["log2_bound_FLOAT"]
                if baseline is None and nu == 1: baseline = value
                if baseline is not None:
                    result["log2_reduction_from_unmarked_FLOAT"] = baseline-value
                    assert value <= baseline+2e-5
                if pure_quiet["status"] == "resolved":
                    qvalue = pure_quiet["log2_component_FLOAT"]
                    assert value >= qvalue-2e-5
                    result["log2_quiet_fraction_FLOAT"] = qvalue-value
            results.append(result)
            progress("selected_return_mark_complete", fulls=fulls, **result)
        progress("selected_decomposition_complete", fulls=fulls,
            markers=[[str(x) for x in pair] for pair in markers],
            exact_source_count_log2_FLOAT=math.log(number)/LN2,
            exact_source_count_sha256=hashlib.sha256(str(number).encode()).hexdigest(),
            pure_quiet=pure_quiet,
            interpretation="nu<1 suppresses paths by their number of returns; it is not a code bound",
            no_return_limit_computed=False)
    guard()
    print(json.dumps({"scope": "selected ordered profiles; FLOAT majorant decomposition only",
        "whole_code_certificate": False, "complete_weight_bound": False, "directed": False,
        "weight": r, "profiles": profiles, "z": str(args.z), "alpha": str(args.alpha),
        "modulus": args.modulus, "return_marks_after_alpha": list(map(str, RETURN_MARKS)),
        "source_factor": {"h": h, "c_h": str(c_h), "eta": str(eta),
            "powered_source_log2_FLOAT": logsource/LN2},
        "finite_updates": 32767, "terminal_packets": 1, "terminal_marked": False,
        "count_policy": "same exact indexed count for all components of each ordered profile",
        "markers": "fixed own-f curve; no proposal or calibration",
        "tests": tests, "calibration_pin": CALIBRATION_PIN, "core_pin": anchor.CORE_PIN,
        "count_pin": screen.COUNT_PIN, "FFT_pin": cells.FFT_PIN,
        "seconds": time.monotonic()-started, "peak_process_bytes": peak(), "generated_bytes": 0,
        "source_sha256": hashlib.sha256(Path(__file__).read_bytes()).hexdigest()}, indent=2), flush=True)


if __name__ == "__main__":
    main()
