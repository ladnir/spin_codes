"""FLOAT fixed-marker-cell r512 screen using one residue128 table per cell.

Nearest-center assignment is deterministic, with ties sent to the smaller
center. Each assigned target(f,512-2f) is extracted from that center's common
P^128/P^127*T residue tables and divided by ITS OWN exact profile/marker
denominator. This is coefficient majorization, not profile interpolation.
No new alpha, marker union, or changing center inside a state path is used.

Only the pinned all-profile screen's build_regions seam is replaced. Its
exact indexed counts, both terminal orientations, and complete source-sum
checks remain unchanged. FLOAT controls/noise policy are imported from the
audited FFT screen. A failed/unresolved marginal aborts: no clipping,
skipping, marker fitting, or silent fallback. Directed replay is separate.
180s/256MiB/20GiB, stdout only, no generated files or numerical scans.
"""

import argparse
from fractions import Fraction as F
import hashlib
import importlib.util
import json
import math
from math import comb
from pathlib import Path
import sys
import time

sys.dont_write_bytecode = True
FFT_PIN = "0b28a10417eec7a1f6669137f2f12fc805e1b694388777fd81040d1c816a1ed9"
DEFAULT_CENTERS = "16,48,80,112,144,176,208,224,240,252"
MODULUS = 128


def load_fft():
    path = Path(__file__).with_name("cheap_lane_cauchy_residue_screen.py")
    if hashlib.sha256(path.read_bytes()).hexdigest() != FFT_PIN:
        raise RuntimeError("frozen FFT screen hash mismatch")
    spec = importlib.util.spec_from_file_location("cell_fft", path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def exact_scale_tests():
    """Exact normalization and positive-alias/min sanity checks."""
    coefficients = {(0, 0): F(1, 3), (1, 2): F(2, 5), (2, 1): F(3, 7),
                    (5, 2): F(5, 11), (1, 6): F(7, 13)}
    target, divisor = (1, 2), 17
    actual = coefficients[target]/divisor
    uppers = []
    for u, v in ((F(1, 3), F(2, 5)), (F(5, 4), F(1, 2))):
        scaled_exact = coefficients[target]*u**target[0]*v**target[1]
        assert scaled_exact/(divisor*u**target[0]*v**target[1]) == actual
        aliases = sum((mass*u**f*v**s for (f, s), mass in coefficients.items()
                       if (f % 4, s % 4) == target), F(0))
        upper = aliases/(divisor*u**target[0]*v**target[1])
        assert upper >= actual
        uppers.append(upper)
    assert min(uppers) >= actual
    return 5


def selection_thresholds(fft, residue, spectrum, np):
    """Compute each entry's reliability threshold once per cached table."""
    thresholds = []
    for i in range(2):
        for j in range(spectrum.shape[3]):
            dc = float(spectrum[0, 0, i, j].real)
            if not dc > 0:
                raise ArithmeticError("nonpositive entry DC mass")
            block = residue[:, :, i, j]
            imaginary = float(np.max(np.abs(block.imag)))/dc
            negative = max(0.0, -float(np.min(block.real)))/dc
            noise = max(imaginary, negative, 32*np.finfo(float).eps)
            thresholds.append((i, j, dc, max(fft.SELECTED_FLOOR, fft.NOISE_MULTIPLIER*noise)))
    return thresholds


def select_entries(residue, target, thresholds):
    """Same mandatory per-entry DC noise test as the frozen FFT screen."""
    selected = []
    minimum = 1.0
    for i, j, dc, threshold in thresholds:
        block = residue[:, :, i, j]
        value = float(block[target].real)
        relative = value/dc
        if relative <= threshold:
            raise ArithmeticError("selected coefficient numerically unresolved")
        selected.append(value)
        minimum = min(minimum, relative)
    return tuple(selected), minimum


def assign_centers(centers):
    """Pure assignment interface shared by a later directed wrapper."""
    assert centers and tuple(sorted(set(centers))) == tuple(centers)
    assert 0 <= centers[0] <= centers[-1] <= 256
    assignment = [min(centers, key=lambda c: (abs(f-c), c)) for f in range(257)]
    assigned = {center: [f for f in range(257) if assignment[f] == center] for center in centers}
    assert all(assigned[c] for c in centers)
    assert sum(map(len, assigned.values())) == 257
    return assignment, assigned


def install(screen, fft, centers, z, alpha, metadata):
    assignment, assigned = assign_centers(centers)
    metadata["assignment"] = [{"center": c, "first_full": min(assigned[c]),
        "last_full": max(assigned[c]), "marginals": len(assigned[c])} for c in centers]

    def build_regions(coefficients, curve_full, curve_single, np, guard):
        started = time.monotonic()
        toy = fft.exact_toy(np, guard)
        scale_tests = exact_scale_tests()
        assignment_test, _ = assign_centers((0, 16, 256))
        assert assignment_test[8] == 0 and assignment_test[9] == 16
        indices = [(f, s) for f in range(33) for s in range(33-f)]
        fd = np.asarray([f for f, _ in indices], dtype=int)
        sd = np.asarray([s for _, s in indices], dtype=int)
        local = np.asarray([coefficients[f, s] for f, s in indices])
        regions = [None]*257
        stats = {"max_imaginary_over_DC": 0.0, "max_negative_over_DC": 0.0,
                 "min_selected_over_DC": 1.0}
        center_reports = []
        for center in centers:
            guard()
            u, v = screen.marker_curve(center, curve_full, curve_single)
            lu = math.log(u.numerator)-math.log(u.denominator)
            lv = math.log(v.numerator)-math.log(v.denominator)
            mass = local*np.exp(fd*lu+sd*lv)[:, None]
            degree = np.zeros((MODULUS, MODULUS, 6))
            np.add.at(degree, (fd % MODULUS, sd % MODULUS), mass)
            p = degree[:, :, :4].reshape(MODULUS, MODULUS, 2, 2)
            t = degree[:, :, 4:].reshape(MODULUS, MODULUS, 2, 1)
            evaluated = np.sum(mass, axis=0)
            assert np.allclose(np.sum(p, axis=(0, 1)).ravel(), evaluated[:4], rtol=2e-13, atol=0)
            assert np.allclose(np.sum(t, axis=(0, 1)).ravel(), evaluated[4:], rtol=2e-13, atol=0)
            ordinary_p = screen.normalized_matrix(tuple(float(x) for x in evaluated[:4]))
            ordinary_t = screen.normalized_matrix(tuple(float(x) for x in evaluated[4:]))
            ordinary_b = screen.matrix_power(ordinary_p, 128)
            ordinary_v = screen.matrix_vector(screen.matrix_power(ordinary_p, 127), ordinary_t)
            spectra = fft.regional_spectra(p, t, np, guard)
            tables = []
            for kind, spectrum, expected in zip(("regular", "terminal"), spectra, (ordinary_b, ordinary_v)):
                raw, scale = spectrum
                try:
                    fft.scaled_check(tuple(float(x.real) for x in raw[0, 0].ravel()), scale, expected)
                    residue, _, _, diagnostic = fft.recover(spectrum, np)
                except ArithmeticError as error:
                    raise ArithmeticError(f"center={center},modulus={MODULUS},object={kind}: {error}") from error
                thresholds = selection_thresholds(fft, residue, raw, np)
                tables.append((residue, scale, kind, thresholds))
                for key in ("max_imaginary_over_DC", "max_negative_over_DC"):
                    stats[key] = max(stats[key], diagnostic[key])
            for full in assigned[center]:
                guard()
                single = screen.WEIGHT-2*full
                target = (full % MODULUS, single % MODULUS)
                projected = []
                for residue, scale, kind, thresholds in tables:
                    try:
                        values, smallest = select_entries(residue, target, thresholds)
                    except ArithmeticError as error:
                        raise ArithmeticError(f"center={center},full={full},single={single},"
                                              f"modulus={MODULUS},object={kind}: {error}") from error
                    stats["min_selected_over_DC"] = min(stats["min_selected_over_DC"], smallest)
                    projected.append(screen.normalized_matrix(values, scale))
                denominator = (math.log(comb(screen.PAIRS, full)*comb(screen.PAIRS-full, single))
                               +full*lu+single*lv)
                assert regions[full] is None
                regions[full] = (projected[0], projected[1], denominator, u, v)
            report = {"center": center, "u": str(u), "v": str(v),
                "first_full": min(assigned[center]), "last_full": max(assigned[center]),
                "seconds_in_marginal_phase": time.monotonic()-started}
            center_reports.append(report)
            print(json.dumps({"phase": "marker_cell_complete", **report}), flush=True)
            del tables, spectra, degree, p, t
        assert all(region is not None for region in regions)
        selected = regions[223]
        chain, reverse = screen.paired_chains(selected, selected)
        assert abs(chain-reverse) < 1e-8
        guard()
        counts = screen.load("cheap_lane_fast_profile_counts.py", screen.COUNT_PIN, "cell_selected_count")
        number = counts.profile_count(4096, 512, 223, 223)
        eta = (1+F(512+comb(512, 2), 255))/(255*F(1, 3))
        source = float(alpha)*(math.log(eta.numerator)-math.log(eta.denominator)
            +512*(160*math.log(256)-32*math.log(255))-734003*math.log(float(z)))
        metadata.update({"center_reports": center_reports,
            "controls": {"exact_scale_and_alias_tests": scale_tests, "exact_FFT_toy_coefficients": toy,
                "nearest_center_tie_test": True,
                "all_257_marginals_retained": True, "bank_packet_and_regional_DC_replays": True,
                "one_inverse_per_regional_matrix_or_terminal_per_cell": True, **stats},
            "selected_f223": {"assigned_center": assignment[223], "u": str(selected[3]),
                "v": str(selected[4]), "z": str(z), "alpha": str(alpha),
                "log2_conditional_chain_FLOAT": chain/math.log(2),
                "log2_weighted_profile_FLOAT": (chain+math.log(number)+source)/math.log(2),
                "directed_replay_required_at_these_marker_values": True}})
        guard()
        print(json.dumps({"phase": "marker_cell_marginals_ready", **metadata}), flush=True)
        return regions

    screen.build_regions = build_regions


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--centers", default=DEFAULT_CENTERS)
    parser.add_argument("--z", default="17447/25000")
    parser.add_argument("--alpha", default="119/1280")
    parser.add_argument("--curve-full", default="1562175/990029")
    parser.add_argument("--curve-single", default="1140210/880327")
    args = parser.parse_args()
    try:
        centers = tuple(sorted(int(x) for x in args.centers.split(",")))
    except ValueError:
        parser.error("centers must be a comma-separated list of integers")
    if not 1 <= len(centers) <= 16 or len(set(centers)) != len(centers) or not 0 <= centers[0] <= centers[-1] <= 256:
        parser.error("require1..16 distinct centers in0..256")
    z, alpha = F(args.z), F(args.alpha)
    fft = load_fft()
    screen = fft.load_screen()
    metadata = {"method": "FLOAT fixed nearest-marker cells; only build_regions replaced",
        "directed": False, "negative_is_not_a_certificate": True, "modulus": MODULUS,
        "centers": centers, "nearest_tie_rule": "smaller center", "fresh_marker_calibration": False,
        "fallback_or_skipping": False, "numerical_policy": {"negative_clipping": False,
            "frequency_dependent_scales": False, "per_entry_DC_relative_tolerance": fft.RELATIVE_TOLERANCE,
            "selected_relative_floor": fft.SELECTED_FLOOR, "selected_noise_multiplier": fft.NOISE_MULTIPLIER},
        "FFT_pin": FFT_PIN, "screen_pin": fft.SCREEN_PIN, "frozen_main_source_sha256": fft.SCREEN_PIN,
        "wrapper_sha256": hashlib.sha256(Path(__file__).read_bytes()).hexdigest()}
    install(screen, fft, centers, z, alpha, metadata)
    print(json.dumps({"phase": "cell_wrapper_start", **metadata}), flush=True)
    sys.argv = [sys.argv[0], "--z", args.z, "--alpha", args.alpha,
                "--curve-full", args.curve_full, "--curve-single", args.curve_single]
    screen.main()
    print(json.dumps({"phase": "cell_wrapper_complete", **metadata}, indent=2), flush=True)


if __name__ == "__main__":
    main()
