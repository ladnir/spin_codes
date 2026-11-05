"""FLOAT all-r512 residue screen; replace only a pinned build_regions seam.

The frozen screen supplies local cap^alpha coefficients, exact indexed
counts, finite scaled2x2 chains, both terminal orientations, and complete
source-count checks. Here each regional numerator is replaced by its
positive bivariate residue projection. FFT2 acts only on route-degree axes.
Seven pointwise matrix doublings yield P^128 and P^127*T. ONE real log scale
per entire frequency array is retained, never a different scale per bin.

This is a numerical diagnostic, NOT a certificate. Negative/imaginary
residues are checked relative to each entry's own DC mass. No clipping is
performed, and unreliable selected coefficients cause a hard failure.
Default modulus64; optional128 is one fixed run, not an automatic scan.
The inherited fixed rational marker
curve/z/alpha options are unchanged.180s/256MiB/20GiB, stdout only.
"""

import argparse
from collections import defaultdict
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
SCREEN_PIN = "6cd1f551da748c461c995f7927144345655cc741156d0550625388fcd6876861"
RELATIVE_TOLERANCE = 5e-11
SELECTED_FLOOR = 1e-10
NOISE_MULTIPLIER = 100.0


def load_screen():
    path = Path(__file__).with_name("cheap_lane_cauchy_profile_screen.py")
    if hashlib.sha256(path.read_bytes()).hexdigest() != SCREEN_PIN:
        raise RuntimeError("frozen profile screen hash mismatch")
    spec = importlib.util.spec_from_file_location("residue_profile_screen", path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def scaled_spectrum(values, scale, np):
    """One scale shared by every frequency, preserving the inverse FFT."""
    dc = values[0, 0]
    largest = float(np.max(dc.real))
    if not largest > 0 or not np.isfinite(values).all():
        raise ArithmeticError("invalid FFT matrix/DC mass")
    if float(np.max(np.abs(dc.imag))) > RELATIVE_TOLERANCE*largest:
        raise ArithmeticError("material imaginary zero-frequency entry")
    if float(np.min(dc.real)) < -RELATIVE_TOLERANCE*largest:
        raise ArithmeticError("negative zero-frequency entry")
    return values/largest, scale+math.log(largest)


def pointwise_product(a, b, np):
    """Explicit state product; route frequencies remain independent."""
    out = np.empty((a.shape[0], a.shape[1], 2, b.shape[3]), dtype=complex)
    for j in range(b.shape[3]):
        out[:, :, 0, j] = a[:, :, 0, 0]*b[:, :, 0, j]+a[:, :, 0, 1]*b[:, :, 1, j]
        out[:, :, 1, j] = a[:, :, 1, 0]*b[:, :, 0, j]+a[:, :, 1, 1]*b[:, :, 1, j]
    return out


def regional_spectra(p, t, np, guard, doublings=7):
    ph, plog = scaled_spectrum(np.fft.fft2(p, axes=(0, 1)), 0.0, np)
    th, tlog = scaled_spectrum(np.fft.fft2(t, axes=(0, 1)), 0.0, np)
    for _ in range(doublings):
        guard()
        th, tlog = scaled_spectrum(pointwise_product(ph, th, np), plog+tlog, np)
        ph, plog = scaled_spectrum(pointwise_product(ph, ph, np), 2*plog, np)
    return (ph, plog), (th, tlog)


def recover(spectrum, np, target=None):
    """Check all residues, then require each selected entry to be resolved."""
    values, logscale = spectrum
    residue = np.fft.ifft2(values, axes=(0, 1))
    selected = []
    max_imaginary = max_negative = 0.0
    min_selected = 1.0
    for i in range(2):
        for j in range(values.shape[3]):
            dc = float(values[0, 0, i, j].real)
            if not dc > 0:
                raise ArithmeticError("nonpositive matrix-entry DC mass")
            block = residue[:, :, i, j]
            imaginary = float(np.max(np.abs(block.imag)))/dc
            negative = max(0.0, -float(np.min(block.real)))/dc
            max_imaginary = max(max_imaginary, imaginary)
            max_negative = max(max_negative, negative)
            if imaginary > RELATIVE_TOLERANCE or negative > RELATIVE_TOLERANCE:
                raise ArithmeticError("material imaginary/negative FFT residue")
            if abs(float(np.sum(block.real))-dc) > RELATIVE_TOLERANCE*dc:
                raise ArithmeticError("sum of residues does not replay DC mass")
            if target is not None:
                x = float(block[target].real)
                relative = x/dc
                # Numerical-noise indicators are diagnostics, not rigorous
                # error bounds. A negative screen still needs directed replay.
                noise = max(imaginary, negative, 32*np.finfo(float).eps)
                if relative <= max(SELECTED_FLOOR, NOISE_MULTIPLIER*noise):
                    raise ArithmeticError("selected FFT coefficient is numerically unresolved")
                selected.append(x)
                min_selected = min(min_selected, relative)
    return residue, tuple(selected), logscale, {
        "max_imaginary_over_DC": max_imaginary,
        "max_negative_over_DC": max_negative,
        "min_selected_over_DC": min_selected}


def scaled_check(actual, scale, expected):
    """Every entry, including small channels, matches the ordinary matrix."""
    values, reference_scale = expected
    for x, y in zip(actual, values):
        if not x > 0 or not y > 0:
            raise ArithmeticError("nonpositive ordinary/DC replay entry")
        if abs(math.log(x)+scale-math.log(y)-reference_scale) > 5e-9:
            raise ArithmeticError("FFT zero-frequency matrix disagrees with ordinary power")


def exact_toy(np, guard):
    """Small rational polynomial matrices versus both FFT regional objects."""
    p = [[{(0, 0): F(1, 3), (2, 1): F(2, 5)}, {(1, 2): F(3, 7)}],
         [{(3, 0): F(5, 11)}, {(0, 1): F(7, 13), (1, 3): F(11, 17)}]]
    t = [[{(4, 2): F(13, 19)}], [{(0, 3): F(17, 23), (2, 1): F(19, 29)}]]

    def product(a, b):
        out = []
        for row in a:
            result = []
            for j in range(len(b[0])):
                item = defaultdict(F)
                for k, left in enumerate(row):
                    for (f, s), mass in left.items():
                        for (ff, ss), other in b[k][j].items():
                            item[f+ff, s+ss] += mass*other
                result.append(dict(item))
            out.append(result)
        return out

    p2 = product(p, p)
    expected = (product(p2, p2), product(product(p2, p), t))
    checks = 0
    for m in (4, 8):
        def pack(table):
            out = np.zeros((m, m, 2, len(table[0])))
            for i, row in enumerate(table):
                for j, item in enumerate(row):
                    for (f, s), mass in item.items():
                        out[f % m, s % m, i, j] += float(mass)
            return out

        spectra = regional_spectra(pack(p), pack(t), np, guard, doublings=2)
        for spectrum, truth in zip(spectra, expected):
            residue, _, scale, _ = recover(spectrum, np)
            exact = np.zeros(residue.shape)
            for i, row in enumerate(truth):
                for j, item in enumerate(row):
                    for (f, s), mass in item.items():
                        exact[f % m, s % m, i, j] += float(mass)
            actual = residue.real*math.exp(scale)
            error = float(np.max(np.abs(actual-exact)))
            assert error <= 2e-12*max(1.0, float(np.max(exact)))
            checks += exact.size
    return checks


def install_seam(screen, modulus, z, alpha, metadata):
    def build_regions(coefficients, curve_full, curve_single, np, guard):
        toy_checks = exact_toy(np, guard)
        indices = [(f, s) for f in range(33) for s in range(33-f)]
        fd = np.asarray([f for f, _ in indices], dtype=int)
        sd = np.asarray([s for _, s in indices], dtype=int)
        local = np.asarray([coefficients[f, s] for f, s in indices])
        regions = []
        diagnostics = {"max_imaginary_over_DC": 0.0, "max_negative_over_DC": 0.0,
                       "min_selected_over_DC": 1.0}
        for full in range(257):
            guard()
            single = screen.WEIGHT-2*full
            u, v = screen.marker_curve(full, curve_full, curve_single)
            lu = math.log(u.numerator)-math.log(u.denominator)
            lv = math.log(v.numerator)-math.log(v.denominator)
            mass = local*np.exp(fd*lu+sd*lv)[:, None]
            degree = np.zeros((modulus, modulus, 6))
            np.add.at(degree, (fd % modulus, sd % modulus), mass)
            p = degree[:, :, :4].reshape(modulus, modulus, 2, 2)
            t = degree[:, :, 4:].reshape(modulus, modulus, 2, 1)
            ordinary_entries = np.sum(mass, axis=0)
            assert np.allclose(np.sum(p, axis=(0, 1)).ravel(), ordinary_entries[:4], rtol=2e-13, atol=0)
            assert np.allclose(np.sum(t, axis=(0, 1)).ravel(), ordinary_entries[4:], rtol=2e-13, atol=0)
            ordinary_p = screen.normalized_matrix(tuple(float(x) for x in ordinary_entries[:4]))
            ordinary_t = screen.normalized_matrix(tuple(float(x) for x in ordinary_entries[4:]))
            ordinary_b = screen.matrix_power(ordinary_p, 128)
            ordinary_v = screen.matrix_vector(screen.matrix_power(ordinary_p, 127), ordinary_t)
            spectra = regional_spectra(p, t, np, guard)
            projected = []
            target = (full % modulus, single % modulus)
            for kind, spectrum, expected in zip(("regular", "terminal"), spectra, (ordinary_b, ordinary_v)):
                raw, scale = spectrum
                try:
                    scaled_check(tuple(float(x.real) for x in raw[0, 0].ravel()), scale, expected)
                    _, selected, scale, stats = recover(spectrum, np, target)
                except ArithmeticError as error:
                    raise ArithmeticError(f"full={full},single={single},modulus={modulus},"
                                          f"object={kind}: {error}") from error
                projected.append(screen.normalized_matrix(selected, scale))
                diagnostics["max_imaginary_over_DC"] = max(diagnostics["max_imaginary_over_DC"], stats["max_imaginary_over_DC"])
                diagnostics["max_negative_over_DC"] = max(diagnostics["max_negative_over_DC"], stats["max_negative_over_DC"])
                diagnostics["min_selected_over_DC"] = min(diagnostics["min_selected_over_DC"], stats["min_selected_over_DC"])
            denominator = (math.log(comb(screen.PAIRS, full)*comb(screen.PAIRS-full, single))
                           +full*lu+single*lv)
            regions.append((projected[0], projected[1], denominator, u, v))
            if full % 64 == 0:
                print(json.dumps({"phase": "FFT_residue_marginals", "full_completed": full,
                                  "modulus": modulus}), flush=True)
        selected = regions[223]
        chain, reverse = screen.paired_chains(selected, selected)
        assert abs(chain-reverse) < 1e-8
        counts = screen.load("cheap_lane_fast_profile_counts.py", screen.COUNT_PIN, "residue_selected_count")
        guard()
        number = counts.profile_count(4096, 512, 223, 223)
        guard()
        eta = (1+F(512+comb(512, 2), 255))/(255*F(1, 3))
        source = float(alpha)*(math.log(eta.numerator)-math.log(eta.denominator)
            +512*(160*math.log(256)-32*math.log(255))-734003*math.log(float(z)))
        weighted = (chain+math.log(number)+source)/math.log(2)
        replay = (modulus == 64 and z == F(17447, 25000) and alpha == F(119, 1280)
                  and selected[3] == F(38707225, 418782267) and selected[4] == F(929060, 41375369))
        if replay:
            # Parent's directed residue64 result at the identical curve markers.
            assert abs(weighted-230.9064716006121906391795875975) < 1e-4
        metadata.update({"FFT_controls": {"exact_toy_coefficients": toy_checks,
            "all_257_packet_DC_and_regional_DC_checks": True,
            "sum_residues_checks": True, **diagnostics},
            "selected_f223_crosscheck": {"full": 223, "single": 66,
                "u": str(selected[3]), "v": str(selected[4]), "z": str(z), "alpha": str(alpha),
                "log2_conditional_chain_FLOAT": chain/math.log(2),
                "indexed_source_count_exact": number,
                "log2_weighted_profile_FLOAT": weighted,
                "directed_curve_marker_reference_checked": replay,
                "directed_reference_log2": 230.9064716006121906391795875975 if replay else None,
                "reference_scope": ("recorded residue64 at identical rational curve markers, not optimized markers"
                    if replay else f"no stored reference for these parameters/modulus{modulus}; directed replay required")}})
        print(json.dumps({"phase": "FFT_selected_profile_crosscheck", **metadata}), flush=True)
        return regions

    screen.build_regions = build_regions


def main():
    parser = argparse.ArgumentParser(add_help=False)
    parser.add_argument("--modulus", type=int, choices=(8, 16, 32, 64, 128), default=64)
    own, remaining = parser.parse_known_args()
    witnesses = argparse.ArgumentParser(add_help=False)
    witnesses.add_argument("--z", default="7/10")
    witnesses.add_argument("--alpha", default="119/1280")
    parsed, _ = witnesses.parse_known_args(remaining)
    metadata = {"wrapper_method": "FLOAT two-axis residue FFT2; only build_regions replaced",
        "modulus": own.modulus, "directed": False, "negative_is_not_a_certificate": True,
        "numerical_policy": {"per_entry_DC_relative_tolerance": RELATIVE_TOLERANCE,
            "selected_relative_floor": SELECTED_FLOOR, "selected_noise_multiplier": NOISE_MULTIPLIER,
            "negative_residue_clipping": False, "frequency_dependent_scales": False},
        "screen_pin": SCREEN_PIN,
        "wrapper_sha256": hashlib.sha256(Path(__file__).read_bytes()).hexdigest()}
    screen = load_screen()
    install_seam(screen, own.modulus, F(parsed.z), F(parsed.alpha), metadata)
    print(json.dumps({"phase": "wrapper_start", **metadata}), flush=True)
    sys.argv = [sys.argv[0]]+remaining
    screen.main()
    print(json.dumps({"phase": "wrapper_complete", **metadata}, indent=2), flush=True)


if __name__ == "__main__":
    main()
