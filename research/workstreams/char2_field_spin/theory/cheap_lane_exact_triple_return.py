"""Exact three-input return moments, followed by a selected FLOAT replay.

Import the frozen actual AES-GF256 matrices and one/two-input refinements.
For each of the560 input triples I, B=S[:,I] has rank3 and C_I is injective.
Each nonzero x has an exact zero-row set F, which is a proper row-span flat.
Enumerate the rank1 flats from normalized rows and the rank2 flats from
cross-products of row pairs. Their closures are16-bit row masks.

A rank2 flat has q-1 vectors with exactly that zero set. A rank1 flat has
q^2-1 minus the counts of containing rank2 flats. The empty flat receives
the remaining q^3-1 vectors. Thus every histogram coefficient is an exact
nonnegative integer; no field-vector enumeration or cancellation is used
to evaluate the resulting positive polynomial. One complete projective
enumeration independently checks a selected actual triple.

Take the uniform maximum of q^-3 sum_w N_w z^w over all I, intersect with
the old D3, and retain the frozen D1/D2 improvements. Only raw return E
changes, before alpha. All other channels, paired route dependencies,
source factors, exact indexed counts and the output-only terminal remain.
Both replays are FLOAT diagnostics for the same two selected r1024 profiles,
not directed or whole-size certificates.180s/256MiB/20GiB; stdout only.
"""

import argparse
from fractions import Fraction as F
import hashlib
import importlib.util
from itertools import combinations, combinations_with_replacement
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
SMALL_PIN = "62eff6f2f831647696f451fdc3e736e63de0484298781725b2ac3fe89f132a3d"
Q, WEIGHT, MODULUS = 256, 1024, 512
PROFILES = ((384, 500), (384, 384))
SECONDS, MEMORY, LN2 = 180, 256*(1 << 20), math.log(2)


def load_small():
    path = Path(__file__).with_name("cheap_lane_exact_small_return.py")
    if hashlib.sha256(path.read_bytes()).hexdigest() != SMALL_PIN:
        raise RuntimeError("frozen small-return helper hash mismatch")
    spec = importlib.util.spec_from_file_location("triple_return_small", path)
    module = importlib.util.module_from_spec(spec); spec.loader.exec_module(module)
    return module


def triple_histogram(rows, mul, inverses):
    """Exact nonzero-vector histogram of an injective16-by3 matrix."""
    assert len(rows) == 16 and all(len(row) == 3 and any(row) for row in rows)

    def canonical(row):
        pivot = next(value for value in row if value)
        inverse = inverses[pivot]
        return tuple(mul[value][inverse] for value in row)

    rank_one = {}
    for i, row in enumerate(rows):
        key = canonical(row)
        rank_one[key] = rank_one.get(key, 0) | (1 << i)
    assert sum(mask.bit_count() for mask in rank_one.values()) == 16
    planes, closures = {}, {}
    for i, j in combinations(range(16), 2):
        a, b = rows[i], rows[j]
        normal = (mul[a[1]][b[2]] ^ mul[a[2]][b[1]],
                  mul[a[2]][b[0]] ^ mul[a[0]][b[2]],
                  mul[a[0]][b[1]] ^ mul[a[1]][b[0]])
        if not any(normal):
            assert canonical(a) == canonical(b)
            continue
        normal = canonical(normal)
        if normal not in planes:
            mask = sum(1 << k for k, row in enumerate(rows)
                       if not (mul[normal[0]][row[0]] ^ mul[normal[1]][row[1]]
                               ^ mul[normal[2]][row[2]]))
            assert mask.bit_count() >= 2 and mask != (1 << 16)-1
            # A closure contains its generating basis, so it identifies its span.
            assert mask not in closures
            planes[normal] = mask; closures[mask] = normal
        assert planes[normal] & (1 << i) and planes[normal] & (1 << j)
    assert planes and len(planes) <= 120
    one_masks, two_masks = tuple(rank_one.values()), tuple(planes.values())
    assert len(set(one_masks)) == len(one_masks)
    assert not set(one_masks).intersection(two_masks)
    histogram = [0]*17
    for mask in two_masks:
        histogram[16-mask.bit_count()] += Q-1
    for mask in one_masks:
        containing = sum(mask & other == mask for other in two_masks)
        number = Q**2-1-containing*(Q-1)
        assert number >= 0
        histogram[16-mask.bit_count()] += number
    empty = Q**3-1-sum(histogram)
    assert empty >= 0
    histogram[16] += empty
    assert histogram[0] == 0 and sum(histogram) == Q**3-1
    assert all(number >= 0 and number % (Q-1) == 0 for number in histogram)
    # Independent first/second zero-incidence identities for every input triple.
    assert sum((16-w)*number for w, number in enumerate(histogram)) == 16*(Q**2-1)
    parallel = sum(math.comb(mask.bit_count(), 2) for mask in one_masks)
    second = parallel*(Q**2-1)+(120-parallel)*(Q-1)
    assert sum(math.comb(16-w, 2)*number for w, number in enumerate(histogram)) == second
    return tuple(histogram), (len(one_masks), len(two_masks))


def triples(c, s, inverses, small, guard):
    """Bounded low-rank geometry; the multiplication table occupies64KiB."""
    mul = []
    for a in range(Q):
        guard(); mul.append(bytes(small.gf_mul(a, b) for b in range(Q)))
    records, flat_counts = [], []
    for points in combinations(range(16), 3):
        guard()
        assert small.gf_rank([[row[i] for i in points] for row in c], inverses) == 3
        rows = [tuple(row[i] for i in points) for row in s]
        assert small.gf_rank(rows, inverses) == 3
        histogram, counts = triple_histogram(rows, mul, inverses)
        records.append((points, histogram)); flat_counts.append(counts)
    assert len(records) == 560
    # Independent full projective enumeration, using(1,b,c),(0,1,c),(0,0,1).
    points = (0, 1, 2)
    rows = [tuple(row[i] for i in points) for row in s]
    direct = [0]*17
    for b in range(Q):
        guard()
        prefix = tuple(row[0] ^ mul[b][row[1]] for row in rows)
        for cc in range(Q):
            weight = sum(bool(value ^ mul[cc][row[2]]) for value, row in zip(prefix, rows))
            direct[weight] += 1
    for cc in range(Q):
        weight = sum(bool(row[1] ^ mul[cc][row[2]]) for row in rows)
        direct[weight] += 1
    direct[sum(bool(row[2]) for row in rows)] += 1
    assert sum(direct) == Q**2+Q+1
    assert tuple((Q-1)*number for number in direct) == dict(records)[points]
    del mul
    return records, {"actual_input_triples": len(records), "independent_C_and_S_triples": 560,
        "closure_uniqueness_checks": 560, "zero_incidence_checks": 1120,
        "maximum_rank1_flats": max(x[0] for x in flat_counts),
        "maximum_rank2_flats": max(x[1] for x in flat_counts),
        "direct_projective_check_triple": points, "direct_projective_directions": Q**2+Q+1,
        "distinct_positive_histograms": len({histogram for _, histogram in records}),
        "histogram_sha256": hashlib.sha256(json.dumps(records, separators=(",", ":")).encode()).hexdigest()}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--z", type=F, default=F(304397809, 625000000))
    parser.add_argument("--alpha", type=F, default=F(119, 1280))
    args = parser.parse_args(); z, alpha = args.z, args.alpha
    if not 0 < z < 1 or not 0 < alpha <= 1:
        parser.error("require0<z<1 and0<alpha<=1")
    if max(n.bit_length() for x in (z, alpha) for n in (x.numerator, x.denominator)) > 64:
        parser.error("rational witnesses limited to64bits")
    if shutil.disk_usage(Path(__file__).resolve().anchor).free < 20*(1 << 30):
        raise RuntimeError("requires20GiB free on source volume")
    started = time.monotonic(); deadline = started+SECONDS
    small = load_small(); calibration = small.load_calibration(); broad = calibration.load_broad()
    neighbor = broad.load("cheap_lane_neighbor_screen.py", broad.NEIGHBOR_PIN, "triple_neighbor")
    cells = neighbor.load_cells(); fft = cells.load_fft(); screen = fft.load_screen()
    old = screen.load("cheap_lane_cauchy_marker_calibration.py", screen.CALIBRATION_PIN, "triple_calibration")
    gradient = old.load_gradient(); anchor = gradient.load_anchor(); peak = anchor.peak_memory_reader()

    def guard():
        if time.monotonic() > deadline: raise TimeoutError("180-second triple-return deadline; no retry")
        if peak() > MEMORY: raise MemoryError("256MiB triple-return budget exceeded")

    def progress(phase, **values):
        print(json.dumps({"phase": phase, "seconds": time.monotonic()-started, **values}), flush=True)

    a, c, s, inverses, weights = small.actual_matrices()
    pairs = small.pair_enumerators(c, s, inverses, guard)
    tests = small.exact_tests(c, s, inverses, weights, pairs, guard)
    records, triple_tests = triples(c, s, inverses, small, guard)
    tests["triple_tests"] = triple_tests
    tests["matrix_sha256"] = hashlib.sha256(bytes(value for row in s for value in row)).hexdigest()
    progress("exact_triple_histograms_ready", tests=tests)
    import numpy as np
    from flint import ctx
    ctx.prec = 256
    tests["FFT_exact_toy_coefficients"] = fft.exact_toy(np, guard)
    adapter = anchor.load_adapter(); core = adapter.load_verifier()
    core.Z, core.ALPHA = z, alpha

    class LocalVerifier(core.Verifier):
        def guard(self): guard()

    verifier = LocalVerifier(2, deadline)
    _, quiet, old_d, _, _ = core.rational_caps()
    baseline_d, block_report = small.refined_blocks(z, quiet, old_d, pairs)
    assert quiet[3] == F(1, Q**3)
    values = [(points, small.polynomial(histogram, z)/Q**3) for points, histogram in records]
    maximum = max(value for _, value in values)
    assert 0 <= maximum <= 1-F(1, Q**3)
    refined_d = list(baseline_d); refined_d[3] = min(baseline_d[3], maximum)
    block_report.update({"old_D3": str(baseline_d[3]), "new_D3": str(refined_d[3]),
        "exact_uniform_D3": str(maximum),
        "D3_maximizing_triples": [points for points, value in values if value == maximum]})
    assert all(0 <= new <= before for new, before in zip(refined_d, baseline_d))
    shapes = list(combinations_with_replacement(range(17), 4))
    caps = verifier.local_caps(shapes); exponent = verifier.exact(alpha)
    powered = {occ: tuple(verifier.positive(value)**exponent for value in row)
               for occ, row in caps.items()}

    def check_float_caps():
        for row in powered.values():
            for value in row:
                approximate = float(value.mid())
                if not math.isfinite(approximate) or approximate < 0 or (value > 0 and approximate == 0):
                    raise ArithmeticError("positive powered cap not representable in FLOAT")

    changed = 0; largest_ratio = F(1)
    for occ in shapes:
        guard()
        adapter.exact_point(caps[occ][2], small.return_cap(occ, quiet, old_d, core.KAPPA), verifier.fmpq)
        row = list(powered[occ])
        row[2] = verifier.exact(small.return_cap(occ, quiet, baseline_d, core.KAPPA))**exponent
        powered[occ] = tuple(row)
    del caps
    check_float_caps()
    baseline = old.float_coefficients(powered, np, guard)
    for occ in shapes:
        guard()
        before = small.return_cap(occ, quiet, baseline_d, core.KAPPA)
        after = small.return_cap(occ, quiet, refined_d, core.KAPPA)
        assert 0 <= after <= before
        if after < before:
            changed += 1
            if after: largest_ratio = max(largest_ratio, before/after)
        row = list(powered[occ]); row[2] = verifier.exact(after)**exponent
        powered[occ] = tuple(row)
    check_float_caps()
    refined = old.float_coefficients(powered, np, guard); del powered
    assert np.array_equal(baseline[:, :, [0, 1, 3, 4, 5]], refined[:, :, [0, 1, 3, 4, 5]])
    assert np.all(refined[:, :, 2] <= baseline[:, :, 2]*(1+2e-13))
    replays = {"ell1_2": calibration.Replay(baseline, WEIGHT, cells, fft, screen, np, guard),
               "ell1_2_3": calibration.Replay(refined, WEIGHT, cells, fft, screen, np, guard)}
    counts = screen.load("cheap_lane_fast_profile_counts.py", screen.COUNT_PIN, "triple_counts")
    counts.DEADLINE = deadline
    h, c_h, eta = anchor.source_factor(WEIGHT)
    log_source = float((exponent*(verifier.exact(eta).log()+WEIGHT*verifier.exact(core.BETA).log()
                               -core.LIMIT*verifier.exact(z).log())).mid())
    cf, cs = F(1562175, 990029), F(1140210, 880327)
    progress("triple_return_caps_ready", block_caps=block_report, changed_shapes=changed,
        max_raw_return_gain_bits=(math.log(largest_ratio.numerator)-math.log(largest_ratio.denominator))/LN2)
    reports = []
    for fulls in PROFILES:
        guard(); number = counts.profile_count(4096, WEIGHT, *fulls)
        assert number > 0
        markers = tuple(broad.marker_curve(WEIGHT, f, cf, cs) for f in fulls)
        report = {"fulls": fulls, "singles": [WEIGHT-2*f for f in fulls],
            "markers": [[str(x) for x in marker] for marker in markers],
            "log2_exact_count": math.log(number)/LN2, "comparison": {}}
        for label, replay in replays.items():
            entries = [replay.marginal(f, marker, MODULUS) for f, marker in zip(fulls, markers)]
            if any(entry["status"] != "resolved" for entry in entries):
                value = {"status": "unresolved", "reasons": [x.get("reason") for x in entries]}
            else:
                chain = screen.paired_chains(entries[0]["region"], entries[1]["region"])[0]
                value = {"status": "resolved", "log2_bound_FLOAT": (log_source+math.log(number)+chain)/LN2,
                    "log2_conditional_chain_FLOAT": chain/LN2,
                    "FFT_controls": [x["FFT_controls"] for x in entries]}
            report["comparison"][label] = value
        before, after = report["comparison"]["ell1_2"], report["comparison"]["ell1_2_3"]
        if before["status"] == after["status"] == "resolved":
            gain = before["log2_bound_FLOAT"]-after["log2_bound_FLOAT"]
            assert gain >= -2e-5
            report["incremental_ell3_gain_bits_FLOAT"] = gain
        reports.append(report); progress("selected_profile_complete", **report)
    guard()
    print(json.dumps({"scope": "exact ell3 return moments; selected-profile FLOAT replay only",
        "directed_global_bound": False, "negative_is_not_a_certificate": True,
        "whole_weight_certificate": False, "whole_code_certificate": False,
        "weight": WEIGHT, "z": str(z), "alpha": str(alpha), "modulus": MODULUS,
        "changed_channel": "raw E only, before alpha; baseline already retains exact ell1/2 caps",
        "uniform_triple_policy": "max over all560 actual input triples, no independent point averaging",
        "source_factor": {"h": h, "c_h_exact": str(c_h), "eta_exact": str(eta),
            "powered_source_log2_FLOAT": log_source/LN2},
        "tests": tests, "block_caps": block_report, "changed_shapes": changed, "results": reports,
        "seconds": time.monotonic()-started, "peak_process_bytes": peak(), "generated_bytes": 0,
        "small_return_pin": SMALL_PIN, "calibration_pin": small.CALIBRATION_PIN,
        "core_pin": anchor.CORE_PIN, "count_pin": screen.COUNT_PIN,
        "source_sha256": hashlib.sha256(Path(__file__).read_bytes()).hexdigest()}, indent=2), flush=True)


if __name__ == "__main__":
    main()
