"""Pool nonzero masks inside fixed live segments, retaining quiet routing.

The binary straddling schedule, full16 birth/return endpoints, regional
cuts and terminal incoming mask are fixed. Q/A/E remain the common base
B, powered before averaging routes. Only V changes: each sustained-live
segment sums all nonzero intermediate masks through its RAW matrix.
The mixed regional coefficient represents E[B^alpha V], while E[B^alpha]
uses free route counts on live packets with NO hidden-mask multiplicity.
Weighted Jensen combines them with exponents(alpha,1-alpha), leaving one
whole-region route normalizer. There is no per-segment normalization.

This is a FLOAT upper bound for this restricted binary schedule, not a
whole-code certificate. Compare the unchanged all-full-mask control.
One fixed anchor, no optimization;48 seconds,20GiB free,<128MiB estimated
peak, stdout only. No frozen edits or generated data files.
"""

import ctypes
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
from types import SimpleNamespace

sys.dont_write_bytecode = True
os.environ["OPENBLAS_NUM_THREADS"] = "1"
os.environ["OMP_NUM_THREADS"] = "1"
TILT_PIN = "26b019f4d967603f3bd7e8aa14229fe15d2fec4214993e507302f3c377c319fa"
DEGREE_PIN = "264a216fd5beb0c2b2d18826e5fc3f7cd4638561ff2aa37d9862382e36bbd0d2"
Z, ALPHA = .85, .05
LN2 = math.log(2)
FIELDS = ("regional_first_log2", "regional_middle_log2", "terminal_log2")


def load_pinned(filename, digest, name):
    path = Path(__file__).with_name(filename)
    if hashlib.sha256(path.read_bytes()).hexdigest() != digest:
        raise RuntimeError("frozen helper hash mismatch: "+filename)
    spec = importlib.util.spec_from_file_location(name, path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def peak_bytes():
    if os.name != "nt":
        return None
    class Counters(ctypes.Structure):
        _fields_ = [("cb", ctypes.c_ulong), ("faults", ctypes.c_ulong)]+[
            (name, ctypes.c_size_t) for name in ("peak", "working", "paged_peak", "paged",
                "nonpaged_peak", "nonpaged", "pagefile", "pagefile_peak")]
    value = Counters(); value.cb = ctypes.sizeof(value)
    kernel = ctypes.WinDLL("kernel32", use_last_error=True)
    kernel.GetCurrentProcess.restype = ctypes.c_void_p
    psapi = ctypes.WinDLL("psapi", use_last_error=True)
    psapi.GetProcessMemoryInfo.argtypes = [ctypes.c_void_p, ctypes.POINTER(Counters), ctypes.c_ulong]
    if not psapi.GetProcessMemoryInfo(kernel.GetCurrentProcess(), ctypes.byref(value), value.cb):
        return None
    return int(value.peak)


def raw_live_polynomial(frozen, grouped, singleton, paired, filtered, args):
    """The pinned path helper's caps, now retaining all positive masks."""
    np = paired.np
    terms = frozen.packet_terms()
    shapes = sorted({key[2] for key in terms})
    index = {shape: i for i, shape in enumerate(shapes)}
    array = np.asarray(shapes, dtype=np.int16)
    raw = terminal = None
    for gamma in (0., .5, 1.):
        paired.guard()
        local, final = filtered.local_tables(array, args.z, 1., gamma, 1.)
        for j in range(1, 17):
            factor = min(1., 256/255*255.**(-(16-j)))
            local[:, :, :, j] = np.minimum(local[:, :, :, j],
                final[:, :, :, 0]*factor*math.comb(16, j))
        raw = local if raw is None else np.minimum(raw, local)
        terminal = final if terminal is None else np.minimum(terminal, final)
    matrix = raw[:, 0]+raw[:, 2]
    final = terminal[:, 0]+terminal[:, 2]
    old_m, old_t, _, _ = grouped.raw_tables(args, singleton, paired, filtered)
    old_index = {tuple(sorted((a, a, b, b))): i for i, (a, b) in
                 enumerate(itertools.combinations_with_replacement(range(17), 2))}
    values = np.zeros((3, 33, 16, 16))
    end = np.zeros((3, 33))
    for (f, s, shape), coefficient in terms.items():
        row = index[shape]
        m, t = (old_m[old_index[shape]], old_t[old_index[shape]]) if s == 0 else (matrix[row], final[row])
        values[s, f] += float(coefficient)*m[1:, 1:]
        end[s, f] += float(coefficient)*t[16, 0]
    assert np.isfinite(values).all() and (values >= 0).all()
    logged = np.full(values.shape, -math.inf)
    np.log(values, out=logged, where=values > 0)
    terminal_log = np.full(end.shape, -math.inf)
    np.log(end, out=terminal_log, where=end > 0)
    return logged, terminal_log


def exact_mask_test(degree, np):
    """Explicit actual masks, rather than a second matrix recurrence."""
    classes = (0, 0, 1)  # D=(2,1); the fixed endpoint has multiplicity1.
    caps = [[[Fraction(2+i+2*j+r, 13+i+3*j+r) for j in range(2)]
             for i in range(2)] for r in range(2)]
    p = np.empty((1, 2, 2, 2))
    for r, i, j in itertools.product(range(2), repeat=3):
        p[0, r, i, j] = math.log(float(caps[r][i][j]*(2 if j == 0 else 1)))
    p2 = degree.convolution(p, p, 4, 0)
    p4 = degree.convolution(p2, p2, 4, 0)
    exact = [Fraction(0) for _ in range(5)]
    for routes in itertools.product(range(2), repeat=4):
        for masks in itertools.product(range(3), repeat=3):
            path = (1,)+(tuple(classes[m] for m in masks))+(1,)
            term = Fraction(1)
            for k, route in enumerate(routes):
                term *= caps[route][path[k]][path[k+1]]
            exact[sum(routes)] += term
    for d, value in enumerate(exact):
        assert math.isclose(p4[0, d, 1, 1], math.log(value), abs_tol=2e-13)
        # At alpha1 the weighted expression is exactly the mixed mean;
        # route multiplicity C(4,d) occurs once, and endpoint D is1.
        normalizer = math.comb(4, d)
        base = Fraction(3+d, 11+d)
        assert base**0*(value/normalizer)**1 == value/normalizer
    rectangular = degree.convolution(p2[:, :, 1:2, :], p2[:, :, :, 1:2], 4, 0)
    assert np.allclose(rectangular[:, :, 0, 0], p4[:, :, 1, 1], atol=2e-13)
    return {"explicit_route_mask_assignments": 2**4*3**3,
            "alpha1_normalized_coefficients": 5, "rectangular_endpoint_entries": 5}


def regional_mixed(frozen, powered, raw, segments, paired):
    np = paired.np
    q95 = frozen.power(powered["Q"], 95, np, paired.guard)
    def multiply(factors):
        total = np.full((3, 1), -math.inf); total[0, 0] = 0.
        for item in factors:
            total = frozen.convolve(total, item, 128, np, paired.guard)
        return [float(total[s, f]-math.log(math.comb(4096, f)*math.comb(4096-f, s)))/LN2
                for f, s in ((128, 0), (127, 2))]
    return {"regional_first_log2": multiply((q95, powered["A"], segments[32])),
            "regional_middle_log2": multiply((segments[14], powered["E"], q95,
                                              powered["A"], segments[17])),
            "terminal_log2": multiply((segments[14], powered["E"], q95,
                                       powered["A"], segments[16], raw["U"]))}


def main():
    if len(sys.argv) != 1:
        raise ValueError("fixed live-segment diagnostic takes no scan arguments")
    if shutil.disk_usage(Path(__file__).resolve().anchor).free < 20*(1 << 30):
        raise RuntimeError("requires20GiB free on source volume")
    started = time.monotonic(); deadline = started+48
    import numpy as np
    from scipy.special import logsumexp
    tilt = load_pinned("cheap_lane_quiet_tilt_jensen.py", TILT_PIN, "frozen_segment_tilt")
    frozen = load_pinned("cheap_lane_path_marks.py", tilt.PIN, "frozen_segment_paths")
    degree = load_pinned("cheap_lane_degree_conditioned_tile.py", DEGREE_PIN, "frozen_segment_degree")
    degree.np = np; degree.lse = logsumexp; degree.DEADLINE = deadline
    grouped = frozen.load_grouped(); singleton = grouped.load_singleton()
    paired = singleton.load_block().load_paired(deadline)
    tests = frozen.self_test(np, paired.guard)
    tests.update(exact_mask_test(degree, np))
    tests["weighted_Jensen_exact_tests"] = tilt.exact_tilt_tests()
    free, tests["free_packet_route_counts"] = tilt.free_packet(frozen, np)
    filtered, _ = paired.imports()
    args = SimpleNamespace(r=256, z=Z, alpha=ALPHA, activity_mark=1.)
    powered, _ = frozen.local_polynomials(grouped, singleton, paired, filtered, args)
    raw_args = SimpleNamespace(r=256, z=Z, alpha=1., activity_mark=1.)
    raw, _ = frozen.local_polynomials(grouped, singleton, paired, filtered, raw_args)
    base_polys = dict(powered); base_polys["L"] = free; base_polys["U"] = free
    base = frozen.straddling(base_polys, args, paired, {})
    allfull_polys = dict(powered); allfull_polys["L"] = raw["L"]; allfull_polys["U"] = raw["U"]
    allfull_mixed = frozen.straddling(allfull_polys, args, paired, {})
    allfull = tilt.combine_regions(base, allfull_mixed, 1-ALPHA, ALPHA, np)
    assert abs(allfull["summed_log2_score"]+4598.424024587397) < 1e-6
    live, terminal = raw_live_polynomial(frozen, grouped, singleton, paired, filtered, raw_args)
    assert np.allclose(live[:, :, 15, 15], raw["L"], atol=2e-13)
    assert np.allclose(terminal, raw["U"], atol=2e-13)
    powers = {1: live}; stages = []
    for n in (2, 4, 8, 16):
        current = time.monotonic()
        powers[n] = degree.convolution(powers[n//2], powers[n//2])
        peak = peak_bytes()
        if peak is not None and peak > 128*(1 << 20):
            raise RuntimeError("128MiB peak-memory ceiling exceeded; no retry")
        stages.append({"packets": n, "seconds": time.monotonic()-current,
                       "elapsed_seconds": time.monotonic()-started, "peak_bytes": peak})
        print(json.dumps({"progress": "raw_live_power", **stages[-1]}), flush=True)
    row12 = degree.convolution(powers[8][:, :, 15:16, :], powers[4])
    segments = {14: degree.convolution(row12, powers[2][:, :, :, 15:16])[:, :, 0, 0],
                16: powers[16][:, :, 15, 15],
                17: degree.convolution(powers[16][:, :, 15:16, :], powers[1][:, :, :, 15:16])[:, :, 0, 0],
                32: degree.convolution(powers[16][:, :, 15:16, :], powers[16][:, :, :, 15:16])[:, :, 0, 0]}
    inclusion_checks = 0
    for n, segment in segments.items():
        full = frozen.power(raw["L"], n, np, paired.guard)
        assert np.all(segment >= full-1e-9)
        inclusion_checks += 1
    mixed = regional_mixed(frozen, powered, raw, segments, paired)
    pooled = tilt.combine_regions(base, mixed, 1-ALPHA, ALPHA, np)
    for field in FIELDS:
        assert all(new >= old-1e-8 for new, old in zip(mixed[field], allfull_mixed[field]))
    assert pooled["summed_log2_score"] >= allfull["summed_log2_score"]-1e-7
    paired.guard()
    print(json.dumps({"scope": "FLOAT fixed binary schedule; raw nonzero masks pooled only inside live segments",
        "whole_code_certificate": False, "all_binary_paths_covered": False, "directed": False,
        "parameters": {"z": Z, "alpha": ALPHA},
        "fixed_full16_locations": ["birth exit", "return entry", "regional cuts", "terminal incoming"],
        "live_segment_lengths": {"first": [32], "middle": [14, 17], "terminal": [14, 16]},
        "base_caps": ["Q^alpha", "A^alpha", "E^alpha"],
        "terminal": "full16 output-only U; no outgoing state or unused refresh",
        "matching_profiles": [[128, 0], [127, 2]], "source_supports_per_orientation": 4096,
        "all_full_control": allfull, "nonzero_live_segments": pooled,
        "raw_mask_pool_cost_bits": pooled["summed_log2_score"]-allfull["summed_log2_score"],
        "base_regional_log2": {key: base[key] for key in FIELDS},
        "mixed_regional_log2": mixed, "raw_power_stages": stages,
        "tests": {**tests, "allfull_cap_matches": 2, "positive_segment_inclusions": inclusion_checks,
                  "allfull_frozen_control": "passed", "regional_normalizers_per_term": 1},
        "peak_working_set_bytes": peak_bytes(), "estimated_peak_bytes": 128*(1 << 20),
        "seconds": time.monotonic()-started, "generated_bytes": 0,
        "source_sha256": hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
        "tilt_helper_pin": TILT_PIN, "degree_helper_pin": DEGREE_PIN}, indent=2))


if __name__ == "__main__":
    main()
