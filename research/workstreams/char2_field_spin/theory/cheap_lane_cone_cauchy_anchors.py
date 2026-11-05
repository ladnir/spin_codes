"""Four fixed, directed large-support Cauchy diagnostics, not a size cover.

For a fixed actual matching profile (f,s), let P(u,v) be the packet
matrix obtained by raising each local cone cap to alpha BEFORE averaging
route placements. Let T(u,v) be the output-only terminal column. Positive
coefficients give the entrywise bounds

 R_f,s <= P(u,v)^128 / (u^f v^s M_4096(f,s)),
 V_f,s <= P(u,v)^127 T(u,v) / (u^f v^s M_4096(f,s)).

The two actual matching profiles coincide at each selected anchor. Their
finite chain is e0 R^255 V, with 32767 updates and one terminal emission.
The exact number of all c-run supports, c=r-2f=s, conservatively bounds
the number of supports of this ONE ordered matching profile. It does not
license a claim about other profiles at the same source weight.

The joint missing-group lemma supplies Omega_r=beta^r eta_r, where
 h=floor(r/256), c_h=1-r/(256(h+1)),
 eta_r=(1+sum_{j=1}^h C(r,j)/255)/(255*c_h).
Apply (Omega_r z^-734003)^alpha once. No field distribution is replaced.

All markers and witnesses are fixed rationals. No optimizer, regional
coefficient arrays, discarded route tails, output files, or new encoder
are used. Four directed evaluations have a common 60-second cooperative
guard, 128MiB process-memory guard on Windows, and 20GiB free-space guard.
Positive results reject these fixed witnesses, not the construction.
"""

import argparse
from collections import defaultdict
from fractions import Fraction as F
import hashlib
import importlib.util
from itertools import combinations_with_replacement, product
import json
from math import comb
import os
from pathlib import Path
import shutil
import sys
import time

sys.dont_write_bytecode = True
os.environ["OPENBLAS_NUM_THREADS"] = "1"
os.environ["OMP_NUM_THREADS"] = "1"
ADAPTER_PIN = "83e00b6f6d148c8058428295b4132d401df08cfccd20cadcda1c1ba82df7b478"
CORE_PIN = "7ccf0cdab823ef15fdd5edd6156d4b49d9a000f2ae28c37d5d43f4b6ea341f01"
ALPHA = F(119, 1280)
ANCHORS = (
    (512, 16, 480, F(1, 225), F(2, 15), F(7, 10)),
    (1024, 64, 896, F(1, 49), F(2, 7), F(1, 2)),
    (2048, 256, 1536, F(1, 9), F(2, 3), F(1, 4)),
    (4096, 1024, 2048, F(1), F(2), F(1, 16)),
)
MEMORY_LIMIT = 128*(1 << 20)


def peak_memory_reader():
    if os.name != "nt":
        return lambda: 0
    import ctypes
    from ctypes import wintypes

    class Counters(ctypes.Structure):
        _fields_ = [("cb", wintypes.DWORD), ("PageFaultCount", wintypes.DWORD),
                    ("PeakWorkingSetSize", ctypes.c_size_t),
                    ("WorkingSetSize", ctypes.c_size_t),
                    ("QuotaPeakPagedPoolUsage", ctypes.c_size_t),
                    ("QuotaPagedPoolUsage", ctypes.c_size_t),
                    ("QuotaPeakNonPagedPoolUsage", ctypes.c_size_t),
                    ("QuotaNonPagedPoolUsage", ctypes.c_size_t),
                    ("PagefileUsage", ctypes.c_size_t),
                    ("PeakPagefileUsage", ctypes.c_size_t)]

    kernel = ctypes.WinDLL("kernel32", use_last_error=True)
    psapi = ctypes.WinDLL("psapi", use_last_error=True)
    kernel.GetCurrentProcess.restype = wintypes.HANDLE
    psapi.GetProcessMemoryInfo.argtypes = [wintypes.HANDLE,
                                          ctypes.POINTER(Counters), wintypes.DWORD]
    psapi.GetProcessMemoryInfo.restype = wintypes.BOOL
    handle = kernel.GetCurrentProcess()

    def read():
        info = Counters(); info.cb = ctypes.sizeof(info)
        if not psapi.GetProcessMemoryInfo(handle, ctypes.byref(info), info.cb):
            raise OSError(ctypes.get_last_error(), "GetProcessMemoryInfo failed")
        return int(info.PeakWorkingSetSize)
    return read


def load_adapter():
    path = Path(__file__).with_name("cheap_lane_cone_interval_tests.py")
    if hashlib.sha256(path.read_bytes()).hexdigest() != ADAPTER_PIN:
        raise RuntimeError("frozen interval-adapter hash mismatch")
    spec = importlib.util.spec_from_file_location("cauchy_anchor_adapter", path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    if module.PIN != CORE_PIN:
        raise RuntimeError("unexpected frozen core pin")
    return module


def half_weights(u, v, slots=16):
    """Ordered block occupancies; each singleton orientation has weight v/2."""
    out = defaultdict(F)
    for f in range(slots+1):
        for s in range(slots-f+1):
            weight = F(comb(slots, f)*comb(slots-f, s), 2**s)*u**f*v**s
            for left in range(s+1):
                out[f+left, f+s-left] += comb(s, left)*weight
    return dict(out)


def shape_weights(u, v, guard):
    half = half_weights(u, v)
    assert sum(half.values(), F(0)) == (1+u+v)**16
    # Local caps are invariant under permutations of their four blocks.
    unordered = [(a, b, value*(1 if a == b else 2))
                 for (a, b), value in half.items() if a <= b]
    out = defaultdict(F)
    for i, (a, b, value) in enumerate(unordered):
        guard()
        for j in range(i, len(unordered)):
            c, d, other = unordered[j]
            out[tuple(sorted((a, b, c, d)))] += value*other*(1 if i == j else 2)
    assert sum(out.values(), F(0)) == (1+u+v)**32
    return dict(out)


def source_factor(r):
    h = r//256
    c = 1-F(r, 256*(h+1))
    assert 0 <= h < r and h <= 256 and c > 0
    eta = (1+F(sum(comb(r, j) for j in range(1, h+1)), 255))/(255*c)
    return h, c, eta


def run_total(pairs, weight, runs):
    assert 0 < weight < 2*pairs and 1 <= runs <= min(weight, 2*pairs-weight)
    numerator = 2*pairs*comb(weight-1, runs-1)*comb(2*pairs-weight-1, runs-1)
    assert numerator % runs == 0
    return numerator//runs


def tests(adapter, verifier):
    checks = 0
    for slots in (1, 2):
        u, v = F(2, 7), F(3, 11)
        direct = defaultdict(F)
        categories = ((0, 0, F(1)), (1, 1, u), (1, 0, v/2), (0, 1, v/2))
        for choices in product(categories, repeat=slots):
            a = sum(item[0] for item in choices)
            b = sum(item[1] for item in choices)
            mass = F(1)
            for item in choices:
                mass *= item[2]
            direct[a, b] += mass
        exact = half_weights(u, v, slots)
        assert dict(direct) == exact
        checks += len(exact)
    # Noncommuting positive matrices exercise the exact final-packet order.
    p = [[F(1, 3), F(2, 5)], [F(3, 7), F(4, 11)]]
    t = [[F(5, 13)], [F(6, 17)]]
    mm = lambda a, b: [[sum((a[i][k]*b[k][j] for k in range(len(b))), F(0))
                        for j in range(len(b[0]))] for i in range(len(a))]
    expected = mm(mm(mm(p, p), p), t)
    pa = [[verifier.exact(x) for x in row] for row in p]
    ta = [[verifier.exact(x) for x in row] for row in t]
    actual = verifier.matmul(verifier.matpow(pa, 3), ta)
    for i in range(2):
        adapter.exact_point(actual[i][0], expected[i][0], verifier.fmpq)
        checks += 1
    for pairs in (3, 4):
        totals = defaultdict(int)
        for word in product((0, 1), repeat=2*pairs):
            weight = sum(word)
            if 0 < weight < 2*pairs:
                runs = sum(word[i] and not word[i-1] for i in range(2*pairs))
                totals[weight, runs] += 1
        for (weight, runs), value in totals.items():
            assert run_total(pairs, weight, runs) == value
            checks += 1
    return checks


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--precision", type=int, choices=(256, 384, 512), default=256)
    args = parser.parse_args()
    if shutil.disk_usage(Path(__file__).resolve().anchor).free < 20*(1 << 30):
        raise RuntimeError("requires20GiB free on the source volume")
    started = time.monotonic(); deadline = started+60
    peak_memory = peak_memory_reader()

    def guard():
        if time.monotonic() > deadline:
            raise TimeoutError("60-second fixed-anchor budget exceeded; no automatic retry")
        if peak_memory() > MEMORY_LIMIT:
            raise MemoryError("128MiB process-memory budget exceeded")

    adapter = load_adapter()
    core = adapter.load_verifier()
    from flint import arb, ctx
    ctx.prec = args.precision

    class CauchyVerifier(core.Verifier):
        def guard(self):
            guard()

    verifier = CauchyVerifier(2, deadline)
    check_count = tests(adapter, verifier)
    shapes = list(combinations_with_replacement(range(17), 4))
    reports = []
    logtwo = arb(2).log()
    for r, f, s, u, v, z in ANCHORS:
        guard()
        assert 2*f+s == r and s == r-2*f and f+s <= 4096
        assert u == F(f, 4096-f-s) and v == F(s, 4096-f-s)
        core.Z, core.ALPHA = z, ALPHA
        weights = shape_weights(u, v, guard)
        assert set(weights) == set(shapes)
        caps = verifier.local_caps(shapes)
        exponent = verifier.exact(ALPHA)
        entries = [arb(0) for _ in range(6)]
        for occ, weight in weights.items():
            guard()
            mass = verifier.exact(weight)
            for j, cap in enumerate(caps[occ]):
                entries[j] += mass*verifier.positive(cap)**exponent
        p = [[entries[0], entries[1]], [entries[2], entries[3]]]
        t = [[entries[4]], [entries[5]]]
        denominator = (verifier.exact(u)**f*verifier.exact(v)**s
                       *comb(4096, f)*comb(4096-f, s))
        rr = [[x/denominator for x in row] for row in verifier.matpow(p, 128)]
        tt = [[x/denominator for x in row]
              for row in verifier.matmul(verifier.matpow(p, 127), t)]
        chain = verifier.matmul(verifier.matpow(rr, 255), tt)[0][0]
        h, c, eta = source_factor(r)
        source = (verifier.exact(eta)*verifier.exact(core.BETA)**r
                  *verifier.exact(z)**(-core.LIMIT))**exponent
        number = run_total(4096, r, s)
        bound = number*source*chain
        guard()
        report = {"source_weight": r, "ordered_profile": [[f, s], [f, s]],
            "runs": s, "z": str(z), "alpha": str(ALPHA),
            "markers": {"full": str(u), "singleton": str(v)},
            "joint_missing": {"maximum_missing_groups": h, "c_h": str(c),
                              "eta_exact": str(eta)},
            "source_count_upper": "all indexed supports with this run count",
            "source_count_log2_interval": str(arb(number).log()/logtwo),
            "powered_source_log2_interval": str(source.log()/logtwo),
            "conditional_chain_log2_interval": str(chain.log()/logtwo),
            "log2_bound_interval": str(bound.log()/logtwo),
            "selected_profile_proved_below_2_pow_minus_40": bool(bound < arb(2)**(-40)),
            "seconds": time.monotonic()-started}
        reports.append(report)
        print(json.dumps({"phase": "fixed_anchor_complete", **report}), flush=True)
        del caps, weights
    print(json.dumps({"scope": "four selected actual ordered matching profiles only",
        "directed": True, "whole_code_certificate": False,
        "all_profiles_at_any_weight_covered": False, "generated_bytes": 0,
        "precision_bits": args.precision, "anchors": len(reports),
        "tests": {"exact_small_tests": check_count,
                  "each_packet_normalizer": "(1+u+v)^32 exactly",
                  "fractional_order": "local caps before route averaging",
                  "terminal_updates": 32767},
        "seconds": time.monotonic()-started, "peak_process_bytes": peak_memory(),
        "nonnegative_intersections": verifier.intersections,
        "core_pin": CORE_PIN, "adapter_pin": ADAPTER_PIN,
        "source_sha256": hashlib.sha256(Path(__file__).read_bytes()).hexdigest()}, indent=2))


if __name__ == "__main__":
    main()
