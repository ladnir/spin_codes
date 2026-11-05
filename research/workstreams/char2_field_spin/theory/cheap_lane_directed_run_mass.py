"""Directed run-mass join for normalized r1024 regional upper bounds.

Input regions[f]=(R_f,V_f) contain four row-major regular entries and two
terminal entries. Each entry already includes its own route coefficient
normalizer. Inputs are nonnegative Arb enclosures at one common z/alpha.
The helper does not apply a source factor, a new alpha, or marker weights.

For an ordered pair(a,b), the conditional chain is
 e0 (R_a R_b)^127 R_a V_b.
For each cyclic run count c, maximize this chain over all degree-feasible
ordered pairs a+b=r-c. Some pairs can have zero actual count; including
them only increases the bound. Multiply the maximum by the exact number
N_c of indexed cyclic supports, then sum over c. There are no BEST queries.

One power supplies both orders of each unordered pair. Positive explicit
2x2 products and Arb.max preserve enclosures. The final result bounds the
source-counted conditional sum only. The caller multiplies its proved
common source factor exactly once and owns the overall job's limits.

Standalone execution runs small exact tests, not the r1024 matrix join.
It has a 180-second/256MiB/20GiB guard, writes only stdout, and creates no
coefficient files. No numerical computation occurs on import.
"""

from collections import Counter
from fractions import Fraction as F
import hashlib
import json
from math import comb
from pathlib import Path
import shutil
import sys
import time

sys.dont_write_bytecode = True
WEIGHT, POINTS, POWER_BITS = 1024, 8192, 7


def product(a, b):
    """Positive row-major 2x2 product, with no generic array machinery."""
    return (a[0]*b[0]+a[1]*b[2], a[0]*b[1]+a[1]*b[3],
            a[2]*b[0]+a[3]*b[2], a[2]*b[1]+a[3]*b[3])


def square(a):
    trace, cross = a[0]+a[3], a[1]*a[2]
    return (a[0]*a[0]+cross, a[1]*trace,
            a[2]*trace, a[3]*a[3]+cross)


def power_ones(a, bits=POWER_BITS):
    """A^(2^bits-1), using six squarings at the production exponent127."""
    if not 1 <= bits <= POWER_BITS:
        raise ValueError("power bits must be1..7")
    result = a
    for _ in range(bits-1):
        result = product(square(result), a)
    return result


def two_orientations(region_a, region_b, bits=POWER_BITS, diagonal=False):
    """Share (A B)^(2^bits-1), retaining each order's own terminal."""
    a, va = region_a
    b, vb = region_b
    power = power_ones(product(a, b), bits)
    x0, x1 = a[0]*vb[0]+a[1]*vb[1], a[2]*vb[0]+a[3]*vb[1]
    forward = power[0]*x0+power[1]*x1
    if diagonal:
        return forward, forward
    x0, x1 = power[0]*va[0]+power[1]*va[1], power[2]*va[0]+power[3]*va[1]
    reverse = b[0]*x0+b[1]*x1
    return forward, reverse


def run_counts(points=POINTS, weight=WEIGHT):
    """Exact indexed cyclic support counts; rotations are not quotiented."""
    if points < 2 or points % 2 or not 0 < weight <= points//2:
        raise ValueError("require even points and1<=weight<=points/2")
    result = [0]*(weight+1)
    value = points
    for c in range(1, weight+1):
        result[c] = value
        if c < weight:
            numerator = value*(weight-c)*(points-weight-c)
            denominator = c*(c+1)
            value, remainder = divmod(numerator, denominator)
            if remainder:
                raise ArithmeticError("nonintegral run-count recurrence")
    if sum(result) != comb(points, weight):
        raise ArithmeticError("run counts do not cover all indexed supports")
    return result


def directed_join(regions, guard, progress=None, weight=WEIGHT, points=POINTS):
    """Return {'bound': Arb, 'run_bounds': tuple, 'report': dict}.

    regions must have floor(weight/2)+1 entries. R and V are flattened
    tuples, not scaled FLOAT tuples. Their denominators must already be
    removed. The caller supplies guard(); progress is optional.
    """
    from flint import arb
    started = time.monotonic()
    counts = run_counts(points, weight)
    full = weight//2
    if len(regions) != full+1:
        raise ValueError("missing regional marginal")
    clean, intersections = [], 0
    for region in regions:
        guard()
        if len(region) != 2 or len(region[0]) != 4 or len(region[1]) != 2:
            raise ValueError("expected flat Arb R4/V2 with no scales or denominator")
        rr = []
        for row in region:
            values = []
            for value in row:
                if not isinstance(value, arb) or not value.is_finite() or value.upper() < 0:
                    raise ArithmeticError("regional entry must be a finite nonnegative Arb enclosure")
                if not value >= 0:
                    intersections += 1
                values.append(value.nonnegative_part())
            rr.append(tuple(values))
        clean.append(tuple(rr))
    total = arb(0)
    per_run = []
    ordered = unordered = 0
    for c in range(1, weight+1):
        guard()
        degree_sum = weight-c
        lo, hi = max(0, degree_sum-full), min(full, degree_sum)
        maximum = arb(0)
        for a in range(lo, min(hi, degree_sum//2)+1):
            b = degree_sum-a
            if unordered % 128 == 0:
                guard()
            left, right = two_orientations(clean[a], clean[b], diagonal=a == b)
            maximum = maximum.max(left)
            if a != b:
                maximum = maximum.max(right)
            ordered += 1 if a == b else 2
            unordered += 1
        if hi < lo:
            raise ArithmeticError("empty run-class profile range")
        term = counts[c]*maximum
        if not term.is_finite() or term.upper() < 0:
            raise ArithmeticError("invalid positive run contribution")
        term = term.nonnegative_part()
        total += term
        per_run.append(term)
        if progress is not None and (c == 1 or c % 64 == 0 or c == weight):
            progress("directed_run_mass", run=c, ordered_profiles=ordered,
                     unordered_powers=unordered, seconds=time.monotonic()-started)
    expected = (full+1)**2-(1 if weight % 2 == 0 else 0)
    if ordered != expected:
        raise ArithmeticError("ordered-profile coverage mismatch")
    if weight == WEIGHT and points == POINTS and ordered != 263168:
        raise ArithmeticError("r1024 profile control mismatch")
    guard()
    logtwo = arb(2).log()
    selected = sorted({1, weight, *(c for c in (2, 8, 16, 32, 64, 128, 256, 512)
                                     if c <= weight)})
    report = {
        "scope": "exact run mass times directed maximum ordered conditional chain",
        "weight": weight, "points": points, "source_factor_applied": False,
        "new_alpha_applied": False, "regional_normalizers_already_removed": True,
        "ordered_profiles": ordered, "unordered_powers": unordered,
        "all_run_counts_sum_to_binomial": True, "exact_count_queries": 0,
        "regional_order": "e0 (Ra Rb)^127 Ra Vb; reverse retains Va",
        "regional_factors": 256, "terminal_region_updates": 127,
        "packet_updates": 32767, "input_nonnegative_intersections": intersections,
        "conditional_source_counted_log2_interval": str(total.log()/logtwo) if total > 0 else "zero",
        "selected_run_log2_intervals": {str(c): str(per_run[c-1].log()/logtwo)
            if per_run[c-1] > 0 else "zero" for c in selected},
        "seconds": time.monotonic()-started,
    }
    return {"bound": total, "run_bounds": tuple(per_run), "report": report}


def exact_tests(guard=lambda: None):
    """Independent Fraction chains and exhaustive small indexed cycles."""
    checks = masks = 0
    a = tuple(map(F, (2, 1, 3, 2)))
    b = (F(3, 2), F(2, 3), F(1, 5), F(7, 6))
    va, vb = (F(2, 7), F(3, 11)), (F(5, 13), F(7, 17))
    def direct(left, right, terminal, repetitions):
        vector = terminal
        for step in range(2*repetitions+1):
            matrix = left if step % 2 == 0 else right
            vector = (matrix[0]*vector[0]+matrix[1]*vector[1],
                      matrix[2]*vector[0]+matrix[3]*vector[1])
        return vector[0]
    for bits in (1, 2, 3, 7):
        guard()
        forward, reverse = two_orientations((a, va), (b, vb), bits)
        assert forward == direct(a, b, vb, (1 << bits)-1)
        assert reverse == direct(b, a, va, (1 << bits)-1)
        assert square(a) == product(a, a)
        checks += 3
    for points in (6, 8, 10):
        for weight in range(1, points//2+1):
            guard()
            counted, profiles = Counter(), Counter()
            for mask in range(1 << points):
                if mask.bit_count() != weight:
                    continue
                bits = [(mask >> i) & 1 for i in range(points)]
                c = sum(bits[i] and not bits[(i-1) % points] for i in range(points))
                f0 = sum(bits[2*i]*bits[2*i+1] for i in range(points//2))
                f1 = sum(bits[2*i+1]*bits[(2*i+2) % points] for i in range(points//2))
                assert c == weight-f0-f1
                counted[c] += 1; profiles[f0, f1] += 1; masks += 1
            counts = run_counts(points, weight)
            assert [counted[c] for c in range(weight+1)] == counts
            regions = [((F(1+f, 9), F(2, 11), F(3, 13), F(4+f, 17)),
                        (F(5+f, 19), F(7, 23+f))) for f in range(weight//2+1)]
            values = {}
            exact = F(0)
            for (f0, f1), number in profiles.items():
                value = direct(regions[f0][0], regions[f1][0], regions[f1][1], 3)
                values[f0, f1] = value
                exact += number*value
            envelope = F(0)
            for c in range(1, weight+1):
                candidates = []
                for f0 in range(weight//2+1):
                    f1 = weight-c-f0
                    if 0 <= f1 <= weight//2:
                        # Includes degree-feasible pairs of actual multiplicity0.
                        candidates.append(direct(regions[f0][0], regions[f1][0], regions[f1][1], 3))
                envelope += counts[c]*max(candidates)
            assert exact <= envelope
            checks += 2
    counts = run_counts()
    assert len(counts) == 1025 and counts[1] == POINTS and sum(counts) == comb(POINTS, WEIGHT)
    return {"fraction_chain_and_cover_checks": checks, "exhaustive_masks": masks,
            "actual1024_run_count_checksum": hashlib.sha256(str(sum(counts)).encode()).hexdigest(),
            "all_tests_exact": True}


def main():
    if shutil.disk_usage(Path(__file__).resolve().anchor).free < 20*(1 << 30):
        raise RuntimeError("requires20GiB free on source volume")
    started = time.monotonic()
    if sys.platform == "win32":
        import ctypes
        from ctypes import wintypes
        class Counters(ctypes.Structure):
            _fields_ = [("cb", wintypes.DWORD), ("PageFaultCount", wintypes.DWORD),
                        *[(name, ctypes.c_size_t) for name in ("PeakWorkingSetSize", "WorkingSetSize",
                          "QuotaPeakPagedPoolUsage", "QuotaPagedPoolUsage", "QuotaPeakNonPagedPoolUsage",
                          "QuotaNonPagedPoolUsage", "PagefileUsage", "PeakPagefileUsage")]]
        process = ctypes.windll.kernel32.GetCurrentProcess
        process.restype = wintypes.HANDLE
        info = ctypes.windll.psapi.GetProcessMemoryInfo
        info.argtypes = [wintypes.HANDLE, ctypes.POINTER(Counters), wintypes.DWORD]
        def peak():
            result = Counters(); result.cb = ctypes.sizeof(result)
            if not info(process(), ctypes.byref(result), result.cb):
                raise OSError("cannot inspect memory budget")
            return result.PeakWorkingSetSize
    else:
        import resource
        def peak():
            value = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
            return value if sys.platform == "darwin" else 1024*value
    def guard():
        if time.monotonic()-started > 180:
            raise TimeoutError("180-second exact-test deadline")
        if peak() > 256*(1 << 20):
            raise MemoryError("256MiB exact-test budget exceeded")
    tests = exact_tests(guard)
    guard()
    print(json.dumps({"scope": "standalone exact tests only; no r1024 matrix join",
        "tests": tests, "seconds": time.monotonic()-started, "peak_process_bytes": peak(),
        "generated_bytes": 0, "source_sha256": hashlib.sha256(Path(__file__).read_bytes()).hexdigest()}), flush=True)


if __name__ == "__main__":
    main()
