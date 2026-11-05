"""Bounded FLOAT g7 source-rank diagnostic; not a distance certificate.

Retain the frozen exact striped-pair geometry and top-32-column cutoff.
Replace only the activation/projected-rank probabilities by those for
independent nonzero GF(2^32) diagonals. Each actual rank is tagged before
averaging; no product of marginal activation and rank bounds is used.

All 16 matching profiles are covered at the specified threshold markers.
Each profile includes the full binomial source-set bound, so their maximum
is used, not a second profile union. All-quiet histories are included.
The terminal sums every discarded syndrome support; it imposes no reset
and charges only the incoming refresh, not the unused outgoing refresh.

Regional coefficients use binary64 positive arithmetic. The source factor,
cutoff weighting, and final 64-column recurrence use logs to avoid large
powers and long-path underflow. Regional underflow can still underestimate:
only a separate directed fixed-witness replay could certify a result.
No files are written. Use Python -B; the root owns all numerical runs.
"""

import argparse
from functools import lru_cache
from fractions import Fraction
import hashlib
import importlib.util
import json
import math
import os
from pathlib import Path
import shutil
import sys
import time

sys.dont_write_bytecode = True
os.environ.setdefault("OPENBLAS_NUM_THREADS", "1")
os.environ.setdefault("OMP_NUM_THREADS", "1")
import numpy as np


GROUPS, SOURCES, STATES, LENGTH = 7, 8192, 17, 1 << 21
WEIGHT_LIMIT = 7 * LENGTH // 20
FROZEN_HASH = "04b62fddce34dcdf60f14d9969573fbf0818a4755eeffb593845e37a955e9053"
FROZEN_PATH = (Path(__file__).resolve().parents[2] / "char2_field_spin" /
               "theory" / "cheap_lane_rank_debit_gate.py")
DEADLINE = None
MEMORY_MIB = 1024


def resident_bytes():
    if sys.platform == "win32":
        import ctypes
        from ctypes import wintypes

        class Counters(ctypes.Structure):
            _fields_ = [("cb", wintypes.DWORD), ("PageFaultCount", wintypes.DWORD)] + [
                (name, ctypes.c_size_t) for name in ("PeakWorkingSetSize", "WorkingSetSize",
                    "QuotaPeakPagedPoolUsage", "QuotaPagedPoolUsage", "QuotaPeakNonPagedPoolUsage",
                    "QuotaNonPagedPoolUsage", "PagefileUsage", "PeakPagefileUsage")]

        counters = Counters()
        counters.cb = ctypes.sizeof(counters)
        kernel = ctypes.WinDLL("kernel32", use_last_error=True)
        psapi = ctypes.WinDLL("psapi", use_last_error=True)
        kernel.GetCurrentProcess.restype = wintypes.HANDLE
        psapi.GetProcessMemoryInfo.argtypes = (wintypes.HANDLE, ctypes.c_void_p, wintypes.DWORD)
        if not psapi.GetProcessMemoryInfo(kernel.GetCurrentProcess(), ctypes.byref(counters), counters.cb):
            raise OSError(ctypes.get_last_error(), "GetProcessMemoryInfo failed")
        return int(counters.PeakWorkingSetSize)
    import resource
    value = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
    return int(value if sys.platform == "darwin" else 1024 * value)


def guard():
    if DEADLINE is not None and time.monotonic() >= DEADLINE:
        raise TimeoutError("bounded direct-rank diagnostic reached its time cap")
    if resident_bytes() > MEMORY_MIB * 1024**2:
        raise MemoryError("bounded direct-rank diagnostic reached its memory cap")


def load_geometry():
    # A bounded remote source-only run stages this exact dependency beside
    # the driver. Either layout must satisfy the same content hash.
    path = FROZEN_PATH if FROZEN_PATH.is_file() else Path(__file__).with_name(FROZEN_PATH.name)
    if hashlib.sha256(path.read_bytes()).hexdigest() != FROZEN_HASH:
        raise RuntimeError("frozen paired-geometry helper hash mismatch")
    spec = importlib.util.spec_from_file_location("small_state_frozen_rank_geometry", path)
    frozen = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(frozen)
    # The separate module object is used only as a geometry/arithmetic library.
    # Its source and local probability formulas remain untouched.
    frozen.guard = guard
    return frozen


@lru_cache(maxsize=2048)
def activation_rank_floor(occupancy, checks):
    """Frozen four-block coordinate-intersection knapsack, field independent."""
    value = [0] + [-100] * checks
    for ell in occupancy:
        dimension = min(4, ell)
        following = [-100] * (checks + 1)
        for used, rank in enumerate(value):
            if rank < 0:
                continue
            for coordinates in range(min(4, checks - used) + 1):
                following[used + coordinates] = max(following[used + coordinates],
                    rank + max(0, dimension - 4 + coordinates))
        value = following
    return max(value)


def log2_projective(dimension, field_bits):
    if dimension < 1:
        raise ValueError("Projective dimension must be positive")
    # log2((q^dimension-1)/(q-1)), without constructing large powers.
    return ((dimension - 1) * field_bits
            + (math.log1p(-math.exp2(-field_bits * dimension))
               - math.log1p(-math.exp2(-field_bits))) / math.log(2.))


def rank_cap(dimension, checks, rank, field_bits):
    if not 0 <= rank <= dimension or checks < 0:
        raise ValueError("Invalid rank-class arguments")
    if rank == dimension:
        return 1.
    log_q_minus_one = field_bits + math.log1p(-math.exp2(-field_bits)) / math.log(2.)
    exponent = (log2_projective(dimension, field_bits) - checks * log_q_minus_one
                - log2_projective(dimension - rank, field_bits))
    return math.exp2(min(0., exponent))


@lru_cache(maxsize=2048)
def local_terms(occupancy, tau, field_bits):
    """Tuples (rank debit, incoming weight, outgoing weight, weight floor, mass)."""
    c = sum(min(4, ell) for ell in occupancy)
    occupied = sum(ell != 0 for ell in occupancy)
    terms = [(c, 0, 0, 0, 1.)]
    for h in range(STATES):
        blocks = (h + 3) // 4
        emitted = max(0, 13 * blocks - sum(occupancy[:blocks]))
        for j in range(STATES):
            m = 16 - j
            if h == 0:
                if not (c and j >= max(1, 17 - 4 * occupied)):
                    continue
                lower = max(0, c - j, activation_rank_floor(occupancy, m))
                ranks = ((d, rank_cap(c, m, d, field_bits))
                         for d in range(lower, min(m, c - 1) + 1))
                weight, coefficient = 1, float(math.comb(16, j))
            else:
                allowed = ((j == 0 and c and h <= 4 * occupied)
                           or (j > 0 and h + j >= 17 - 4 * occupied))
                if not allowed:
                    continue
                k = max(0, c - h)
                if j == 0:
                    ranks = ((k, 1.),)
                elif h >= m or k == 0:
                    ranks = ((0, 1.),)
                else:
                    checks, slots = m - h, 16 - h
                    ranks = ((d, rank_cap(k, checks, d, field_bits))
                             for d in range(max(0, checks + k - slots), min(checks, k) + 1))
                weight = emitted
                coefficient = math.comb(16, j) * math.exp2(-min(m, h) * math.log2(tau))
            for d, cap in ranks:
                terms.append((d, h, j, weight, coefficient * cap))
    return tuple(terms)


def regional(singles, full, tau, z, field_bits, frozen):
    base = {}
    for a in range(singles + 1):
        for b in range(full + 1):
            for occupancy, count in frozen.occupancy_counts(a, b):
                for d, h, j, weight, mass in local_terms(occupancy, tau, field_bits):
                    key = (a, b, d)
                    if key not in base:
                        base[key] = np.zeros((STATES, STATES))
                    base[key][h, j] += count * z**weight * mass
    value = base
    for _ in range(7):
        guard()
        value = frozen.product(value, value, singles, full)
    normalizer = math.comb(4096, singles) * math.comb(4096 - singles, full) * 2**singles
    result = [np.zeros((STATES, STATES)) for _ in range(singles + 2 * full + 1)]
    for (a, b, d), matrix in value.items():
        if a == singles and b == full:
            result[d] += matrix / normalizer
    if any(np.any(~np.isfinite(matrix)) or np.any(matrix < 0) for matrix in result):
        raise ArithmeticError("Regional positive coefficient overflow or invalid sign")
    return result


def logarithms(value):
    with np.errstate(divide="ignore"):
        return np.log2(value)


def advance(log_state, log_matrix):
    return np.logaddexp2.reduce(log_state[..., :, None] + log_matrix, axis=-2)


def threshold_bound(coefficients, tau, z, threshold, log2_x):
    guard()
    limit = 128 * GROUPS - 32 * threshold
    log_tau = math.log2(tau)
    if threshold == 4 * GROUPS - 1:
        log_zero = np.logaddexp2.reduce(np.asarray([logarithms(m) for m in coefficients[:-1]]))
        log_one = logarithms(coefficients[-1]) - log_tau
        state = np.full((32, STATES), -math.inf)
        state[0, 0] = 0.
        for _ in range(64):
            following = advance(state, log_zero)
            following[1:] = np.logaddexp2(following[1:], advance(state[:-1], log_one))
            state = following
        factors = np.asarray([math.log1p(-math.exp2(-(32-e)*log_tau))/math.log(2.)
                              for e in range(32)])
        moment = float(np.logaddexp2.reduce((state + factors[:, None]).ravel()))
    else:
        terms = np.asarray([logarithms(matrix)
                            + (log2_x-log_tau)*max(d-threshold, 0)
                            for d, matrix in enumerate(coefficients)])
        matrix = np.logaddexp2.reduce(terms)
        state = np.full(STATES, -math.inf)
        state[0] = 0.
        for _ in range(64):
            state = advance(state, matrix)
        moment = float(np.logaddexp2.reduce(state)) - (limit-1)*log2_x
    return (math.log2(math.comb(SOURCES, GROUPS)) + limit*log_tau
            - math.log2(tau-1) - WEIGHT_LIMIT*math.log2(z) + moment)


def selected_markers(central_z=None):
    # Reuse the retained finite threshold schedule, not a marker grid.
    for v in range(4 * GROUPS - 1):
        if v == 2 * GROUPS and central_z is not None:
            yield v, central_z, 0.
        else:
            yield v, 1 - Fraction(27, 2000000) * (128*GROUPS - 32*v), -.5
    for mark in (".9993", ".9995", ".9997"):
        yield 4 * GROUPS - 1, Fraction(mark), None


def screen(field_bits=32, tau=65280., seconds=50., memory_mib=1024, profiles=None,
           central_z=None):
    if (field_bits not in (8, 16, 32) or not 1 < tau <= (1 << field_bits)-1
            or not 0 < seconds <= 50 or not 64 <= memory_mib <= 1024):
        raise ValueError("Require field bits8/16/32, 1<tau<=q-1, <=50 seconds and <=1024 MiB")
    central_z = None if central_z is None else Fraction(central_z)
    if central_z is not None and not 0 < central_z < 1:
        raise ValueError("Require a central output marker strictly between zero and one")
    all_profiles = tuple((a, b) for a in range(4) for b in range(4))
    profiles = all_profiles if profiles is None else tuple(profiles)
    if not profiles or len(set(profiles)) != len(profiles) or any(p not in all_profiles for p in profiles):
        raise ValueError("Invalid g7 paired profiles")
    free = shutil.disk_usage(Path(__file__).resolve().parent).free
    if free < 20 * 1024**3:
        raise RuntimeError("At least20 GiB free space is required")
    global DEADLINE, MEMORY_MIB
    started = time.monotonic()
    DEADLINE, MEMORY_MIB = started + seconds, memory_mib
    frozen = load_geometry()
    best = {p: {} for p in profiles}
    used = sorted({f for p in profiles for f in p})
    evaluated = 0
    status = "complete-selected-markers"
    try:
        for v, exact_z, log2_x in selected_markers(central_z):
            guard()
            z = float(exact_z)
            regions = {f: regional(GROUPS-2*f, f, tau, z, field_bits, frozen) for f in used}
            for p in profiles:
                guard()
                half = frozen.column_product(regions[p[0]], regions[p[1]])
                coefficients = frozen.column_product(half, half)
                value = threshold_bound(coefficients, tau, z, v, log2_x)
                if not math.isfinite(value):
                    raise ArithmeticError("Nonfinite floating threshold score")
                evaluated += 1
                if v not in best[p] or value < best[p][v]["log2_bound"]:
                    best[p][v] = dict(v=v, z=str(exact_z), log2_x=log2_x,
                                      exact_32_counter=v == 27, log2_bound=value)
    except (TimeoutError, MemoryError) as error:
        status = str(error)
    rows = []
    for p in profiles:
        terms = [best[p][v] for v in sorted(best[p])]
        complete = len(terms) == 28
        rows.append(dict(even_full=p[0], odd_full=p[1], complete_threshold_cover=complete,
            checked_thresholds=len(terms),
            log2_union_bound=(float(np.logaddexp2.reduce([t["log2_bound"] for t in terms]))
                              if complete else None),
            worst_threshold=max(terms, key=lambda t: t["log2_bound"]) if terms else None,
            thresholds=terms))
    complete = all(row["complete_threshold_cover"] for row in rows)
    return dict(schema="small-state-direct-rank-float-v1", status=status,
        groups=GROUPS, coefficient_field_bits=field_bits, tau=tau, distance="7/20",
        source_dimension=128*GROUPS, top_columns=32, threshold_count=28,
        central_threshold_override=(dict(v=2*GROUPS, z=str(central_z), log2_x=0.)
                                    if central_z is not None else None),
        evaluated_profile_markers=evaluated, profiles=rows,
        complete_selected_profile_cover=complete,
        all_16_profiles=len(profiles) == 16,
        max_log2_union_bound=max(row["log2_union_bound"] for row in rows) if complete else None,
        all_quiet_histories="included", terminal="all final syndrome supports; no final reset",
        geometry="exact striped pair slots and fair orientations",
        arithmetic="binary64 positive coefficients; log2 cutoff and final recurrence; not directed",
        frozen_dependency=str(FROZEN_PATH.name), frozen_sha256=FROZEN_HASH,
        seconds=time.monotonic()-started, peak_resident_bytes=resident_bytes(), free_bytes_at_start=free,
        directed=False, whole_code_certificate=False, generated_file_bytes=0)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--field-bits", type=int, choices=(8, 16, 32), default=32)
    parser.add_argument("--tau", type=float, default=65280.)
    parser.add_argument("--seconds", type=float, default=50.)
    parser.add_argument("--memory-mib", type=int, default=1024)
    parser.add_argument("--profile", type=int, nargs=2, help="Restrict to one diagnostic profile")
    parser.add_argument("--central-z", type=Fraction,
                        help="Override only threshold14 with this rational z and x=1")
    parser.add_argument("--details", action="store_true")
    args = parser.parse_args()
    profiles = [tuple(args.profile)] if args.profile is not None else None
    result = screen(args.field_bits, args.tau, args.seconds, args.memory_mib, profiles,
                    args.central_z)
    if not args.details:
        for row in result["profiles"]:
            del row["thresholds"]
    print(json.dumps(result, allow_nan=False, sort_keys=True), flush=True)


if __name__ == "__main__":
    main()
