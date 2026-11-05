"""Directed fixed-witness replay for one small-state source count7..32.

Default to g7,z=9634151/9689600. For larger g, interpolate rationally to
z16=2238127/2263291 through g16; above g16 an explicit rational --z
is required. An explicit marker may also replace the interpolation.
Default to alpha=1; an explicit rational --alpha may lie in(0,1].
Fix missing-group cutoff0 and c=(256-g)/256. Apply alpha entrywise before
route averaging and to the conditional source factor, not the source count.
Compute the preserved local rational majorants without float conversion.
Arb then encloses every positive coefficient operation, the exact regional
normalization, and the complete source-prefactored bound for every allowed
profile0<=f0,f1<=floor(g/2), f0+f1<g.

Acceptance compares each complete bound directly with the exact value
2^-target_bits. Displayed log intervals are diagnostic, not the test.
The binomial factor already counts all source sets; use the maximum over
profiles, not a second profile union. This is a fixed-source-count class bound,
not a whole-code certificate. No output files. Run with Python -B.
"""

import argparse
from fractions import Fraction
from functools import lru_cache
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


GROUPS, LENGTH, Q0 = 7, 1 << 21, 256
Z = Fraction(9634151, 9689600)
Z16 = Fraction(2238127, 2263291)
WEIGHT_LIMIT = 7 * LENGTH // 20
FROZEN_HASH = "04b62fddce34dcdf60f14d9969573fbf0818a4755eeffb593845e37a955e9053"
OLD_THEORY = Path(__file__).resolve().parents[2] / "char2_field_spin" / "theory"
DEPENDENCIES = (
    "hd_common_refresh_packets", "hd_common_refresh_principal", "hd_bucket_feedback",
    "hd_systematic_mds", "hd_blockwise_shared_transfer", "hd_sumrank_maps",
    "hd_blockwise_maps", "hd_blockwise_inner",
)


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


def verify(precision=128, target_bits=900, seconds=50., memory_mib=1024,
           groups=GROUPS, z=None, alpha=Fraction(1)):
    if (precision not in (128, 256) or not 1 <= target_bits <= 10000
            or not isinstance(groups, int) or not 7 <= groups <= 32
            or not 0 < seconds <= 50 or not 64 <= memory_mib <= 1024):
        raise ValueError("Require groups7..32, precision128/256, bounded target, <=50 seconds and <=1024 MiB")
    if groups > 16 and z is None:
        raise ValueError("Require an explicit rational output marker above g16")
    explicit_z = z is not None
    z = Z + (Z16-Z)*Fraction(groups-7, 9) if z is None else Fraction(z)
    if not 0 < z < 1:
        raise ValueError("Require a rational output marker strictly between zero and one")
    alpha = Fraction(alpha)
    if not 0 < alpha <= 1:
        raise ValueError("Require a rational fractional marker in(0,1]")
    free = shutil.disk_usage(Path(__file__).resolve().parent).free
    if free < 20 * 1024**3:
        raise RuntimeError("At least20 GiB free space is required")
    started = time.monotonic()

    def guard():
        if time.monotonic()-started >= seconds:
            raise TimeoutError("directed fixed-witness replay reached its time cap")
        if resident_bytes() > memory_mib * 1024**2:
            raise MemoryError("directed fixed-witness replay reached its memory cap")

    directory = OLD_THEORY if OLD_THEORY.is_dir() else Path(__file__).resolve().parent
    sys.path.insert(0, str(directory))
    geometry_path = directory / "cheap_lane_rank_debit_gate.py"
    if hashlib.sha256(geometry_path.read_bytes()).hexdigest() != FROZEN_HASH:
        raise RuntimeError("frozen striped-geometry hash mismatch")
    spec = importlib.util.spec_from_file_location("small_state_arb_route_geometry", geometry_path)
    geometry = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(geometry)
    from hd_common_refresh_packets import local_moments, matrix
    from flint import arb, arb_mat, ctx
    ctx.prec = precision

    def rational(value):
        value = Fraction(value)
        return arb(value.numerator) / value.denominator

    def zero():
        return arb_mat(2, 2)

    def power(value, exponent):
        result = arb_mat([[1, 0], [0, 1]])
        while exponent:
            guard()
            if exponent & 1:
                result = result * value
            exponent >>= 1
            if exponent:
                value = value * value
        return result

    local = local_moments(z, q=Q0, h=4)
    epsilon = Fraction(1, Q0**16)
    guard()

    @lru_cache(maxsize=1024)
    def packet(occupancy):
        # All min operations in the provider compare exact Fractions.
        quiet, birth, ret, live = matrix(occupancy, z, local, q=Q0, h=4)
        total = quiet + birth
        all_live = epsilon*total + (1-epsilon)*live
        entries = ((quiet, total, epsilon*quiet+(1-epsilon)*ret, all_live)
                   if sum(occupancy) else (Fraction(1), Fraction(0), Fraction(0), all_live))
        if min(entries) < 0 or max(entries) > 1:
            raise ArithmeticError("Exact rational local majorant left [0,1]")
        result = zero()
        for index, value in enumerate(entries):
            entry = rational(value)
            if alpha != 1 and value:
                entry = entry**rational(alpha)
            result[index//2, index%2] = entry
        return result

    def product(left, right, singles, full):
        guard()
        result = {}
        for (a, b), aa in left.items():
            for (c, d), bb in right.items():
                if a+c <= singles and b+d <= full:
                    key = (a+c, b+d)
                    value = aa*bb
                    if key in result:
                        result[key] += value
                    else:
                        result[key] = value
        return result

    normalization_checks = []

    def regional(full):
        singles = groups-2*full
        base, counts = {}, {}
        for a in range(singles+1):
            for b in range(full+1):
                value, mass = zero(), 0
                for occupancy, count in geometry.occupancy_counts(a, b):
                    # The local formula is symmetric in the four blocks.
                    value += packet(tuple(sorted(occupancy))) * count
                    mass += count
                base[(a, b)], counts[(a, b)] = value, mass
        value = base
        for _ in range(7):
            value = product(value, value, singles, full)
            following = {}
            for (a, b), left in counts.items():
                for (c, d), right in counts.items():
                    if a+c <= singles and b+d <= full:
                        key = (a+c, b+d)
                        following[key] = following.get(key, 0) + left*right
            counts = following
        normalizer = (math.comb(4096, singles) * math.comb(4096-singles, full)
                      * 2**singles)
        if counts[(singles, full)] != normalizer:
            raise ArithmeticError("Exact regional route normalization failed")
        normalization_checks.append(dict(singles=singles, full=full, normalizer=str(normalizer)))
        return value[(singles, full)] / arb(normalizer)

    # For q=2^32, Gamma <=(1+1/(q-1))^(16*n)<2 follows by
    # bounding binomial coefficients by powers and summing a geometric series.
    refresh_exponent = 16 * (LENGTH//64)
    if not 2*refresh_exponent < (1 << 32)-1:
        raise ArithmeticError("The uniform Gamma<2 receipt is not applicable")
    retained = Fraction(1)-Fraction(groups, Q0)
    if retained != Fraction(Q0-groups, Q0) or retained <= 0:
        raise ArithmeticError("Incorrect fixed missing-group cutoff")
    beta = arb(Q0)**160 / arb(Q0-1)**32
    if alpha == 1:
        # Preserve the prior fixed-witness arithmetic exactly.
        prefactor = (math.comb(8192, groups) * 2 * beta**groups
                     / ((Q0-1)*rational(retained)) * rational(z)**(-WEIGHT_LIMIT))
    else:
        conditional = (2 * beta**groups / ((Q0-1)*rational(retained))
                       * rational(z)**(-WEIGHT_LIMIT))
        prefactor = math.comb(8192, groups) * conditional**rational(alpha)
    regions = {full: regional(full) for full in range(groups//2+1)}
    target = arb(2)**(-target_bits)
    reports = []
    for even in range(groups//2+1):
        for odd in range(groups//2+1):
            # A proper subset of the8192-cycle has at most g-1 internal
            # matching edges across the two alternating pairings.
            if even+odd >= groups:
                continue
            guard()
            run = power(regions[even]*regions[odd], 128)
            bound = prefactor * (run[0, 0]+run[0, 1])
            finite = bool(bound.is_finite())
            passed = bool(finite and bound >= 0 and bound < target)
            reports.append(dict(even_full=even, odd_full=odd, passed=passed,
                log2_bound_interval=str(bound.log()/arb(2).log()) if finite else "nonfinite"))
    guard()
    dependencies = {name + ".py": hashlib.sha256(Path(sys.modules[name].__file__).read_bytes()).hexdigest()
                    for name in DEPENDENCIES}
    dependencies[geometry_path.name] = FROZEN_HASH
    return dict(schema=f"small-state-exact-route-arb-g{groups}-v1", status="directed-fixed-witness-replay",
        passed=all(row["passed"] for row in reports), precision_bits=precision,
        target_bits=target_bits, acceptance="every full profile bound < exact2^(-target_bits)",
        groups=groups, coefficient_field_bits=32, auxiliary_field_order=Q0,
        z=str(z), marker_rule="explicit override" if explicit_z else "rational interpolation z7..z16",
        alpha=str(alpha), source_mode="cutoff", missing_groups=0,
        retained_fraction=str(retained), gamma_upper="2",
        source_prefactor=(f"binom(8192,{groups})*2*(256^160/255^32)^{groups}"
                          f"/[255*({retained})]*z^(-{WEIGHT_LIMIT})" if alpha == 1 else
                          f"binom(8192,{groups})*(2*(256^160/255^32)^{groups}"
                          f"/[255*({retained})]*z^(-{WEIGHT_LIMIT}))^({alpha})"),
        output_limit=WEIGHT_LIMIT, output_length=LENGTH, profiles=reports,
        exact_integer_normalization_checks=normalization_checks,
        terminal="positive f<=M*ones; no physical final refresh or reset charge",
        profile_accounting="full source-set factor in each profile; maximum, no extra profile sum",
        scope=f"all nonzero messages on exactly{groups} of8192 groups; ideal setup",
        seconds=time.monotonic()-started, peak_resident_bytes=resident_bytes(),
        free_bytes_at_start=free, dependencies=dependencies,
        directed=True, whole_code_certificate=False, generated_file_bytes=0)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--precision", type=int, choices=(128, 256), default=128)
    parser.add_argument("--target-bits", type=int, default=900)
    parser.add_argument("--seconds", type=float, default=50.)
    parser.add_argument("--memory-mib", type=int, default=1024)
    parser.add_argument("--groups", type=int, choices=range(7, 33), default=GROUPS)
    parser.add_argument("--z", type=Fraction, help="Override the fixed rational interpolation marker")
    parser.add_argument("--alpha", type=Fraction, default=Fraction(1),
                        help="Rational entrywise fractional marker in(0,1]")
    args = parser.parse_args()
    result = verify(args.precision, args.target_bits, args.seconds, args.memory_mib,
                    args.groups, args.z, args.alpha)
    print(json.dumps(result, allow_nan=False, sort_keys=True), flush=True)
    if not result["passed"]:
        raise SystemExit(1)


if __name__ == "__main__":
    main()
