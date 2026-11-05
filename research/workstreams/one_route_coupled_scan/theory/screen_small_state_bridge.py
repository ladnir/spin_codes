"""Bounded FLOAT diagnostic for the whole-state GF32 bridge, not a certificate.

Each printed row concerns one actual source count and one fixed pair profile.
The binomial source factor overcounts its possible indexed source sets.
No interpolation between these profiles, weights, or markers is claimed.

The local rational majorants come from the preserved GF256 RS[16,4] maps.
Only their matrix arithmetic is converted to binary64. The regional bound
conditions independent pair types on their actual without-replacement counts.
Its coefficient bound is deliberately coarse and requires no FFT or search.

The final vector is bounded by M*ones, including after entrywise alpha.
This permits n matrices in the diagnostic; it is not a probability charge
for the unused final refresh. The uniform refresh-density factor Gamma<2
is included once in the source prefactor.

No output files or caches are written. Run with Python -B. The root owns runs.
"""

import argparse
from collections import Counter
from fractions import Fraction
from functools import lru_cache
from itertools import combinations_with_replacement, product
import json
import math
from pathlib import Path
import shutil
import sys
import time

import numpy as np


OLD_THEORY = Path(__file__).resolve().parents[2] / "char2_field_spin" / "theory"
sys.path.insert(0, str(OLD_THEORY))
from hd_common_refresh_packets import local_moments, matrix  # noqa: E402


GROUPS = 8192
PAIRS = 4096
REGIONS = 256
BATCHES_PER_REGION = 128
LIMIT = (1 << 21) * 7 // 20
Q0 = 256
DEFAULT_GROUPS = (7, 16, 64, 256, 1024, 4096, 8192)
# Sparse sources need marks near one; full support needs a mark near
# delta/((q0-1)*(1-delta)) ~= .0021. The quick grid is only a smoke test.
DEFAULT_MARKERS = (".99999", ".9999", ".999", ".99", ".9", ".5",
                   ".1", ".03", ".01", ".003", ".002")
DEFAULT_ALPHAS = (.002, .01, .03, .1, .25, .5, 1.)


def log2_sum(values):
    values = np.asarray(tuple(values), dtype=float)
    if not len(values):
        return -math.inf
    return float(np.logaddexp2.reduce(values))


def log2_binomial(n, k):
    if not 0 <= k <= n:
        return -math.inf
    return (math.lgamma(n + 1) - math.lgamma(k + 1)
            - math.lgamma(n - k + 1)) / math.log(2.)


def log_matrix_product(left, right):
    """Ordered 2x2 matrix multiplication in base-two log arithmetic."""
    left, right = np.asarray(left, dtype=float), np.asarray(right, dtype=float)
    if left.shape != (2, 2) or right.shape != (2, 2):
        raise ValueError("Expected two 2x2 log matrices")
    return np.logaddexp2(left[:, 0, None] + right[None, 0, :],
                        left[:, 1, None] + right[None, 1, :])


def log_matrix_power(value, exponent):
    if not isinstance(exponent, int) or exponent < 0:
        raise ValueError("Expected a nonnegative integer power")
    value = np.asarray(value, dtype=float)
    if value.shape != (2, 2):
        raise ValueError("Expected a 2x2 log matrix")
    result = np.array(((0., -math.inf), (-math.inf, 0.)))
    while exponent:
        if exponent & 1:
            result = log_matrix_product(result, value)
        exponent >>= 1
        if exponent:
            value = log_matrix_product(value, value)
    return result


@lru_cache(maxsize=1)
def occupancy_profiles():
    """Sorted four-block occupancies and their block-order multiplicities."""
    rows = tuple(combinations_with_replacement(range(17), 4))
    multiplicities = tuple(24 // math.prod(math.factorial(v)
                                         for v in Counter(row).values())
                           for row in rows)
    return np.asarray(rows, dtype=np.int16), np.asarray(multiplicities, dtype=float)


def occupancy_distribution(probabilities, pairs=16):
    """Striped half-packet: (p0+p1*(x+y)/2+p2*x*y)^pairs.

    A physical pair feeds two DIFFERENT RS16 blocks. The two coordinates
    of a singleton are equiprobable; their occupancies are correlated.
    This is hd_shifted_alias_certificate.odd_coefficients' actual layout.
    """
    probabilities = np.asarray(probabilities, dtype=float)
    if (probabilities.shape != (3,) or np.any(probabilities < 0.)
            or not abs(float(probabilities.sum()) - 1.) < 1e-12
            or not isinstance(pairs, int) or pairs < 0):
        raise ValueError("Expected pair-type probabilities and a nonnegative pair count")
    result = np.ones((1, 1))
    for _ in range(pairs):
        extended = np.zeros((len(result)+1, len(result)+1))
        extended[:-1, :-1] += probabilities[0]*result
        extended[1:, :-1] += (probabilities[1]/2)*result
        extended[:-1, 1:] += (probabilities[1]/2)*result
        extended[1:, 1:] += probabilities[2]*result
        result = extended
    return result


@lru_cache(maxsize=1)
def sorted_profile_indices():
    profiles, _ = occupancy_profiles()
    powers = np.array((17**3, 17**2, 17, 1), dtype=np.int64)
    ordered = np.indices((17,)*4, dtype=np.int16).reshape(4, -1).T
    keys = np.sort(ordered, axis=1) @ powers
    indices = np.searchsorted(profiles @ powers, keys)
    assert np.array_equal((profiles @ powers)[indices], keys)
    return indices


def profile_weights(probabilities):
    half = occupancy_distribution(probabilities)
    ordered_weights = (half[:, :, None, None]*half[None, None, :, :]).ravel()
    return np.bincount(sorted_profile_indices(), weights=ordered_weights,
                       minlength=len(occupancy_profiles()[0]))


@lru_cache(maxsize=1)
def local_table(z):
    """Exact rational local caps converted to 4845 floating log matrices.

    The old matrix returns (k,T-k,J_nonzero/(q0^16-1),V_nonzero).
    Retain its proved return majorant and add the relaxed zero-state mass.
    Do not use the old beta-multiplied or inverse-quiet-minimized table.
    """
    z = Fraction(z)
    if not 0 < z < 1:
        raise ValueError("Expected a rational output marker strictly between zero and one")
    profiles, _ = occupancy_profiles()
    local = local_moments(z, q=Q0, h=4)
    epsilon = Fraction(1, Q0 ** 16)
    logs = np.full((len(profiles), 2, 2), -math.inf)
    nonzero_live = np.full((len(profiles), 2), -math.inf)
    zero_live = np.full((len(profiles), 2), -math.inf)
    for index, row in enumerate(profiles):
        occupancy = tuple(map(int, row))
        quiet, birth, ret, live = matrix(occupancy, z, local, q=Q0, h=4)
        total = quiet + birth
        all_live = epsilon * total + (1 - epsilon) * live
        if not sum(occupancy):
            entries = (Fraction(1), Fraction(0), Fraction(0), all_live)
        else:
            entries = (quiet, total, epsilon * quiet + (1 - epsilon) * ret, all_live)
        if min(entries) < 0 or max(entries) > 1:
            raise ArithmeticError("Local rational moment left [0,1]")
        logs[index] = np.asarray(tuple(math.log2(float(v)) if v else -math.inf
                                      for v in entries)).reshape((2, 2))
        # Retain exact positive components. Subtracting floating majorants
        # would lose or mis-sign the small nonzero-state contribution.
        for target, values in ((nonzero_live, ((1-epsilon)*ret, (1-epsilon)*live)),
                               (zero_live, (epsilon*quiet, epsilon*total))):
            target[index] = tuple(math.log2(float(v)) if v else -math.inf for v in values)
    return dict(base=logs, live_nonzero=nonzero_live, live_zero=zero_live)


def penalized_table(table, log2_penalty, retained_fraction):
    if log2_penalty < 0 or not 0 < retained_fraction <= 1:
        raise ValueError("Invalid zero-state penalty or retained source fraction")
    if not log2_penalty:
        return table["base"]
    value = table["base"].copy()
    value[:, 1, :] = (np.logaddexp2(table["live_nonzero"], table["live_zero"]-log2_penalty)
                     + log2_penalty/(256.*retained_fraction))
    value[0, 1, 0] = -math.inf  # Actual-empty batches still forbid a reset.
    return value


def pair_histogram(group_count, full_pairs):
    if (not isinstance(group_count, int) or not isinstance(full_pairs, int)
            or not 1 <= group_count <= GROUPS
            or not max(0, group_count - PAIRS) <= full_pairs <= group_count // 2):
        raise ValueError("Invalid actual group count or full-pair count")
    return (PAIRS - group_count + full_pairs, group_count - 2 * full_pairs, full_pairs)


def histogram_log_probability(histogram, probabilities=None):
    """Probability of a fixed histogram under specified independent pair types."""
    total = sum(histogram)
    if total <= 0 or min(histogram) < 0:
        raise ValueError("Invalid multinomial histogram")
    value = (math.lgamma(total + 1)
             - sum(math.lgamma(n + 1) for n in histogram)) / math.log(2.)
    if probabilities is None:
        probabilities = tuple(n / total for n in histogram)
    if (len(probabilities) != len(histogram)
            or any(p < 0 or (n and p == 0) for n, p in zip(histogram, probabilities))
            or abs(sum(probabilities)-1) > 1e-12):
        raise ValueError("Invalid IID proposal for the required histogram")
    value += sum(n * math.log2(p) for n, p in zip(histogram, probabilities) if n)
    return min(0., value)  # Binary64 diagnostic only; not directed rounding.


def regional_log_matrix(table, group_count, full_pairs, alpha, fractional="packet",
                        tilt=False):
    """IID positive matrix power divided by the exact-histogram probability.

    Conditioning on that histogram gives the actual regional pair placement
    distribution. Fair singleton orientation disappears from count-only caps.
    Fractional powers are taken entrywise BEFORE the pair-type average.
    """
    if not 0 < alpha <= 1 or fractional not in ("packet", "region"):
        raise ValueError("Expected alpha in (0,1] and a supported fractional grouping")
    histogram = pair_histogram(group_count, full_pairs)
    # The region version first conditions the unpowered regional matrix,
    # then uses E[R^alpha] <= (E[R])^alpha entrywise. It loses rare-route
    # protection within a region but avoids fractionalization at every batch.
    local_alpha = alpha if fractional == "packet" else 1.
    result = np.full((2, 2), math.inf)
    for exponent in ((0., .5, 1., 1.5, 2.) if tilt else (0.,)):
        probabilities = np.asarray(tuple(n / PAIRS for n in histogram))
        probabilities *= np.exp2(8. * local_alpha * exponent * np.arange(3))
        probabilities /= probabilities.sum()
        weights = profile_weights(probabilities)
        if not abs(float(weights.sum()) - 1.) < 1e-10:
            raise ArithmeticError("Four-block occupancy weights do not normalize")
        positive = weights > 0
        terms = np.log2(weights[positive, None, None]) + local_alpha * table[positive]
        average = np.logaddexp2.reduce(terms, axis=0)
        candidate = log_matrix_power(average, BATCHES_PER_REGION)
        candidate -= histogram_log_probability(histogram, probabilities)
        if fractional == "region":
            candidate *= alpha
        # Each proposal bounds the SAME conditional regional entry. Taking
        # their entrywise minimum remains a valid positive upper matrix.
        result = np.minimum(result, candidate)
    return result


def source_factor(group_count, source_mode="cutoff"):
    """Return a FLOAT source log using the cutoff or joint-Jensen bound.

    beta is the sharper invertible-Toeplitz density, not the old default.
    Gamma<2 costs one bit. The binomial choice of actual groups is separate.
    """
    if (not isinstance(group_count, int) or not 1 <= group_count <= GROUPS
            or source_mode not in ("cutoff", "joint")):
        raise ValueError("Invalid actual group count")
    log_beta = 1280. - 32. * math.log2(255.)
    ratio = 2048. - log_beta
    if source_mode == "joint":
        # Section 8 cancels kappa^J on each exact zero-group stratum.
        # There is no missing-group cutoff or retained-fraction denominator.
        correction = (group_count / Q0) * ratio
        return dict(source_mode=source_mode,
                    log2_prefactor=1. + group_count * log_beta + correction - math.log2(255.),
                    log2_missing_factor=correction, density_mode="invertible-toeplitz",
                    missingfactor_mode="joint-jensen-exact-zero-group-strata",
                    origin_eta="1/256", log2_outer_density=log_beta, gamma_log2_cap=1.)
    first = group_count // 256
    candidates = sorted({first, first + 1, group_count // 128,
                         (group_count + 63) // 64})
    records = []
    for missing in candidates:
        if not 0 <= missing < group_count:
            continue
        fraction = Fraction(1) - Fraction(group_count, Q0 * (missing + 1))
        if fraction <= 0:
            continue
        # Exact missing-group sets are disjoint zero-group strata in the
        # completed input space. Their densities are bounded pointwise by
        # beta^g*kappa^h. The older binomial sum is valid but unnecessary.
        log_correction = missing * ratio
        value = (1. + group_count * log_beta + log_correction
                 - math.log2(255.) - math.log2(float(fraction)))
        records.append((value, missing, str(fraction), log_correction))
    if not records:
        raise ArithmeticError("No legal missing-shadow witness")
    value, missing, fraction, correction = min(records)
    return dict(source_mode=source_mode, log2_prefactor=value,
                missing_groups=missing, retained_fraction=fraction,
                log2_missing_factor=correction, density_mode="invertible-toeplitz",
                missingfactor_mode="exact-zero-group-strata",
                log2_outer_density=log_beta, gamma_log2_cap=1.,
                candidate_missing_groups=[row[1] for row in records])


def selected_profiles(group_count):
    low, high = max(0, group_count - PAIRS), group_count // 2
    counts = sorted({low, (low + high) // 2, high})
    return tuple((left, right) for left, right in product(counts, repeat=2)
                 if group_count == GROUPS or left + right < group_count)


def selfcheck():
    """Tiny coefficient and ordered-matrix arithmetic checks, not a proof run."""
    coefficient_cases = 0
    probability_cases = ((Fraction(1, 2), Fraction(1, 3), Fraction(1, 6)),
                         (Fraction(0), Fraction(1, 2), Fraction(1, 2)),
                         (Fraction(1), Fraction(0), Fraction(0)))
    for probabilities in probability_cases:
        for pairs in range(5):
            exact = [[Fraction(0) for _ in range(pairs+1)] for _ in range(pairs+1)]
            oriented = (probabilities[0], probabilities[1]/2,
                        probabilities[1]/2, probabilities[2])
            for sequence in product(range(4), repeat=pairs):
                left = sum(j in (1, 3) for j in sequence)
                right = sum(j in (2, 3) for j in sequence)
                exact[left][right] += math.prod(oriented[j] for j in sequence)
            got = occupancy_distribution(tuple(map(float, probabilities)), pairs)
            expected = np.asarray([[float(x) for x in row] for row in exact])
            if sum(map(sum, exact)) != 1 or not np.allclose(got, expected, rtol=2e-14, atol=2e-15):
                raise ArithmeticError("Tiny IID pair polynomial check failed")
            coefficient_cases += 1
    matrices = (np.array(((.8, .2), (.1, .9))),
                np.array(((1., 0.), (.125, .5))),
                np.array(((.1, .7), (.3, .2))))
    matrix_cases = 0
    for left in matrices:
        with np.errstate(divide="ignore"):
            logs = np.log2(left)
        for exponent in range(9):
            got = np.exp2(log_matrix_power(logs, exponent))
            if not np.allclose(got, np.linalg.matrix_power(left, exponent),
                               rtol=2e-13, atol=2e-15):
                raise ArithmeticError("Tiny log-matrix power check failed")
            matrix_cases += 1
        for right in matrices:
            with np.errstate(divide="ignore"):
                got = np.exp2(log_matrix_product(logs, np.log2(right)))
            if not np.allclose(got, left @ right, rtol=2e-14, atol=2e-15):
                raise ArithmeticError("Tiny ordered log-matrix product check failed")
            matrix_cases += 1
    for probabilities in probability_cases:
        mass = float(profile_weights(tuple(map(float, probabilities))).sum())
        if not abs(mass - 1.) < 2e-13:
            raise ArithmeticError("Sorted occupancy multiplicities do not normalize")
    singleton = occupancy_distribution((0., 1., 0.))
    for left, right in product(range(17), repeat=2):
        expected = math.comb(16, left)/2**16 if left+right == 16 else 0.
        if singleton[left, right] != expected:
            raise ArithmeticError("All-singleton striped-layout regression failed")
    profiles, _ = occupancy_profiles()
    balanced = np.flatnonzero(np.all(profiles == (8, 8, 8, 8), axis=1))[0]
    if profile_weights((0., 1., 0.))[balanced] != (math.comb(16, 8)/2**16)**2:
        raise ArithmeticError("Incorrect independent-whole-pair layout was reintroduced")
    return dict(exact_pair_polynomials=coefficient_cases,
                ordered_matrix_arithmetic=matrix_cases, profile_mass_checks=len(probability_cases),
                striped_singleton_cells=289, wrong_layout_negative_control=True)


def _resident_bytes():
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
    maximum = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
    return int(maximum if sys.platform == "darwin" else 1024 * maximum)


def screen(groups=DEFAULT_GROUPS, quick=False, seconds=60., memory_mib=1024,
           fractional="packet", tilt=False, origin_penalty=False, source_mode="cutoff"):
    if (not 0 < seconds <= 60 or not 32 <= memory_mib <= 1024
            or not groups or len(groups) > len(DEFAULT_GROUPS)
            or len(set(groups)) != len(groups)
            or any(not isinstance(g, int) or not 1 <= g <= GROUPS for g in groups)
            or fractional not in ("packet", "region")
            or source_mode not in ("cutoff", "joint")):
        raise ValueError("Use at most seven distinct sizes, at most 60 seconds and 1024 MiB")
    if shutil.disk_usage(Path(__file__).resolve().parent).free < 20 * 1024 ** 3:
        raise RuntimeError("At least 20 GiB free space is required")
    started = time.monotonic()
    deadline = started + seconds
    checks = selfcheck()

    def guard():
        if time.monotonic() >= deadline:
            raise TimeoutError("Bounded floating screen reached its time cap")
        if _resident_bytes() > memory_mib * 1024 ** 2:
            raise MemoryError("Bounded floating screen reached its memory cap")

    markers = (".9", ".5", ".1") if quick else DEFAULT_MARKERS
    alphas = (.1, .5, 1.) if quick else DEFAULT_ALPHAS
    sources = {g: source_factor(g, source_mode) for g in groups}
    profiles = {g: selected_profiles(g) for g in groups}
    best = {}
    evaluations = 0
    status = "complete-selected-grid"
    try:
        for mark in markers:
            guard()
            z = Fraction(mark)
            table = local_table(z)
            guard()
            for alpha in alphas:
                for g in groups:
                    guard()
                    base = (log2_binomial(GROUPS, g)
                            + alpha * (sources[g]["log2_prefactor"]
                                       - LIMIT * math.log2(float(z))))
                    for penalty in ((0., 8., 16., 32.) if origin_penalty else (0.,)):
                        counts = sorted({f for pair in profiles[g] for f in pair})
                        # Joint Jensen uses eta=1/256. The helper's effective
                        # fraction 1 is not a retained-source-fraction claim.
                        penalty_fraction = (1. if source_mode == "joint" else
                                            float(Fraction(sources[g]["retained_fraction"])))
                        weighted = penalized_table(table, penalty, penalty_fraction)
                        regions = {f: regional_log_matrix(weighted, g, f, alpha, fractional, tilt)
                                   for f in counts}
                        for f0, f1 in profiles[g]:
                            run = log_matrix_power(log_matrix_product(regions[f0], regions[f1]), 128)
                            operator = log2_sum(run[0])
                            value = base + operator
                            evaluations += 1
                            key = (g, f0, f1)
                            if key not in best or value < best[key]["log2_uncapped_bound"]:
                                best[key] = dict(source_mode=source_mode, actual_groups=g,
                                    even_full_pairs=f0, odd_full_pairs=f1,
                                    z=str(z), alpha=alpha, log2_uncapped_bound=value,
                                    log2_operator=operator, log2_prefactor_and_source_choices=base,
                                    origin_log2_penalty=penalty,
                                    missing_witness=sources[g],
                                    regional_histogram_log_probabilities=[histogram_log_probability(
                                        pair_histogram(g, f)) for f in (f0, f1)],
                                    directed=False, whole_code_certificate=False,
                                    scope="one selected actual pair profile; binomial source-set overcount")
                        guard()
                    guard()
    except (TimeoutError, MemoryError) as error:
        status = str(error)
    rows = [best[key] for key in sorted(best)]
    summary = dict(schema="small-state-bridge-float-v2-striped", status=status,
        source_mode=source_mode,
        evaluations=evaluations, returned_profiles=len(rows),
        requested_profiles=sum(map(len, profiles.values())), seconds=time.monotonic() - started,
        peak_resident_bytes=_resident_bytes(), source_groups=list(groups),
        output_limit=LIMIT, output_length=1 << 21,
        coefficient_field_bits=32, auxiliary_field_order=Q0,
        terminal="positive f<=M*ones relaxation; no actual final refresh charge",
        fractional_grouping=fractional,
        iid_tilt_exponents=[0., .5, 1., 1.5, 2.] if tilt else [0.],
        origin_log2_penalties=[0., 8., 16., 32.] if origin_penalty else [0.],
        route="IID pair types conditioned on actual regional histogram",
        selfchecks=checks,
        directed=False, whole_code_certificate=False, generated_file_bytes=0)
    return rows, summary


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--groups", nargs="+", type=int, default=list(DEFAULT_GROUPS))
    parser.add_argument("--quick", action="store_true", help="Use only three z and three alpha marks")
    parser.add_argument("--seconds", type=float, default=60.)
    parser.add_argument("--memory-mib", type=int, default=1024)
    parser.add_argument("--fractional", choices=("packet", "region"), default="packet")
    parser.add_argument("--tilt", action="store_true", help="Try five positive coefficient tilts")
    parser.add_argument("--origin-penalty", action="store_true",
                        help="Try the conditional-Jensen penalty on auxiliary zero live states")
    parser.add_argument("--source-mode", choices=("cutoff", "joint"), default="cutoff",
                        help="Use the missing-group cutoff or joint-Jensen source bound")
    args = parser.parse_args()
    rows, summary = screen(tuple(args.groups), args.quick, args.seconds, args.memory_mib,
                           args.fractional, args.tilt, args.origin_penalty, args.source_mode)
    for row in rows:
        print(json.dumps(row, allow_nan=False), flush=True)
    print(json.dumps(summary, allow_nan=False), flush=True)


if __name__ == "__main__":
    main()
