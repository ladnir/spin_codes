#!/usr/bin/env python3
"""Outward/rational certificate for every RM/EBCH post-prefix row family.

All combinatorial inputs and monotonicity checks are exact rational or integer
computations. Transcendental reporting is isolated in the standard-library
outward interval layer.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
from dataclasses import asdict, dataclass
from decimal import Decimal
from fractions import Fraction
from pathlib import Path

from bound_fullsplit_episode_gaps import load_spectrum
from certify_rm_outer_prefix_exact import direct_sum_coefficients, load_local_spectrum
from outward_log2 import Interval, log2_binom, log2_fraction, log2_int, self_check


ROOT = Path(__file__).resolve().parent
DEFAULT_ARTIFACT = ROOT / "fullsplit_postprefix_rational.json"
N = 2**21
B = 64
OUTER_BLOCKS = 4096
LATE_BLOCKS = 5949
DISTANCE = (9 * N) // 100


def exact_split_data(inner_spectrum: Path) -> tuple[list[Fraction], list[tuple[int, int, Fraction]]]:
    spectrum = load_spectrum(inner_spectrum)
    nonzero_total = sum(count for weight, count in spectrum if weight > 0)
    q_law = [Fraction(0) for _ in range(B + 1)]
    entries: list[tuple[int, int, Fraction]] = []
    for weight, count in spectrum:
        if weight <= 0 or count <= 0:
            continue
        den = nonzero_total * math.comb(2 * B, weight)
        for q in range(max(0, weight - B), min(B, weight) + 1):
            j = weight - q
            probability = Fraction(count * math.comb(B, q) * math.comb(B, j), den)
            q_law[q] += probability
            entries.append((j, q, probability))
    if sum(q_law, Fraction(0)) != 1:
        raise SystemExit("postprefix rational: split law does not sum to one")
    return q_law, entries


def certify_global_turnoff(q_law: list[Fraction]) -> Fraction:
    for q in range(1, B + 1):
        if q_law[q] * math.comb(B, B) > q_law[B] * math.comb(B, q):
            raise SystemExit("postprefix rational: Bernstein maximum is not q=b")
    n_min = B * LATE_BLOCKS
    correction = Fraction(n_min, n_min - (B - 1)) ** B
    envelope = q_law[0] + q_law[B] * correction
    if envelope.numerator << 62 > envelope.denominator:
        raise SystemExit("postprefix rational: global turnoff envelope exceeds 2^-62")
    return envelope


def exact_mgf(entries: list[tuple[int, int, Fraction]], pole: Fraction) -> Fraction:
    value = sum(
        (probability * pole**j for j, q, probability in entries if q > 0),
        Fraction(0),
    )
    if not 0 < value < 1:
        raise SystemExit("postprefix rational: invalid rational MGF")
    return value


def load_full_local_spectrum(path: Path) -> list[tuple[int, int]]:
    with path.open(newline="") as handle:
        return [(int(row["weight"]), int(row["count"])) for row in csv.DictReader(handle)]


def local_weight_enumerator(spectrum: list[tuple[int, int]], z: Fraction) -> Fraction:
    return sum((count * z**weight for weight, count in spectrum), Fraction(0))


def upper_sum_log2(max_term: Decimal, count: int) -> Decimal:
    return max_term + log2_int(count).hi


@dataclass(frozen=True)
class PoleData:
    pole: Fraction
    rho: Fraction
    mgf: Fraction
    G: Fraction


def pole_data(entries: list[tuple[int, int, Fraction]], pole: Fraction, rho: Fraction) -> PoleData:
    mgf = exact_mgf(entries, pole)
    G = ((1 + rho) ** B - 1) * mgf / (1 - mgf)
    if G <= 1:
        raise SystemExit("postprefix rational: selected occupancy pole requires G>1")
    return PoleData(pole, rho, mgf, G)


def inner_branch_upper_logs(
    T: int, H: int, pole: PoleData, *, tail_T: int | None = None
) -> tuple[Decimal, Decimal]:
    pref = log2_fraction(pole.pole).times_int(-DISTANCE)
    e0 = pref + log2_fraction(pole.mgf).times_int(T)

    x_max = min(H, T)
    # sum_{x<=x_max}(x+1)G^x <= (x_max+1)G^(x_max+1)/(G-1), G>1.
    log_s = (
        log2_int(x_max + 1)
        + log2_fraction(pole.G).times_int(x_max + 1)
        - log2_fraction(pole.G - 1)
    )
    tail_T = T if tail_T is None else tail_T
    q = Fraction(1, 1 << 62) * (1 - pole.mgf)
    first_tail = Fraction(H * tail_T, 2) * q
    next_ratio = q * Fraction(max(0, H - 1) * max(0, tail_T - 1), 6)
    if next_ratio >= 1:
        raise SystemExit("postprefix rational: H,T tail ratio is not geometric")
    tail_factor = 1 + first_tail / (1 - next_ratio)
    ege1 = (
        Interval.exact(-62)
        + pref
        - log2_binom(B * T, H)
        - log2_fraction(pole.rho).times_int(H)
        + log_s
        + log2_fraction(tail_factor)
    )
    return e0.hi, ege1.hi


@dataclass(frozen=True)
class CriticalWindowCertificate:
    h_min: int
    h_max: int
    live_outer_weights: int
    cap_exact_log2_upper: float
    episode_bucket_log2_upper: tuple[float, ...]
    total_log2_upper: float
    threshold_shift: int
    status: str


def certify_critical_window(
    *, local_spectrum_csv: Path, inner_spectrum: Path, threshold_shift: int = 180
) -> CriticalWindowCertificate:
    h_min, h_max = 501, 2000
    outer = direct_sum_coefficients(
        load_local_spectrum(local_spectrum_csv, h_max), blocks=OUTER_BLOCKS, h_max=h_max
    )
    live_h = [h for h in range(h_min, h_max + 1) if outer[h]]
    q_law, entries = exact_split_data(inner_spectrum)
    certify_global_turnoff(q_law)
    low_pole = pole_data(entries, Fraction(997, 1000), Fraction(3, 1000))
    high_pole = pole_data(entries, Fraction(99, 100), Fraction(1, 100))

    cap_total = Fraction(0)
    for h in live_h:
        numerator = 0
        for r in range(1, B + 1):
            H = h - r
            numerator += math.comb(B, r) * 4000 * math.comb(B * 9948, H)
        cap_total += Fraction(outer[h] * numerator, math.comb(N, h))
    cap_upper = log2_fraction(cap_total).hi

    bucket_bounds: list[Decimal] = []
    episode_buckets = (
        (9949, 13948, 4000),
        (13949, 17948, 4000),
        (17949, 22948, 5000),
        (22949, 27948, 5000),
        (27949, 32767, 4819),
    )
    for T_inner, T_placement, gap_size in episode_buckets:
        peak = Decimal("-Infinity")
        for h in live_h:
            outer_log = log2_int(outer[h])
            denominator = log2_binom(N, h)
            for r in range(1, B + 1):
                H = h - r
                selected = low_pole if H <= 1100 else high_pole
                e0, ege1 = inner_branch_upper_logs(
                    T_inner, H, selected, tail_T=T_placement
                )
                placement = (
                    log2_int(gap_size)
                    + log2_int(math.comb(B, r))
                    + log2_binom(B * T_placement, H)
                    - denominator
                )
                base = outer_log + placement
                peak = max(peak, base.hi + e0, base.hi + ege1)
        bucket_bounds.append(upper_sum_log2(peak, len(live_h) * B * 2))

    total_upper = max(cap_upper, *bucket_bounds) + log2_int(1 + len(bucket_bounds)).hi
    if total_upper > Decimal(-threshold_shift):
        raise SystemExit(
            f"postprefix rational: critical total {total_upper} exceeds 2^-{threshold_shift}"
        )
    return CriticalWindowCertificate(
        h_min=h_min,
        h_max=h_max,
        live_outer_weights=len(live_h),
        cap_exact_log2_upper=float(cap_upper),
        episode_bucket_log2_upper=tuple(float(value) for value in bucket_bounds),
        total_log2_upper=float(total_upper),
        threshold_shift=threshold_shift,
        status="PASS",
    )


LATE_INTERVALS = (
    (501, 2000),
    (2001, 7858),
    (7859, 20550),
    (20551, 75000),
    (75001, 150000),
    (150001, 250000),
    (250001, 380736),
)


def ratio_above_one(h: int, *, m: int, z: Fraction) -> bool:
    return z.denominator * (m - h) > z.numerator * (N - h)


def unimodal_peak(lo: int, hi: int, *, m: int, z: Fraction) -> int:
    if lo >= hi:
        return lo
    left, right = lo, hi
    while left < right:
        mid = (left + right) // 2
        if ratio_above_one(mid, m=m, z=z):
            left = mid + 1
        else:
            right = mid
    return left


@dataclass(frozen=True)
class LatePlacementCertificate:
    interval_log2_upper: tuple[float, ...]
    total_log2_upper: float
    threshold_shift: int
    z_numerator: int
    z_denominator: int
    status: str


def certify_late_placement(
    *, local_spectrum_csv: Path, threshold_shift: int = 500
) -> LatePlacementCertificate:
    z = Fraction(39605985943459426, 10**17)
    local = local_weight_enumerator(load_full_local_spectrum(local_spectrum_csv), z)
    global_outer = log2_fraction(local).times_int(OUTER_BLOCKS)
    log_z = log2_fraction(z)
    m = B * LATE_BLOCKS
    interval_bounds: list[Decimal] = []
    for lo, hi in LATE_INTERVALS:
        peak_h = unimodal_peak(lo, hi, m=m, z=z)
        row = (
            global_outer
            - log_z.times_int(peak_h)
            + log2_binom(m, peak_h)
            - log2_binom(N, peak_h)
        )
        interval_bounds.append(upper_sum_log2(row.hi, hi - lo + 1))
    total_upper = max(interval_bounds) + log2_int(len(interval_bounds)).hi
    if total_upper > Decimal(-threshold_shift):
        raise SystemExit(
            f"postprefix rational: late-placement total {total_upper} exceeds 2^-{threshold_shift}"
        )
    return LatePlacementCertificate(
        interval_log2_upper=tuple(float(value) for value in interval_bounds),
        total_log2_upper=float(total_upper),
        threshold_shift=threshold_shift,
        z_numerator=z.numerator,
        z_denominator=z.denominator,
        status="PASS",
    )


COMPLEMENT_INTERVALS = (
    # h_lo, h_hi, rho, z_e0, z_e>=1
    (1048577, 1148736, "1", "0.999", "0.918412512474"),
    (1148737, 1300000, "1", "0.844325868938", "0.713599772731"),
    (1300001, 1500000, "1", "0.603113862047", "0.509734370012"),
    (1500001, 1700000, "1", "0.430812727615", "0.430812727615"),
    (1700001, 1900000, "10", "0.430812727615", "0.430812727615"),
    (1900001, 2097089, "10", "0.396059859435", "0.364110441035"),
    (2097090, 2097152, "10", "0.364110441035", "0.000113389878138"),
)


def decimal_fraction(text: str) -> Fraction:
    return Fraction(Decimal(text))


def outer_gf_log(
    spectrum: list[tuple[int, int]], z: Fraction, cache: dict[Fraction, Interval]
) -> Interval:
    cached = cache.get(z)
    if cached is None:
        cached = log2_fraction(local_weight_enumerator(spectrum, z)).times_int(OUTER_BLOCKS)
        cache[z] = cached
    return cached


@dataclass(frozen=True)
class ComplementHighCertificate:
    interval_log2_upper: tuple[float, ...]
    total_log2_upper: float
    threshold_shift: int
    pole_numerator: int
    pole_denominator: int
    status: str


def certify_complement_high(
    *,
    local_spectrum_csv: Path,
    inner_spectrum: Path,
    threshold_shift: int = 30000,
) -> ComplementHighCertificate:
    spectrum = load_full_local_spectrum(local_spectrum_csv)
    spectrum_map = dict(spectrum)
    for weight, count in spectrum:
        if spectrum_map.get(512 - weight) != count:
            raise SystemExit("postprefix rational: local RM spectrum is not complement symmetric")
    _q_law, entries = exact_split_data(inner_spectrum)
    pole = Fraction(1, 10)
    mgf = exact_mgf(entries, pole)
    pref = log2_fraction(pole).times_int(-DISTANCE)
    log_mgf = log2_fraction(mgf)
    pterm = Fraction(1, 1 << 62)
    T_max = 32767
    n_endpoint = B * T_max
    outer_cache: dict[Fraction, Interval] = {}
    interval_bounds: list[Decimal] = []

    for h_lo, h_hi, rho_text, z0_text, z1_text in COMPLEMENT_INTERVALS:
        rho = decimal_fraction(rho_text)
        z0 = decimal_fraction(z0_text)
        z1 = decimal_fraction(z1_text)
        count_log = log2_int(h_hi - h_lo + 1)
        H_min = h_lo - B
        H_max = h_hi - 1

        outer_e0 = outer_gf_log(spectrum, z0, outer_cache) - log2_fraction(z0).times_int(N - h_lo)
        ratio_e0 = Fraction(1, 1) / mgf * Fraction(n_endpoint - H_min, n_endpoint) ** B
        if ratio_e0 >= 1:
            raise SystemExit("postprefix rational: complement e0 endpoint ratio is not below one")
        e0 = (
            count_log
            + outer_e0
            + pref
            + log_mgf.times_int(T_max)
            + log2_fraction(Fraction(1, 1) / (1 - ratio_e0))
        )

        outer_base = outer_gf_log(spectrum, z1, outer_cache)
        log_z1 = log2_fraction(z1)
        log_rho = log2_fraction(rho)
        endpoint_values = []
        for h in (h_lo, h_hi):
            endpoint_values.append(
                outer_base
                - log_z1.times_int(N - h)
                - log2_binom(N, h)
                - log_rho.times_int(h)
            )
        outer_volume_hi = max(value.hi for value in endpoint_values)

        G = ((1 + rho) ** B - 1) * mgf / (1 - mgf)
        if G <= 2:
            raise SystemExit("postprefix rational: complement paired ratio requires G>2")
        t_overhead = (G - 1) / (G - 2)
        r_sum = (1 + rho) ** B - 1
        q = pterm * (1 - mgf)
        first_tail = Fraction(H_max * T_max, 2) * q
        next_ratio = q * Fraction((H_max - 1) * (T_max - 1), 6)
        if next_ratio >= 1:
            raise SystemExit("postprefix rational: complement tail is not geometric")
        tail_factor = 1 + first_tail / (1 - next_ratio)
        ege1_hi = (
            count_log.hi
            + outer_volume_hi
            + log2_fraction(r_sum).hi
            - Decimal(62)
            + pref.hi
            + log2_int(T_max + 1).hi
            + (T_max + 1) * log2_fraction(G).hi
            - log2_fraction(G - 1).lo
            + log2_fraction(tail_factor).hi
            + log2_fraction(t_overhead).hi
        )
        interval_bounds.append(max(e0.hi, ege1_hi) + log2_int(2).hi)

    total_upper = max(interval_bounds) + log2_int(len(interval_bounds)).hi
    if total_upper > Decimal(-threshold_shift):
        raise SystemExit(
            f"postprefix rational: complement total {total_upper} exceeds 2^-{threshold_shift}"
        )
    return ComplementHighCertificate(
        interval_log2_upper=tuple(float(value) for value in interval_bounds),
        total_log2_upper=float(total_upper),
        threshold_shift=threshold_shift,
        pole_numerator=pole.numerator,
        pole_denominator=pole.denominator,
        status="PASS",
    )


EARLY_ENDPOINT_INTERVALS = (
    # h_lo, h_hi, rational Chernoff pole, rho, outer z
    (2001, 7858, Fraction(9801986733, 10**10), Fraction(1, 100), "0.396059859435"),
    (7859, 20550, Fraction(9512294245, 10**10), Fraction(3, 100), "0.396059859435"),
    (20551, 30000, Fraction(8187307531, 10**10), Fraction(1, 10), "0.396059859435"),
    (30001, 50000, Fraction(8187307531, 10**10), Fraction(1, 10), "0.396059859435"),
    (50001, 75000, Fraction(8187307531, 10**10), Fraction(1, 10), "0.396059859435"),
    (75001, 90000, Fraction(8187307531, 10**10), Fraction(1, 10), "0.396059859435"),
    (90001, 100000, Fraction(8187307531, 10**10), Fraction(1, 10), "0.396059859435"),
    (100001, 110000, Fraction(6065306597, 10**10), Fraction(3, 10), "0.396059859435"),
    (110001, 125000, Fraction(6065306597, 10**10), Fraction(3, 10), "0.396059859435"),
    (125001, 160000, Fraction(4493289641, 10**10), Fraction(3, 10), "0.396059859435"),
    (160001, 250000, Fraction(3011942119, 10**10), Fraction(1, 2), "0.396059859435"),
    (250001, 350000, Fraction(3011942119, 10**10), Fraction(1, 2), "0.430812727615"),
)

EARLY_GAPS = (
    (17949, 22948, 5000),
    (22949, 27948, 5000),
    (27949, 32767, 4819),
)


def common_ratio_above_one(h: int, *, r: int, n_end: int, z: Fraction) -> bool:
    H = h - r
    return (
        z.denominator * (n_end - H) * (h + 1)
        > z.numerator * (H + 1) * (N - h)
    )


def common_peak(lo: int, hi: int, *, r: int, n_end: int, z: Fraction) -> int:
    left, right = lo, hi
    while left < right:
        mid = (left + right) // 2
        if common_ratio_above_one(mid, r=r, n_end=n_end, z=z):
            left = mid + 1
        else:
            right = mid
    return left


def common_row_log(
    *,
    h: int,
    r: int,
    n_end: int,
    gap_size: int,
    z: Fraction,
    outer_log: Interval,
) -> Interval:
    H = h - r
    return (
        outer_log
        - log2_fraction(z).times_int(h)
        - log2_binom(N, h)
        + log2_int(math.comb(B, r))
        + log2_int(gap_size)
        + log2_binom(n_end, H)
    )


@dataclass(frozen=True)
class EarlyEndpointCertificate:
    interval_log2_upper: tuple[float, ...]
    total_log2_upper: float
    threshold_shift: int
    status: str


def certify_early_endpoint(
    *,
    local_spectrum_csv: Path,
    inner_spectrum: Path,
    threshold_shift: int = 400,
    intervals: tuple[tuple[int, int, Fraction, Fraction, str], ...] = EARLY_ENDPOINT_INTERVALS,
    gaps: tuple[tuple[int, int, int], ...] = EARLY_GAPS,
) -> EarlyEndpointCertificate:
    spectrum = load_full_local_spectrum(local_spectrum_csv)
    _q_law, entries = exact_split_data(inner_spectrum)
    outer_cache: dict[Fraction, Interval] = {}
    pole_cache: dict[tuple[Fraction, Fraction], PoleData] = {}
    interval_bounds: list[Decimal] = []

    for h_lo, h_hi, pole_value, rho, z_text in intervals:
        z = decimal_fraction(z_text)
        outer_log = outer_gf_log(spectrum, z, outer_cache)
        pdata = pole_cache.get((pole_value, rho))
        if pdata is None:
            mgf = exact_mgf(entries, pole_value)
            G = ((1 + rho) ** B - 1) * mgf / (1 - mgf)
            pdata = PoleData(pole_value, rho, mgf, G)
            pole_cache[(pole_value, rho)] = pdata
        pref = log2_fraction(pdata.pole).times_int(-DISTANCE)
        peak = Decimal("-Infinity")

        for T_min, T_max, gap_size in gaps:
            n_min, n_end = B * T_min, B * T_max
            e0_inner = pref + log2_fraction(pdata.mgf).times_int(T_min)
            H_max = h_hi - 1
            q = Fraction(1, 1 << 62) * (1 - pdata.mgf)
            first_tail = Fraction(H_max * T_max, 2) * q
            next_ratio = q * Fraction((H_max - 1) * (T_max - 1), 6)
            if next_ratio >= 1:
                raise SystemExit("postprefix rational: early endpoint tail is not geometric")
            tail_factor = 1 + first_tail / (1 - next_ratio)
            inner_base = (
                Interval.exact(-62)
                + pref
                + log2_fraction(tail_factor)
            )

            for r in range(1, B + 1):
                peak_h = common_peak(h_lo, h_hi, r=r, n_end=n_end, z=z)
                common_e0 = common_row_log(
                    h=peak_h,
                    r=r,
                    n_end=n_end,
                    gap_size=gap_size,
                    z=z,
                    outer_log=outer_log,
                )
                peak = max(peak, common_e0.hi + e0_inner.hi)

                # After division by C(n_min,H) and rho^H the adjacent row
                # ratio is increasing, so the maximum is at an h endpoint.
                for h in (h_lo, h_hi):
                    H = h - r
                    if pdata.G > 1:
                        x_max = min(H, T_min)
                        log_s = (
                            log2_int(x_max + 1)
                            + log2_fraction(pdata.G).times_int(x_max + 1)
                            - log2_fraction(pdata.G - 1)
                        )
                    else:
                        log_s = log2_fraction(1 - pdata.G).times_int(-2)
                    common = common_row_log(
                        h=h,
                        r=r,
                        n_end=n_end,
                        gap_size=gap_size,
                        z=z,
                        outer_log=outer_log,
                    )
                    ege1 = (
                        common
                        + inner_base
                        + log_s
                        - log2_binom(n_min, H)
                        - log2_fraction(rho).times_int(H)
                    )
                    peak = max(peak, ege1.hi)

        count = (h_hi - h_lo + 1) * B * len(gaps) * 2
        interval_bounds.append(upper_sum_log2(peak, count))

    total_upper = max(interval_bounds) + log2_int(len(interval_bounds)).hi
    if total_upper > Decimal(-threshold_shift):
        raise SystemExit(
            f"postprefix rational: early endpoint total {total_upper} exceeds 2^-{threshold_shift}; "
            f"intervals={interval_bounds}"
        )
    return EarlyEndpointCertificate(
        interval_log2_upper=tuple(float(value) for value in interval_bounds),
        total_log2_upper=float(total_upper),
        threshold_shift=threshold_shift,
        status="PASS",
    )


EARLY_PAIRED_INTERVALS = (
    # h_lo, h_hi, z_e0, z_e>=1
    (350001, 524288, "0.430812727615", "0.430812727615"),
    (524289, 750000, "0.55446177913", "0.509734370012"),
    (750001, 1048576, "0.999", "0.918412512474"),
)

PREFIX_ENDPOINT_INTERVALS = (
    (2001, 7858, Fraction(9801986733, 10**10), Fraction(1, 100), "0.396059859435"),
    (7859, 20550, Fraction(9512294245, 10**10), Fraction(3, 100), "0.396059859435"),
    (20551, 75000, Fraction(8187307531, 10**10), Fraction(1, 10), "0.396059859435"),
    (75001, 250000, Fraction(6065306597, 10**10), Fraction(3, 10), "0.396059859435"),
    (250001, 350000, Fraction(3011942119, 10**10), Fraction(3, 4), "0.396059859435"),
)

PREFIX_EPISODE_GAPS = (
    (9949, 13948, 4000),
    (13949, 17948, 4000),
)

PREFIX_CAP_INTERVALS = (
    (2001, 7858, "0.396059859435"),
    (7859, 20550, "0.396059859435"),
    (20551, 75000, "0.396059859435"),
    (75001, 250000, "0.396059859435"),
    (250001, 350000, "0.396059859435"),
    (350001, 400000, "0.430812727615"),
    (400001, 450000, "0.430812727615"),
    (450001, 550000, "0.430812727615"),
    (550001, 636736, "0.430812727615"),
)

PREFIX_PAIRED_RHO1_INTERVALS = (
    (350001, 400000, "0.430812727615", "0.430812727615"),
    (400001, 450000, "0.430812727615", "0.430812727615"),
    (450001, 550000, "0.430812727615", "0.430812727615"),
    (550001, 650000, "0.430812727615", "0.430812727615"),
)
PREFIX_PAIRED_RHO2_INTERVALS = (
    (650001, 725000, "0.430812727615", "0.430812727615"),
)
PREFIX_PAIRED_RHO3_INTERVALS = (
    (725001, 750000, "0.509734370012", "0.509734370012"),
    (750001, 850000, "0.55446177913", "0.55446177913"),
    (850001, 950000, "0.65603499517", "0.65603499517"),
    (950001, 1048576, "0.844325868938", "0.844325868938"),
)


@dataclass(frozen=True)
class PrefixCapCertificate:
    interval_log2_upper: tuple[float, ...]
    total_log2_upper: float
    threshold_shift: int
    status: str


def certify_prefix_cap(
    *, local_spectrum_csv: Path, threshold_shift: int = 700
) -> PrefixCapCertificate:
    spectrum = load_full_local_spectrum(local_spectrum_csv)
    outer_cache: dict[Fraction, Interval] = {}
    m = B * 9948 + B
    interval_bounds: list[Decimal] = []
    for h_lo, h_hi, z_text in PREFIX_CAP_INTERVALS:
        z = decimal_fraction(z_text)
        peak_h = unimodal_peak(h_lo, min(h_hi, m), m=m, z=z)
        row = (
            outer_gf_log(spectrum, z, outer_cache)
            - log2_fraction(z).times_int(peak_h)
            + log2_int(4000)
            + log2_binom(m, peak_h)
            - log2_binom(N, peak_h)
        )
        interval_bounds.append(upper_sum_log2(row.hi, h_hi - h_lo + 1))
    total_upper = max(interval_bounds) + log2_int(len(interval_bounds)).hi
    if total_upper > Decimal(-threshold_shift):
        raise SystemExit(
            f"postprefix rational: prefix cap total {total_upper} exceeds 2^-{threshold_shift}"
        )
    return PrefixCapCertificate(
        interval_log2_upper=tuple(float(value) for value in interval_bounds),
        total_log2_upper=float(total_upper),
        threshold_shift=threshold_shift,
        status="PASS",
    )


def paired_h_peak(lo: int, hi: int, *, m: int, z: Fraction) -> int:
    left, right = lo, min(hi, m)
    if left > right:
        raise SystemExit("postprefix rational: empty paired h range")
    while left < right:
        mid = (left + right) // 2
        if z.denominator * (m - mid) > z.numerator * (N - mid):
            left = mid + 1
        else:
            right = mid
    return left


@dataclass(frozen=True)
class EarlyPairedCertificate:
    interval_log2_upper: tuple[float, ...]
    total_log2_upper: float
    threshold_shift: int
    pole_numerator: int
    pole_denominator: int
    status: str


def certify_early_paired(
    *,
    local_spectrum_csv: Path,
    inner_spectrum: Path,
    threshold_shift: int = 80000,
    intervals: tuple[tuple[int, int, str, str], ...] = EARLY_PAIRED_INTERVALS,
    pole: Fraction = Fraction(1353352832, 10**10),
    rho: Fraction = Fraction(4, 5),
    T_min: int = 17949,
    T_max: int = 32767,
) -> EarlyPairedCertificate:
    spectrum = load_full_local_spectrum(local_spectrum_csv)
    _q_law, entries = exact_split_data(inner_spectrum)
    mgf = exact_mgf(entries, pole)
    G = ((1 + rho) ** B - 1) * mgf / (1 - mgf)
    if G <= 2:
        raise SystemExit("postprefix rational: early paired ratio requires G>2")
    pref = log2_fraction(pole).times_int(-DISTANCE)
    log_mgf = log2_fraction(mgf)
    outer_cache: dict[Fraction, Interval] = {}
    pterm = Fraction(1, 1 << 62)
    r_sum = (1 + rho) ** B - 1
    log_s_full = (
        log2_int(T_max + 1)
        + log2_fraction(G).times_int(T_max + 1)
        - log2_fraction(G - 1)
    )
    t_overhead = (G - 1) / (G - 2)
    interval_bounds: list[Decimal] = []

    for h_lo, h_hi, z0_text, z1_text in intervals:
        z0 = decimal_fraction(z0_text)
        z1 = decimal_fraction(z1_text)
        outer0 = outer_gf_log(spectrum, z0, outer_cache)
        peak_e0 = Decimal("-Infinity")
        for T in range(T_min, T_max + 1):
            m = B * T + B
            if m < h_lo:
                continue
            peak_h = paired_h_peak(h_lo, h_hi, m=m, z=z0)
            row = (
                outer0
                - log2_fraction(z0).times_int(peak_h)
                + log2_binom(m, peak_h)
                - log2_binom(N, peak_h)
                + pref
                + log_mgf.times_int(T)
            )
            peak_e0 = max(peak_e0, row.hi)
        e0_upper = upper_sum_log2(
            peak_e0, (h_hi - h_lo + 1) * (T_max - T_min + 1)
        )

        outer1 = outer_gf_log(spectrum, z1, outer_cache)
        log_z1 = log2_fraction(z1)
        log_rho = log2_fraction(rho)
        outer_volume_hi = Decimal("-Infinity")
        for h in (h_lo, h_hi):
            row = (
                outer1
                - log_z1.times_int(h)
                - log2_binom(N, h)
                - log_rho.times_int(h)
            )
            outer_volume_hi = max(outer_volume_hi, row.hi)
        H_max = h_hi - 1
        q = pterm * (1 - mgf)
        first_tail = Fraction(H_max * T_max, 2) * q
        next_ratio = q * Fraction((H_max - 1) * (T_max - 1), 6)
        if next_ratio >= 1:
            raise SystemExit("postprefix rational: early paired tail is not geometric")
        tail_factor = 1 + first_tail / (1 - next_ratio)
        ege1_upper = (
            log2_int(h_hi - h_lo + 1).hi
            + outer_volume_hi
            + log2_fraction(r_sum).hi
            - Decimal(62)
            + pref.hi
            + log_s_full.hi
            + log2_fraction(tail_factor).hi
            + log2_fraction(t_overhead).hi
        )
        interval_bounds.append(max(e0_upper, ege1_upper) + log2_int(2).hi)

    total_upper = max(interval_bounds) + log2_int(len(interval_bounds)).hi
    if total_upper > Decimal(-threshold_shift):
        raise SystemExit(
            f"postprefix rational: early paired total {total_upper} exceeds 2^-{threshold_shift}"
        )
    return EarlyPairedCertificate(
        interval_log2_upper=tuple(float(value) for value in interval_bounds),
        total_log2_upper=float(total_upper),
        threshold_shift=threshold_shift,
        pole_numerator=pole.numerator,
        pole_denominator=pole.denominator,
        status="PASS",
    )


def denominator_step_ratio(*, T: int, H: int) -> Fraction:
    """C(B(T+1),H)/C(BT,H) as a 64-factor exact rational."""

    n = B * T
    numerator = denominator = 1
    for i in range(1, B + 1):
        numerator *= n + i
        denominator *= n + i - H
    return Fraction(numerator, denominator)


def geometric_tail_factor(*, H: int, T: int, mgf: Fraction) -> Fraction:
    q = Fraction(1, 1 << 62) * (1 - mgf)
    first = Fraction(H * T, 2) * q
    ratio = q * Fraction(max(0, H - 1) * max(0, T - 1), 6)
    if ratio >= 1:
        raise SystemExit("postprefix rational: monotonicity tail is not geometric")
    return 1 + first / (1 - ratio)


def certify_endpoint_t_monotonicity(
    *,
    intervals: tuple[tuple[int, int, Fraction, Fraction, str], ...],
    gaps: tuple[tuple[int, int, int], ...],
    entries: list[tuple[int, int, Fraction]],
) -> int:
    """Exact version of the two-case fixed-pole T reduction."""

    rows = 0
    cache: dict[tuple[Fraction, Fraction], tuple[Fraction, Fraction]] = {}
    for h_lo, h_hi, pole, rho, _z in intervals:
        pair = cache.get((pole, rho))
        if pair is None:
            mgf = exact_mgf(entries, pole)
            G = ((1 + rho) ** B - 1) * mgf / (1 - mgf)
            pair = (mgf, G)
            cache[(pole, rho)] = pair
        mgf, G = pair
        H_min = max(1, h_lo - B)
        H_max = h_hi - 1
        for T_min, T_max, _size in gaps:
            tail_factor = geometric_tail_factor(H=H_max, T=T_max, mgf=mgf)
            T_den = T_max - 1
            if H_min <= min(H_max, T_den):
                if denominator_step_ratio(T=T_den, H=H_min) < tail_factor:
                    raise SystemExit("postprefix rational: no-new-occupancy T reduction failed")
            for T in range(T_min, T_max):
                H = max(H_min, T + 1)
                if H > H_max or H > B * T:
                    continue
                occupancy_growth = 1 + Fraction(T + 2, T + 1) * G
                if denominator_step_ratio(T=T, H=H) < occupancy_growth * tail_factor:
                    raise SystemExit("postprefix rational: new-occupancy T reduction failed")
            rows += 1
    return rows


def certify_critical_t_monotonicity(
    *, entries: list[tuple[int, int, Fraction]]
) -> int:
    poles = (
        exact_mgf(entries, Fraction(997, 1000)),
        exact_mgf(entries, Fraction(99, 100)),
    )
    rows = 0
    for mgf in poles:
        for T_min, T_max, _size in (
            (9949, 13948, 4000),
            (13949, 17948, 4000),
            (17949, 22948, 5000),
            (22949, 27948, 5000),
            (27949, 32767, 4819),
        ):
            tail_factor = geometric_tail_factor(H=1999, T=T_max, mgf=mgf)
            # Here every H is below T_min, so the occupancy sum is unchanged;
            # H=1 is the smallest denominator growth.
            if denominator_step_ratio(T=T_max - 1, H=1) < tail_factor:
                raise SystemExit("postprefix rational: critical T reduction failed")
            rows += 1
    return rows


@dataclass(frozen=True)
class CompletePostprefixCertificate:
    h_min: int
    h_max: int
    component_count: int
    total_log2_upper: float
    threshold_shift: int
    coverage: tuple[str, ...]
    status: str


def combine_complete_postprefix(
    components: list[object], *, threshold_shift: int = 180
) -> CompletePostprefixCertificate:
    totals = [Decimal.from_float(float(getattr(item, "total_log2_upper"))) for item in components]
    # One nanobit absorbs the binary64 serialization of already outward
    # Decimal endpoints; the actual component margins are multiple bits.
    total_upper = max(totals) + Decimal("1e-9") + log2_int(len(totals)).hi
    if total_upper > Decimal(-threshold_shift):
        raise SystemExit(
            f"postprefix rational: complete total {total_upper} exceeds 2^-{threshold_shift}"
        )
    return CompletePostprefixCertificate(
        h_min=501,
        h_max=N,
        component_count=len(components),
        total_log2_upper=float(total_upper),
        threshold_shift=threshold_shift,
        coverage=(
            "501--2000: all non-ultra-late first-active positions",
            "501--380736: ultra-late placement",
            "2001--636736: prefix cap bucket",
            "2001--350000: prefix fixed-T episode buckets",
            "350001--1048576: prefix paired-T episode buckets",
            "2001--350000: early fixed-T buckets",
            "350001--1048576: early paired-T buckets",
            "1048577--2097152: complement-high all-T cover",
        ),
        status="PASS",
    )


def recompute_all(*, local_spectrum_csv: Path, inner_spectrum: Path) -> dict[str, object]:
    self_check()
    q_law, exact_entries = exact_split_data(inner_spectrum)
    global_turnoff = certify_global_turnoff(q_law)
    monotonic_rows = certify_critical_t_monotonicity(entries=exact_entries)
    monotonic_rows += certify_endpoint_t_monotonicity(
        intervals=EARLY_ENDPOINT_INTERVALS,
        gaps=EARLY_GAPS,
        entries=exact_entries,
    )
    monotonic_rows += certify_endpoint_t_monotonicity(
        intervals=PREFIX_ENDPOINT_INTERVALS,
        gaps=PREFIX_EPISODE_GAPS,
        entries=exact_entries,
    )
    critical = certify_critical_window(
        local_spectrum_csv=local_spectrum_csv, inner_spectrum=inner_spectrum
    )
    late = certify_late_placement(local_spectrum_csv=local_spectrum_csv)
    complement = certify_complement_high(
        local_spectrum_csv=local_spectrum_csv,
        inner_spectrum=inner_spectrum,
    )
    early_endpoint = certify_early_endpoint(
        local_spectrum_csv=local_spectrum_csv,
        inner_spectrum=inner_spectrum,
    )
    early_paired = certify_early_paired(
        local_spectrum_csv=local_spectrum_csv,
        inner_spectrum=inner_spectrum,
    )
    prefix_cap = certify_prefix_cap(local_spectrum_csv=local_spectrum_csv)
    prefix_endpoint = certify_early_endpoint(
        local_spectrum_csv=local_spectrum_csv,
        inner_spectrum=inner_spectrum,
        threshold_shift=500,
        intervals=PREFIX_ENDPOINT_INTERVALS,
        gaps=PREFIX_EPISODE_GAPS,
    )
    prefix_paired_rho1 = certify_early_paired(
        local_spectrum_csv=local_spectrum_csv,
        inner_spectrum=inner_spectrum,
        threshold_shift=100000,
        intervals=PREFIX_PAIRED_RHO1_INTERVALS,
        pole=Fraction(3011942119, 10**10),
        rho=Fraction(1, 1),
        T_min=9949,
        T_max=17948,
    )
    prefix_paired_rho2 = certify_early_paired(
        local_spectrum_csv=local_spectrum_csv,
        inner_spectrum=inner_spectrum,
        threshold_shift=200000,
        intervals=PREFIX_PAIRED_RHO2_INTERVALS,
        pole=Fraction(3011942119, 10**10),
        rho=Fraction(2, 1),
        T_min=9949,
        T_max=17948,
    )
    prefix_paired_rho3 = certify_early_paired(
        local_spectrum_csv=local_spectrum_csv,
        inner_spectrum=inner_spectrum,
        threshold_shift=200000,
        intervals=PREFIX_PAIRED_RHO3_INTERVALS,
        pole=Fraction(3011942119, 10**10),
        rho=Fraction(3, 1),
        T_min=9949,
        T_max=17948,
    )
    components = [
        critical,
        late,
        complement,
        early_endpoint,
        early_paired,
        prefix_cap,
        prefix_endpoint,
        prefix_paired_rho1,
        prefix_paired_rho2,
        prefix_paired_rho3,
    ]
    complete = combine_complete_postprefix(components)
    names = (
        "critical_501_2000",
        "late_placement",
        "complement_high",
        "early_endpoint",
        "early_paired",
        "prefix_cap",
        "prefix_endpoint",
        "prefix_paired_rho1",
        "prefix_paired_rho2",
        "prefix_paired_rho3",
    )
    return {
        "schema": 1,
        "arithmetic": "exact rationals plus outward Decimal log2 intervals",
        "decimal_precision": 100,
        "exact_global_turnoff": {
            "status": "PASS",
            "threshold": "2^-62",
            "log2_for_display": float(log2_fraction(global_turnoff).hi),
        },
        "exact_t_monotonicity": {"status": "PASS", "rows": monotonic_rows},
        "components": {name: asdict(item) for name, item in zip(names, components, strict=True)},
        "complete": asdict(complete),
    }


def artifact_payload(data: dict[str, object]) -> dict[str, object]:
    # Round-trip once so regenerated dataclass tuples have the same JSON-native
    # list representation as a payload loaded from disk.
    normalized = json.loads(json.dumps(data))
    canonical = json.dumps(normalized, sort_keys=True, separators=(",", ":")).encode()
    return {**normalized, "sha256": hashlib.sha256(canonical).hexdigest()}


def write_artifact(path: Path, data: dict[str, object]) -> None:
    path.write_text(json.dumps(artifact_payload(data), indent=2) + "\n")


def load_artifact(path: Path) -> dict[str, object]:
    payload = json.loads(path.read_text())
    claimed = payload.pop("sha256")
    canonical = json.dumps(payload, sort_keys=True, separators=(",", ":")).encode()
    actual = hashlib.sha256(canonical).hexdigest()
    if claimed != actual:
        raise SystemExit("postprefix rational: artifact SHA-256 mismatch")
    if payload.get("schema") != 1 or payload.get("decimal_precision") != 100:
        raise SystemExit("postprefix rational: artifact metadata mismatch")
    complete = payload.get("complete")
    if not isinstance(complete, dict) or complete.get("status") != "PASS":
        raise SystemExit("postprefix rational: artifact complete status is not PASS")
    if Decimal(str(complete["total_log2_upper"])) > -Decimal(int(complete["threshold_shift"])):
        raise SystemExit("postprefix rational: artifact complete threshold failed")
    monotonicity = payload.get("exact_t_monotonicity")
    if not isinstance(monotonicity, dict) or monotonicity.get("status") != "PASS":
        raise SystemExit("postprefix rational: artifact T-monotonicity status failed")
    turnoff = payload.get("exact_global_turnoff")
    if not isinstance(turnoff, dict) or turnoff.get("status") != "PASS":
        raise SystemExit("postprefix rational: artifact global turnoff status failed")
    components = payload.get("components")
    if not isinstance(components, dict) or len(components) != 10:
        raise SystemExit("postprefix rational: artifact component coverage mismatch")
    for name, item in components.items():
        if not isinstance(item, dict) or item.get("status") != "PASS":
            raise SystemExit(f"postprefix rational: artifact component {name} is not PASS")
        if Decimal(str(item["total_log2_upper"])) > -Decimal(int(item["threshold_shift"])):
            raise SystemExit(f"postprefix rational: artifact component {name} threshold failed")
    payload["sha256"] = actual
    return payload


def artifact_differences(stored: object, regenerated: object, path: str = "root") -> list[str]:
    """Return compact paths for reproducibility failures."""

    if type(stored) is not type(regenerated):
        return [f"{path}: type {type(stored).__name__} != {type(regenerated).__name__}"]
    if isinstance(stored, dict):
        differences: list[str] = []
        for key in sorted(stored.keys() | regenerated.keys()):
            if key not in stored or key not in regenerated:
                differences.append(f"{path}.{key}: missing from one artifact")
            else:
                differences.extend(artifact_differences(stored[key], regenerated[key], f"{path}.{key}"))
            if len(differences) >= 20:
                break
        return differences
    if isinstance(stored, list):
        if len(stored) != len(regenerated):
            return [f"{path}: length {len(stored)} != {len(regenerated)}"]
        differences = []
        for index, (left, right) in enumerate(zip(stored, regenerated, strict=True)):
            differences.extend(artifact_differences(left, right, f"{path}[{index}]"))
            if len(differences) >= 20:
                break
        return differences
    return [] if stored == regenerated else [f"{path}: {stored!r} != {regenerated!r}"]


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--local-spectrum-csv", type=Path, default=ROOT / "rm512_256_spectrum.csv")
    parser.add_argument("--inner-spectrum", type=Path, default=ROOT / "EBCH128_64.wd")
    parser.add_argument("--artifact", type=Path, default=DEFAULT_ARTIFACT)
    parser.add_argument("--generate-artifact", action="store_true")
    parser.add_argument("--recompute-artifact", action="store_true")
    args = parser.parse_args()
    if args.generate_artifact:
        data = recompute_all(
            local_spectrum_csv=args.local_spectrum_csv,
            inner_spectrum=args.inner_spectrum,
        )
        write_artifact(args.artifact, data)
        print(f"postprefix_rational_artifact_written,{args.artifact}")
        return 0
    if args.recompute_artifact:
        stored = load_artifact(args.artifact)
        regenerated = artifact_payload(
            recompute_all(
                local_spectrum_csv=args.local_spectrum_csv,
                inner_spectrum=args.inner_spectrum,
            )
        )
        if stored != regenerated:
            for difference in artifact_differences(stored, regenerated):
                print(f"postprefix_rational_artifact_difference,{difference}")
            raise SystemExit("postprefix rational: regenerated artifact differs")
        print("postprefix_rational_artifact_recompute_status,PASS")
        print(f"postprefix_rational_artifact_sha256,{stored['sha256']}")
        return 0
    payload = load_artifact(args.artifact)
    complete = payload["complete"]
    print("fullsplit_postprefix_rational_status,PASS")
    print(f"fullsplit_postprefix_rational_log2_upper,{complete['total_log2_upper']:.12f}")
    print(f"fullsplit_postprefix_rational_threshold,2^-{complete['threshold_shift']}")
    print(f"fullsplit_postprefix_rational_artifact_sha256,{payload['sha256']}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
