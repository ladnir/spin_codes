#!/usr/bin/env python3
"""Rational certificate for the complete RM/EBCH ``h <= 500`` family.

The first-active position is split into the ultra-late prefix and six gap
buckets.  The three named prefix buckets retain ``e=0,...,8`` and the three
early buckets retain ``e=0,...,16``, using exact combinatorics and an
outward-rounded rational Chernoff envelope.  The corresponding ``e>=9`` and
``e>=17`` tails are covered.  The final comparison uses integers only;
logarithms are diagnostics.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
from dataclasses import asdict, dataclass
from fractions import Fraction
from pathlib import Path

from certificate_spectra import file_sha256, load_ebch128_spectrum, load_rm512_spectrum
from certify_rm_outer_prefix_exact import direct_sum_coefficients, load_local_spectrum


ROOT = Path(__file__).resolve().parent
DEFAULT_GAP_SUMS = ROOT / "fullsplit_h500_gap_sums_exact.json"
DEFAULT_INNER_BOUNDS = ROOT / "fullsplit_h500_inner_bounds_dyadic.json"
BUCKETS = (
    (1, 4000, 5949),
    (4001, 8000, 9949),
    (8001, 12000, 13949),
    (12001, 17000, 17949),
    (17001, 22000, 22949),
    (22001, 26819, 27949),
)


def ceil_div(num: int, den: int) -> int:
    return -(-num // den)


def log2_int(value: int) -> float:
    if value <= 0:
        return float("-inf")
    shift = max(0, value.bit_length() - 53)
    return math.log2(value >> shift) + shift


def log2_fraction(value: Fraction) -> float:
    return log2_int(value.numerator) - log2_int(value.denominator)


def load_outer_coefficients(local_spectrum_csv: Path, *, blocks: int, h_max: int) -> list[int]:
    local = load_local_spectrum(local_spectrum_csv, h_max)
    return direct_sum_coefficients(local, blocks=blocks, h_max=h_max)


def exact_split_law(inner_spectrum: Path, b: int) -> tuple[list[Fraction], Fraction]:
    spectrum = load_ebch128_spectrum(inner_spectrum)
    nonzero_total = sum(count for weight, count in spectrum if weight > 0)
    split_law = [Fraction(0) for _ in range(b + 1)]
    mgf_terms: list[tuple[int, Fraction]] = []
    z = Fraction(2333, 2373)
    for weight, count in spectrum:
        if weight <= 0 or count <= 0:
            continue
        den = nonzero_total * math.comb(2 * b, weight)
        for state_weight in range(max(0, weight - b), min(b, weight) + 1):
            output_weight = weight - state_weight
            probability = Fraction(
                count * math.comb(b, state_weight) * math.comb(b, output_weight), den
            )
            split_law[state_weight] += probability
            if state_weight > 0:
                mgf_terms.append((output_weight, probability))
    if sum(split_law, Fraction(0)) != 1:
        raise SystemExit("rational h<=500: exact EBCH split law does not sum to one")
    mgf = sum((probability * z**output for output, probability in mgf_terms), Fraction(0))
    if not 0 < mgf < 1:
        raise SystemExit("rational h<=500: invalid fixed-pole MGF")
    return split_law, mgf


def certify_termination_cap(
    split_law: list[Fraction], *, b: int, T_min: int, h_max: int, cap_shift: int
) -> Fraction:
    n = b * T_min
    if n + 1 < b * (h_max + 1):
        raise SystemExit("rational h<=500: termination endpoint monotonicity check failed")
    den = math.comb(n, h_max)
    atom = split_law[0]
    for q in range(1, min(b, h_max) + 1):
        atom += split_law[q] * Fraction(math.comb(n - b, h_max - q), den)
    if atom.numerator << cap_shift > atom.denominator:
        raise SystemExit(f"rational h<=500: termination atom exceeds 2^-{cap_shift}")
    return atom


@dataclass(frozen=True)
class ChernoffEnvelope:
    z: Fraction
    mgf: Fraction
    distance: int
    crossing: int
    dyadic_bits: int
    tail_ratio: Fraction
    anchor_scaled: int

    def raw_num_den(self, live: int) -> tuple[int, int]:
        # z^{-d} M(z)^live
        return (
            pow(self.mgf.numerator, live) * pow(self.z.denominator, self.distance),
            pow(self.mgf.denominator, live) * pow(self.z.numerator, self.distance),
        )

    def scaled_upper(self, live: int) -> int:
        scale = 1 << self.dyadic_bits
        if live <= self.crossing:
            return scale
        # The exact Chernoff value at crossing+1 is below one.  Since each
        # following live block multiplies it by M(z) < 2909/5000, the small-base
        # geometric envelope below avoids repeated exponentiation of the very
        # large exact MGF numerator and denominator.
        exponent = live - self.crossing - 1
        num = pow(self.tail_ratio.numerator, exponent)
        den = pow(self.tail_ratio.denominator, exponent)
        return min(scale, ceil_div(self.anchor_scaled * num, den))


def build_chernoff_envelope(
    mgf: Fraction, *, distance: int, dyadic_bits: int, crossing: int = 5923
) -> ChernoffEnvelope:
    z = Fraction(2333, 2373)
    # Freeze the precomputed candidate but prove both adjacent inequalities
    # exactly.  Reusing the powers at `crossing` avoids a costly big-integer
    # binary search in the default verifier.
    left = pow(mgf.numerator, crossing) * pow(z.denominator, distance)
    right = pow(mgf.denominator, crossing) * pow(z.numerator, distance)
    if left <= right or left * mgf.numerator > right * mgf.denominator:
        raise SystemExit("rational h<=500: exact Chernoff crossing check failed")
    tail_ratio = Fraction(2909, 5000)
    if mgf > tail_ratio:
        raise SystemExit("rational h<=500: 2909/5000 does not upper-bound the exact MGF")
    anchor_scaled = ceil_div(
        left * mgf.numerator << dyadic_bits,
        right * mgf.denominator,
    )
    if not 0 < anchor_scaled < (1 << dyadic_bits):
        raise SystemExit("rational h<=500: invalid outward Chernoff anchor")
    return ChernoffEnvelope(
        z, mgf, distance, crossing, dyadic_bits, tail_ratio, anchor_scaled
    )


def nonempty_coefficients(*, b: int, h_max: int) -> list[list[int]]:
    """Return c[x][H] = [z^H] ((1+z)^b-1)^x."""

    rows = [[0] * (h_max + 1) for _ in range(h_max + 1)]
    rows[0][0] = 1
    block = [math.comb(b, r) for r in range(b + 1)]
    for x in range(1, h_max + 1):
        prev = rows[x - 1]
        row = rows[x]
        for H in range(x, h_max + 1):
            row[H] = sum(prev[H - r] * block[r] for r in range(1, min(b, H) + 1))
    return rows


def adjacent_ratio_num_den(T: int, x: int, e: int, live: int) -> tuple[int, int]:
    return (T - live) * (live + 1 - e), (T - live + e - 1) * (live + 1 - x)


def peak_live(
    T: int,
    x: int,
    e: int,
    lo: int,
    hi: int,
    *,
    multiplier_num: int = 1,
    multiplier_den: int = 1,
) -> int:
    """Locate a maximum of a_l multiplier^l on an integer interval."""

    if lo >= hi:
        return lo

    def grows(live: int) -> bool:
        num, den = adjacent_ratio_num_den(T, x, e, live)
        return num * multiplier_num > den * multiplier_den

    left, right = lo, hi
    while left < right:
        mid = (left + right) // 2
        if grows(mid):
            left = mid + 1
        else:
            right = mid
    return left


def a_term(T: int, x: int, e: int, live: int) -> int:
    return math.comb(T - live + e - 1, e - 1) * math.comb(live - e, x - e)


def self_check_peak_live() -> None:
    """Exhaust small cases for both log-concave peak searches."""

    for T in range(2, 18):
        for x in range(1, T + 1):
            for e in range(1, min(8, x) + 1):
                for multiplier in (Fraction(1), Fraction(2909, 5000)):
                    peak = peak_live(
                        T,
                        x,
                        e,
                        x,
                        T,
                        multiplier_num=multiplier.numerator,
                        multiplier_den=multiplier.denominator,
                    )
                    peak_value = a_term(T, x, e, peak) * multiplier**peak
                    brute = max(a_term(T, x, e, live) * multiplier**live for live in range(x, T + 1))
                    if peak_value != brute:
                        raise SystemExit("rational h<=500: suffix peak self-check failed")


def episode_suffix_scaled(
    T: int, x: int, e: int, envelope: ChernoffEnvelope, survival_cache: dict[int, int]
) -> int:
    """Upper-bound the skipped-gap sum, scaled by 2^Q."""

    scale = 1 << envelope.dyadic_bits
    if e == x + 1:
        # All x+1 gaps are selected, so the only skipped count is T-x.
        return math.comb(T, x) * scale
    if not 1 <= e <= x:
        return 0

    total = 0
    low_hi = min(T, envelope.crossing)
    if x <= low_hi:
        peak = peak_live(T, x, e, x, low_hi)
        total += (low_hi - x + 1) * a_term(T, x, e, peak) * scale

    high_lo = max(x, envelope.crossing + 1)
    if high_lo <= T:
        peak = peak_live(
            T,
            x,
            e,
            high_lo,
            T,
            multiplier_num=envelope.tail_ratio.numerator,
            multiplier_den=envelope.tail_ratio.denominator,
        )
        survival = survival_cache.get(peak)
        if survival is None:
            survival = envelope.scaled_upper(peak)
            survival_cache[peak] = survival
        total += (T - high_lo + 1) * a_term(T, x, e, peak) * survival
    return total


def inner_bounds_scaled(
    T: int,
    *,
    b: int,
    h_max: int,
    e_max: int,
    termination_shift: int,
    envelope: ChernoffEnvelope,
    coeffs: list[list[int]],
) -> list[int]:
    """Outward dyadic bounds for the e=0,...,e_max inner probabilities."""

    scale = 1 << envelope.dyadic_bits
    survival_cache: dict[int, int] = {}
    suffix = [[0] * (h_max + 1) for _ in range(e_max + 1)]
    for e in range(1, e_max + 1):
        for x in range(e - 1, h_max + 1):
            suffix[e][x] = episode_suffix_scaled(
                T, x, e, envelope, survival_cache
            )
    out = [0] * (h_max + 1)
    e0 = envelope.scaled_upper(T)
    for H in range(h_max + 1):
        denom = math.comb(b * T, H)
        inner = e0
        for e in range(1, min(e_max, H + 1) + 1):
            numerator = 0
            for x in range(max(0, e - 1), min(T, H) + 1):
                coeff = coeffs[x][H]
                if coeff:
                    numerator += (
                        coeff
                        * math.comb(x + 1, e)
                        * suffix[e][x]
                    )
            inner += ceil_div(numerator, denom << (termination_shift * e))
        out[H] = min(scale, inner)
    return out


def generate_gap_sums(*, b: int, h_max: int) -> list[list[int]]:
    """Exact placement sums, one [H=0..h_max] row per bucket."""

    totals = [[0] * (h_max + 1) for _ in BUCKETS]
    for bucket_index, (gap_min, gap_max, _T) in enumerate(BUCKETS):
        row = totals[bucket_index]
        for gap in range(gap_min, gap_max + 1):
            n = b * (5949 + gap - 1)
            choose = 1
            row[0] += 1
            for H in range(1, h_max + 1):
                choose = choose * (n - H + 1) // H
                row[H] += choose
    return totals


def gap_payload(sums: list[list[int]], *, b: int, h_max: int) -> dict[str, object]:
    data = {
        "schema": 1,
        "block_bits": b,
        "h_max": h_max,
        "buckets": [
            {"gap_min": lo, "gap_max": hi, "T": T, "sums": [str(value) for value in row]}
            for (lo, hi, T), row in zip(BUCKETS, sums)
        ],
    }
    canonical = json.dumps(data, sort_keys=True, separators=(",", ":")).encode()
    data["sha256"] = hashlib.sha256(canonical).hexdigest()
    return data


def write_gap_sums(path: Path, *, b: int, h_max: int) -> None:
    payload = gap_payload(generate_gap_sums(b=b, h_max=h_max), b=b, h_max=h_max)
    path.write_text(json.dumps(payload, indent=2) + "\n")


def load_gap_sums(path: Path, *, b: int, h_max: int) -> tuple[list[list[int]], str]:
    payload = json.loads(path.read_text())
    claimed = payload.pop("sha256")
    canonical = json.dumps(payload, sort_keys=True, separators=(",", ":")).encode()
    actual = hashlib.sha256(canonical).hexdigest()
    if claimed != actual:
        raise SystemExit("rational h<=500: gap-sum artifact SHA-256 mismatch")
    if payload["schema"] != 1 or payload["block_bits"] != b or payload["h_max"] != h_max:
        raise SystemExit("rational h<=500: gap-sum artifact metadata mismatch")
    rows: list[list[int]] = []
    for expected, item in zip(BUCKETS, payload["buckets"], strict=True):
        if (item["gap_min"], item["gap_max"], item["T"]) != expected:
            raise SystemExit("rational h<=500: gap-sum bucket metadata mismatch")
        row = [int(value) for value in item["sums"]]
        if len(row) != h_max + 1 or row[0] != expected[1] - expected[0] + 1:
            raise SystemExit("rational h<=500: malformed gap-sum row")
        rows.append(row)
    if len(rows) != len(BUCKETS):
        raise SystemExit("rational h<=500: wrong gap-sum bucket count")
    return rows, actual


def inner_payload(
    rows: list[list[int]],
    *,
    envelope: ChernoffEnvelope,
    h_max: int,
    bucket_episode_max: tuple[int, ...],
) -> dict[str, object]:
    data = {
        "schema": 1,
        "h_max": h_max,
        "dyadic_bits": envelope.dyadic_bits,
        "crossing": envelope.crossing,
        "mgf_numerator": str(envelope.mgf.numerator),
        "mgf_denominator": str(envelope.mgf.denominator),
        "tail_ratio_numerator": envelope.tail_ratio.numerator,
        "tail_ratio_denominator": envelope.tail_ratio.denominator,
        "anchor_scaled": str(envelope.anchor_scaled),
        "buckets": [
            {
                "gap_min": lo,
                "gap_max": hi,
                "T": T,
                "episode_max": e_max,
                "inner_scaled": [str(value) for value in row],
            }
            for (lo, hi, T), e_max, row in zip(
                BUCKETS, bucket_episode_max, rows, strict=True
            )
        ],
    }
    canonical = json.dumps(data, sort_keys=True, separators=(",", ":")).encode()
    data["sha256"] = hashlib.sha256(canonical).hexdigest()
    return data


def write_inner_bounds(
    path: Path,
    *,
    envelope: ChernoffEnvelope,
    b: int,
    h_max: int,
    termination_shift: int,
    bucket_episode_max: tuple[int, ...],
) -> None:
    coeffs = nonempty_coefficients(b=b, h_max=h_max)
    rows = [
        inner_bounds_scaled(
            T,
            b=b,
            h_max=h_max,
            e_max=e_max,
            termination_shift=termination_shift,
            envelope=envelope,
            coeffs=coeffs,
        )
        for (_lo, _hi, T), e_max in zip(BUCKETS, bucket_episode_max, strict=True)
    ]
    path.write_text(
        json.dumps(
            inner_payload(
                rows,
                envelope=envelope,
                h_max=h_max,
                bucket_episode_max=bucket_episode_max,
            ),
            indent=2,
        )
        + "\n"
    )


def load_inner_bounds(
    path: Path,
    *,
    envelope: ChernoffEnvelope,
    h_max: int,
    bucket_episode_max: tuple[int, ...],
) -> tuple[list[list[int]], str]:
    payload = json.loads(path.read_text())
    claimed = payload.pop("sha256")
    canonical = json.dumps(payload, sort_keys=True, separators=(",", ":")).encode()
    actual = hashlib.sha256(canonical).hexdigest()
    if claimed != actual:
        raise SystemExit("rational h<=500: inner-bound artifact SHA-256 mismatch")
    expected_metadata = (
        payload["schema"] == 1
        and payload["h_max"] == h_max
        and payload["dyadic_bits"] == envelope.dyadic_bits
        and payload["crossing"] == envelope.crossing
        and int(payload["mgf_numerator"]) == envelope.mgf.numerator
        and int(payload["mgf_denominator"]) == envelope.mgf.denominator
        and payload["tail_ratio_numerator"] == envelope.tail_ratio.numerator
        and payload["tail_ratio_denominator"] == envelope.tail_ratio.denominator
        and int(payload["anchor_scaled"]) == envelope.anchor_scaled
    )
    if not expected_metadata:
        raise SystemExit("rational h<=500: inner-bound artifact metadata mismatch")
    rows: list[list[int]] = []
    for expected, e_max, item in zip(
        BUCKETS, bucket_episode_max, payload["buckets"], strict=True
    ):
        if (
            (item["gap_min"], item["gap_max"], item["T"]) != expected
            or item["episode_max"] != e_max
        ):
            raise SystemExit("rational h<=500: inner-bound bucket metadata mismatch")
        row = [int(value) for value in item["inner_scaled"]]
        if len(row) != h_max + 1:
            raise SystemExit("rational h<=500: malformed inner-bound row")
        rows.append(row)
    if len(rows) != len(BUCKETS):
        raise SystemExit("rational h<=500: wrong inner-bound bucket count")
    return rows, actual


@dataclass(frozen=True)
class RationalH500Certificate:
    h_max: int
    live_outer_weights: int
    chernoff_crossing: int
    chernoff_mgf_log2: float
    termination_exact_log2: float
    termination_cap_shift: int
    dyadic_bits: int
    episode_max: int
    early_episode_max: int
    gap_sums_sha256: str
    inner_bounds_sha256: str
    local_spectrum_sha256: str
    inner_spectrum_sha256: str
    bucket_log2: tuple[float, ...]
    explicit_episode_log2: float
    episode_tail_log2: float
    ultra_late_log2: float
    total_log2: float
    threshold_num: int
    threshold_shift: int
    threshold_log2: float
    threshold_bits_exact_status: str
    cross_multiply_slack_bits: int


def certify(
    *,
    inner_spectrum: Path,
    local_spectrum_csv: Path,
    gap_sums_path: Path,
    inner_bounds_path: Path = DEFAULT_INNER_BOUNDS,
    n: int = 2**21,
    b: int = 64,
    outer_blocks: int = 4096,
    late_blocks: int = 5949,
    h_max: int = 500,
    episode_max: int = 8,
    early_episode_max: int = 16,
    termination_shift: int = 63,
    dyadic_bits: int = 1024,
    threshold_num: int = 6793,
    threshold_shift: int = 50,
) -> RationalH500Certificate:
    if (
        h_max != 500
        or b != 64
        or late_blocks != 5949
        or episode_max != 8
        or early_episode_max != 16
    ):
        raise SystemExit(
            "rational h<=500: this frozen certificate expects h_max=500, "
            "b=64, late=5949, prefix e_max=8, early e_max=16"
        )
    self_check_peak_live()
    load_rm512_spectrum(local_spectrum_csv)
    load_ebch128_spectrum(inner_spectrum)
    outer = load_outer_coefficients(local_spectrum_csv, blocks=outer_blocks, h_max=h_max)
    split_law, mgf = exact_split_law(inner_spectrum, b)
    termination = certify_termination_cap(
        split_law, b=b, T_min=late_blocks, h_max=h_max - 1, cap_shift=termination_shift
    )
    envelope = build_chernoff_envelope(mgf, distance=(9 * n) // 100, dyadic_bits=dyadic_bits)
    gap_sums, artifact_sha = load_gap_sums(gap_sums_path, b=b, h_max=h_max - 1)

    bucket_episode_max = (episode_max,) * 3 + (early_episode_max,) * 3
    inner_by_bucket, inner_artifact_sha = load_inner_bounds(
        inner_bounds_path,
        envelope=envelope,
        h_max=h_max - 1,
        bucket_episode_max=bucket_episode_max,
    )
    scale = 1 << dyadic_bits
    bucket_totals = [Fraction(0) for _ in BUCKETS]
    tail_total = Fraction(0)
    total_gap_sums = [sum(row[H] for row in gap_sums) for H in range(h_max)]
    # The prefix explicitly retains e<=8; the early buckets retain e<=16.
    # Round their tiny exact tails upward once at a finer dyadic precision.
    # At Q=2048 the accumulated rounding cover is vastly below the tail.
    tail_bits = 2048
    tail_scale = 1 << tail_bits
    prefix_tail_scaled = []
    early_tail_scaled = []

    def tail_scaled(H: int, e_min: int) -> int:
        total = 0
        for e in range(e_min, H + 2):
            coefficient = math.comb(H + 1, e)
            shift = tail_bits - termination_shift * e
            total += coefficient << shift if shift >= 0 else ceil_div(coefficient, 1 << -shift)
        return total

    for H in range(h_max):
        prefix_tail_scaled.append(tail_scaled(H, episode_max + 1))
        early_tail_scaled.append(tail_scaled(H, early_episode_max + 1))
    live_weights = 0
    for h in range(1, h_max + 1):
        A_h = outer[h]
        if not A_h:
            continue
        live_weights += 1
        den_h = math.comb(n, h)
        for bucket_index, gap_row in enumerate(gap_sums):
            placement_inner_num = 0
            for r in range(1, min(b, h) + 1):
                H = h - r
                placement_inner_num += (
                    math.comb(b, r) * gap_row[H] * inner_by_bucket[bucket_index][H]
                )
            bucket_totals[bucket_index] += Fraction(A_h * placement_inner_num, den_h * scale)

        tail_inner_scaled = 0
        for r in range(1, min(b, h) + 1):
            H = h - r
            first = math.comb(b, r)
            prefix_placement = sum(gap_sums[index][H] for index in range(3))
            early_placement = total_gap_sums[H] - prefix_placement
            tail_inner_scaled += first * (
                prefix_placement * prefix_tail_scaled[H]
                + early_placement * early_tail_scaled[H]
            )
        tail_total += Fraction(A_h * tail_inner_scaled, den_h * tail_scale)

    episode_total = sum(bucket_totals, Fraction(0))
    ultra_late = sum(
        (
            Fraction(outer[h] * math.comb(b * late_blocks, h), math.comb(n, h))
            for h in range(1, h_max + 1)
            if outer[h]
        ),
        Fraction(0),
    )
    total = episode_total + tail_total + ultra_late
    threshold = Fraction(threshold_num, 1 << threshold_shift)
    if total > threshold:
        raise SystemExit(
            f"rational h<=500: total {log2_fraction(total):.12f} exceeds threshold "
            f"{threshold_num}/2^{threshold_shift}; buckets="
            f"{[log2_fraction(value) for value in bucket_totals]}, "
            f"tail={log2_fraction(tail_total):.12f}, ultra={log2_fraction(ultra_late):.12f}"
        )
    # Prove threshold <= 2^(-37.27) without transcendental arithmetic:
    # (threshold)^100 <= 2^-3727.
    if pow(threshold_num, 100) > 1 << (100 * threshold_shift - 3727):
        raise SystemExit("rational h<=500: threshold does not certify 37.27 bits")

    return RationalH500Certificate(
        h_max=h_max,
        live_outer_weights=live_weights,
        chernoff_crossing=envelope.crossing,
        chernoff_mgf_log2=log2_fraction(mgf),
        termination_exact_log2=log2_fraction(termination),
        termination_cap_shift=termination_shift,
        dyadic_bits=dyadic_bits,
        episode_max=episode_max,
        early_episode_max=early_episode_max,
        gap_sums_sha256=artifact_sha,
        inner_bounds_sha256=inner_artifact_sha,
        local_spectrum_sha256=file_sha256(local_spectrum_csv),
        inner_spectrum_sha256=file_sha256(inner_spectrum),
        bucket_log2=tuple(log2_fraction(value) for value in bucket_totals),
        explicit_episode_log2=log2_fraction(episode_total),
        episode_tail_log2=log2_fraction(tail_total),
        ultra_late_log2=log2_fraction(ultra_late),
        total_log2=log2_fraction(total),
        threshold_num=threshold_num,
        threshold_shift=threshold_shift,
        threshold_log2=log2_fraction(threshold),
        threshold_bits_exact_status="PASS",
        cross_multiply_slack_bits=(threshold.numerator * total.denominator - total.numerator * threshold.denominator).bit_length(),
    )


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--inner-spectrum", type=Path, default=ROOT / "EBCH128_64.wd")
    parser.add_argument("--local-spectrum-csv", type=Path, default=ROOT / "rm512_256_spectrum.csv")
    parser.add_argument("--gap-sums", type=Path, default=DEFAULT_GAP_SUMS)
    parser.add_argument("--inner-bounds", type=Path, default=DEFAULT_INNER_BOUNDS)
    parser.add_argument("--generate-gap-sums", action="store_true")
    parser.add_argument("--recompute-gap-sums", action="store_true")
    parser.add_argument("--generate-inner-bounds", action="store_true")
    parser.add_argument("--recompute-inner-bounds", action="store_true")
    args = parser.parse_args()
    if args.generate_gap_sums:
        write_gap_sums(args.gap_sums, b=64, h_max=499)
        print(f"gap_sums_written,{args.gap_sums}")
        return 0
    if args.recompute_gap_sums:
        stored, stored_sha = load_gap_sums(args.gap_sums, b=64, h_max=499)
        regenerated = generate_gap_sums(b=64, h_max=499)
        if stored != regenerated:
            raise SystemExit("rational h<=500: regenerated gap sums differ from artifact")
        print("gap_sums_recompute_status,PASS")
        print(f"gap_sums_sha256,{stored_sha}")
        return 0
    if args.generate_inner_bounds or args.recompute_inner_bounds:
        _split_law, mgf = exact_split_law(args.inner_spectrum, 64)
        envelope = build_chernoff_envelope(
            mgf, distance=(9 * (2**21)) // 100, dyadic_bits=1024
        )
        episode_limits = (8, 8, 8, 16, 16, 16)
        if args.recompute_inner_bounds:
            stored, stored_sha = load_inner_bounds(
                args.inner_bounds,
                envelope=envelope,
                h_max=499,
                bucket_episode_max=episode_limits,
            )
            coeffs = nonempty_coefficients(b=64, h_max=499)
            regenerated = [
                inner_bounds_scaled(
                    T,
                    b=64,
                    h_max=499,
                    e_max=e_max,
                    termination_shift=63,
                    envelope=envelope,
                    coeffs=coeffs,
                )
                for (_lo, _hi, T), e_max in zip(BUCKETS, episode_limits, strict=True)
            ]
            if stored != regenerated:
                raise SystemExit("rational h<=500: regenerated inner bounds differ from artifact")
            print("inner_bounds_recompute_status,PASS")
            print(f"inner_bounds_sha256,{stored_sha}")
            return 0
        write_inner_bounds(
            args.inner_bounds,
            envelope=envelope,
            b=64,
            h_max=499,
            termination_shift=63,
            bucket_episode_max=episode_limits,
        )
        print(f"inner_bounds_written,{args.inner_bounds}")
        return 0
    result = certify(
        inner_spectrum=args.inner_spectrum,
        local_spectrum_csv=args.local_spectrum_csv,
        gap_sums_path=args.gap_sums,
        inner_bounds_path=args.inner_bounds,
    )
    print("fullsplit_h500_rational_status,PASS")
    for key, value in asdict(result).items():
        if isinstance(value, tuple):
            for index, item in enumerate(value, 1):
                print(f"{key}_{index},{item:.12f}")
        elif isinstance(value, float):
            print(f"{key},{value:.12f}")
        else:
            print(f"{key},{value}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
