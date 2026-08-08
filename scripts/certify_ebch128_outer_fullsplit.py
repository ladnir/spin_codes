#!/usr/bin/env python3
"""Certify the full-split construction with an EBCH [128,64,22] block outer.

This is a separate theorem lane from the frozen RM(4,9) [512,256,32]
certificate.  It reuses the exact low-weight artifacts and the rational/
outward post-prefix lemmas, with outer-block metadata and rational poles tuned
for 16384 copies of the committed extended-BCH spectrum.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import math
from dataclasses import asdict, dataclass
from decimal import Decimal
from fractions import Fraction
from pathlib import Path

from certificate_spectra import file_sha256, load_csv_spectrum, load_ebch128_spectrum
from certify_fullsplit_h500_rational import (
    DEFAULT_GAP_SUMS,
    DEFAULT_INNER_BOUNDS,
    certify as certify_h500,
)
import certify_fullsplit_postprefix_rational as post


ROOT = Path(__file__).resolve().parent
DEFAULT_ARTIFACT = ROOT / "ebch128_outer_fullsplit_rational.json"
LOCAL_SPECTRUM = ROOT / "ebch128_64_spectrum.csv"
INNER_SPECTRUM = ROOT / "EBCH128_64.wd"
LOCAL_NAME = "extended BCH [128,64,22] outer spectrum"
OUTER_BLOCKS = 16384

LOCAL_KWARGS = {
    "outer_blocks": OUTER_BLOCKS,
    "local_name": LOCAL_NAME,
    "local_length": 128,
    "local_dimension": 64,
    "local_distance": 22,
}

CRITICAL_POLE_SCHEDULE = (
    (650, Fraction(997095156, 10**9), Fraction(1509346, 10**9)),
    (850, Fraction(996039186, 10**9), Fraction(2025006, 10**9)),
    (1050, Fraction(994984292, 10**9), Fraction(2546350, 10**9)),
    (1250, Fraction(993930504, 10**9), Fraction(3070553, 10**9)),
    (1450, Fraction(992877832, 10**9), Fraction(3596385, 10**9)),
    (1650, Fraction(991826276, 10**9), Fraction(4123229, 10**9)),
    (1850, Fraction(990775845, 10**9), Fraction(4650736, 10**9)),
    (2000, Fraction(989726533, 10**9), Fraction(5178708, 10**9)),
)

LATE_SEGMENTS = (
    (501, 2000, "0.348"),
    (2001, 7858, "0.36"),
    (7859, 20550, "0.375"),
    (20551, 75000, "0.39"),
    (75001, 150000, "0.4"),
    (150001, 250000, "0.407"),
    (250001, 380736, "0.415"),
)

EARLY_ENDPOINT_INTERVALS = (
    (2001, 7858, Fraction(4900255817, 5 * 10**9), Fraction(1, 100), "0.355"),
    (7859, 20550, Fraction(4707442269, 5 * 10**9), Fraction(3, 100), "0.375"),
    (20551, 30000, Fraction(8178969571, 10**10), Fraction(1, 10), "0.396059859435"),
    (30001, 50000, Fraction(8178969571, 10**10), Fraction(1, 10), "0.396059859435"),
    (50001, 75000, Fraction(8178969571, 10**10), Fraction(1, 10), "0.396059859435"),
    (75001, 90000, Fraction(8178969571, 10**10), Fraction(1, 10), "0.396059859435"),
    (90001, 100000, Fraction(8178969571, 10**10), Fraction(1, 10), "0.396059859435"),
    (100001, 110000, Fraction(5382199631, 10**10), Fraction(3, 10), "0.396059859435"),
    (110001, 125000, Fraction(5382199631, 10**10), Fraction(3, 10), "0.396059859435"),
    (125001, 160000, Fraction(5382199631, 10**10), Fraction(3, 10), "0.396059859435"),
    (160001, 250000, Fraction(3331239681, 10**10), Fraction(1, 2), "0.396059859435"),
    (250001, 350000, Fraction(3331239681, 10**10), Fraction(1, 2), "0.430812727615"),
)

PREFIX_ENDPOINT_INTERVALS = (
    (2001, 7858, Fraction(4900255817, 5 * 10**9), Fraction(1, 100), "0.355"),
    (7859, 20550, Fraction(4707442269, 5 * 10**9), Fraction(3, 100), "0.375"),
    (20551, 75000, Fraction(8178969571, 10**10), Fraction(1, 10), "0.396059859435"),
    (75001, 250000, Fraction(5382199631, 10**10), Fraction(3, 10), "0.396059859435"),
    (
        250001,
        350000,
        Fraction(3204075791888, 10**13),
        Fraction(7079866209348, 10**13),
        "0.4126497310152",
    ),
)

PREFIX_CAP_INTERVALS = (
    (2001, 7858, "0.355"),
    (7859, 20550, "0.375"),
    (20551, 75000, "0.39"),
    (75001, 250000, "0.405"),
    (250001, 350000, "0.415"),
    (350001, 400000, "0.42"),
    (400001, 450000, "0.425"),
    (450001, 550000, "0.43"),
    (550001, 636736, "0.44"),
)


@dataclass(frozen=True)
class LateAggregateCertificate:
    segments: tuple[post.LatePlacementCertificate, ...]
    total_log2_upper: float
    threshold_shift: int
    exact_union_gate: str
    status: str


@dataclass(frozen=True)
class CompleteCertificate:
    low_weight_log2_upper: float
    postprefix_log2_upper: float
    diagnostic_total_log2_upper: float
    exact_threshold: str
    certified_bits_hundredths: int
    status: str


def fraction_text(value: Fraction) -> str:
    return f"{value.numerator}/{value.denominator}"


def build_late_aggregate(local_spectrum_csv: Path) -> LateAggregateCertificate:
    segments = tuple(
        post.certify_late_placement(
            local_spectrum_csv=local_spectrum_csv,
            threshold_shift=405,
            z=Fraction(Decimal(z_text)),
            intervals=((h_lo, h_hi),),
            **LOCAL_KWARGS,
        )
        for h_lo, h_hi, z_text in LATE_SEGMENTS
    )
    total = max(item.total_log2_upper for item in segments) + math.log2(len(segments)) + 1e-9
    # Each of seven segments is below 2^-405, hence their union is below
    # 7*2^-405 < 2^-402 < 2^-400.  The looser stored gate is integer-exact.
    if total > -400:
        raise SystemExit("EBCH128 outer: late aggregate exceeds 2^-400")
    return LateAggregateCertificate(
        segments=segments,
        total_log2_upper=total,
        threshold_shift=400,
        exact_union_gate="7 * 2^-405 < 2^-400",
        status="PASS",
    )


def recompute(
    *, local_spectrum_csv: Path = LOCAL_SPECTRUM, inner_spectrum: Path = INNER_SPECTRUM
) -> dict[str, object]:
    post.self_check()
    load_csv_spectrum(
        local_spectrum_csv,
        name=LOCAL_NAME,
        length=128,
        dimension=64,
        minimum_distance=22,
    )
    load_ebch128_spectrum(inner_spectrum)
    q_law, entries = post.exact_split_data(inner_spectrum)
    turnoff = post.certify_global_turnoff(q_law)
    monotonic_rows = post.certify_critical_t_monotonicity(
        entries=entries,
        poles=tuple(pole for _h, pole, _rho in CRITICAL_POLE_SCHEDULE),
    )
    monotonic_rows += post.certify_endpoint_t_monotonicity(
        intervals=EARLY_ENDPOINT_INTERVALS,
        gaps=post.EARLY_GAPS,
        entries=entries,
    )
    monotonic_rows += post.certify_endpoint_t_monotonicity(
        intervals=PREFIX_ENDPOINT_INTERVALS,
        gaps=post.PREFIX_EPISODE_GAPS,
        entries=entries,
    )

    low = certify_h500(
        inner_spectrum=inner_spectrum,
        local_spectrum_csv=local_spectrum_csv,
        gap_sums_path=DEFAULT_GAP_SUMS,
        inner_bounds_path=DEFAULT_INNER_BOUNDS,
        threshold_num=19,
        threshold_shift=24,
        exact_target_bits_hundredths=1975,
        **LOCAL_KWARGS,
    )
    critical = post.certify_critical_window(
        local_spectrum_csv=local_spectrum_csv,
        inner_spectrum=inner_spectrum,
        threshold_shift=50,
        pole_schedule=CRITICAL_POLE_SCHEDULE,
        **LOCAL_KWARGS,
    )
    late = build_late_aggregate(local_spectrum_csv)
    complement = post.certify_complement_high(
        local_spectrum_csv=local_spectrum_csv,
        inner_spectrum=inner_spectrum,
        threshold_shift=30000,
        **LOCAL_KWARGS,
    )
    early_endpoint = post.certify_early_endpoint(
        local_spectrum_csv=local_spectrum_csv,
        inner_spectrum=inner_spectrum,
        threshold_shift=380,
        intervals=EARLY_ENDPOINT_INTERVALS,
        gaps=post.EARLY_GAPS,
        **LOCAL_KWARGS,
    )
    early_paired = post.certify_early_paired(
        local_spectrum_csv=local_spectrum_csv,
        inner_spectrum=inner_spectrum,
        threshold_shift=80000,
        **LOCAL_KWARGS,
    )
    prefix_cap = post.certify_prefix_cap(
        local_spectrum_csv=local_spectrum_csv,
        threshold_shift=320,
        intervals=PREFIX_CAP_INTERVALS,
        **LOCAL_KWARGS,
    )
    prefix_endpoint = post.certify_early_endpoint(
        local_spectrum_csv=local_spectrum_csv,
        inner_spectrum=inner_spectrum,
        threshold_shift=100,
        intervals=PREFIX_ENDPOINT_INTERVALS,
        gaps=post.PREFIX_EPISODE_GAPS,
        **LOCAL_KWARGS,
    )
    prefix_paired_rho1 = post.certify_early_paired(
        local_spectrum_csv=local_spectrum_csv,
        inner_spectrum=inner_spectrum,
        threshold_shift=130000,
        intervals=post.PREFIX_PAIRED_RHO1_INTERVALS,
        pole=Fraction(3011942119, 10**10),
        rho=Fraction(1, 1),
        T_min=9949,
        T_max=17948,
        **LOCAL_KWARGS,
    )
    prefix_paired_rho2 = post.certify_early_paired(
        local_spectrum_csv=local_spectrum_csv,
        inner_spectrum=inner_spectrum,
        threshold_shift=250000,
        intervals=post.PREFIX_PAIRED_RHO2_INTERVALS,
        pole=Fraction(3011942119, 10**10),
        rho=Fraction(2, 1),
        T_min=9949,
        T_max=17948,
        **LOCAL_KWARGS,
    )
    prefix_paired_rho3 = post.certify_early_paired(
        local_spectrum_csv=local_spectrum_csv,
        inner_spectrum=inner_spectrum,
        threshold_shift=280000,
        intervals=post.PREFIX_PAIRED_RHO3_INTERVALS,
        pole=Fraction(3011942119, 10**10),
        rho=Fraction(3, 1),
        T_min=9949,
        T_max=17948,
        **LOCAL_KWARGS,
    )

    post_components = [
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
    post_complete = post.combine_complete_postprefix(post_components, threshold_shift=50)

    full_threshold = Fraction(19, 1 << 24) + Fraction(1, 1 << 50)
    if (full_threshold.numerator**100) << 1975 > full_threshold.denominator**100:
        raise SystemExit("EBCH128 outer: combined threshold does not certify 19.75 bits")
    diagnostic_total = max(low.total_log2, post_complete.total_log2_upper)
    diagnostic_total += math.log2(
        1 + 2 ** (-abs(low.total_log2 - post_complete.total_log2_upper))
    )
    complete = CompleteCertificate(
        low_weight_log2_upper=low.total_log2,
        postprefix_log2_upper=post_complete.total_log2_upper,
        diagnostic_total_log2_upper=diagnostic_total,
        exact_threshold=f"{full_threshold.numerator}/{full_threshold.denominator}",
        certified_bits_hundredths=1975,
        status="PASS",
    )

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
        "construction": {
            "N": 2**21,
            "K": 2**20,
            "outer": "direct sum of 16384 extended BCH [128,64,22] blocks",
            "inner": "full-split extended BCH [128,64,22]",
            "target_distance": 188743,
        },
        "inputs": {
            "local_spectrum": {
                "path": local_spectrum_csv.name,
                "sha256": file_sha256(local_spectrum_csv),
                "structural_validation": "extended BCH [128,64,22] PASS",
            },
            "inner_spectrum": {
                "path": inner_spectrum.name,
                "sha256": file_sha256(inner_spectrum),
                "structural_validation": "extended BCH [128,64,22] PASS",
            },
            "gap_sums_sha256": file_sha256(DEFAULT_GAP_SUMS),
            "inner_bounds_sha256": file_sha256(DEFAULT_INNER_BOUNDS),
        },
        "pole_data": {
            "critical_schedule": [
                [h, fraction_text(pole), fraction_text(rho)]
                for h, pole, rho in CRITICAL_POLE_SCHEDULE
            ],
            "late_segments": [list(row) for row in LATE_SEGMENTS],
            "early_endpoint_intervals": [
                [lo, hi, fraction_text(pole), fraction_text(rho), z]
                for lo, hi, pole, rho, z in EARLY_ENDPOINT_INTERVALS
            ],
            "prefix_endpoint_intervals": [
                [lo, hi, fraction_text(pole), fraction_text(rho), z]
                for lo, hi, pole, rho, z in PREFIX_ENDPOINT_INTERVALS
            ],
            "prefix_cap_intervals": [list(row) for row in PREFIX_CAP_INTERVALS],
        },
        "exact_global_turnoff": {
            "status": "PASS",
            "threshold": "2^-62",
            "log2_for_display": float(post.log2_fraction(turnoff).hi),
        },
        "exact_t_monotonicity": {"status": "PASS", "rows": monotonic_rows},
        "low_weight": asdict(low),
        "postprefix_components": {
            name: asdict(item) for name, item in zip(names, post_components, strict=True)
        },
        "postprefix_complete": asdict(post_complete),
        "complete": asdict(complete),
    }


def artifact_payload(data: dict[str, object]) -> dict[str, object]:
    normalized = json.loads(json.dumps(data))
    canonical = json.dumps(normalized, sort_keys=True, separators=(",", ":")).encode()
    return {**normalized, "sha256": hashlib.sha256(canonical).hexdigest()}


def write_artifact(path: Path, data: dict[str, object]) -> None:
    path.write_text(json.dumps(artifact_payload(data), indent=2) + "\n")


def load_artifact(path: Path) -> dict[str, object]:
    payload = json.loads(path.read_text())
    claimed = payload.pop("sha256", None)
    canonical = json.dumps(payload, sort_keys=True, separators=(",", ":")).encode()
    actual = hashlib.sha256(canonical).hexdigest()
    if claimed != actual:
        raise SystemExit("EBCH128 outer: artifact SHA-256 mismatch")
    if payload.get("schema") != 1:
        raise SystemExit("EBCH128 outer: artifact schema mismatch")
    inputs = payload.get("inputs", {})
    if inputs.get("local_spectrum", {}).get("sha256") != file_sha256(LOCAL_SPECTRUM):
        raise SystemExit("EBCH128 outer: local spectrum fingerprint mismatch")
    if inputs.get("inner_spectrum", {}).get("sha256") != file_sha256(INNER_SPECTRUM):
        raise SystemExit("EBCH128 outer: inner spectrum fingerprint mismatch")
    if inputs.get("gap_sums_sha256") != file_sha256(DEFAULT_GAP_SUMS):
        raise SystemExit("EBCH128 outer: gap-sums artifact fingerprint mismatch")
    if inputs.get("inner_bounds_sha256") != file_sha256(DEFAULT_INNER_BOUNDS):
        raise SystemExit("EBCH128 outer: inner-bounds artifact fingerprint mismatch")
    load_csv_spectrum(
        LOCAL_SPECTRUM,
        name=LOCAL_NAME,
        length=128,
        dimension=64,
        minimum_distance=22,
    )
    load_ebch128_spectrum(INNER_SPECTRUM)
    return {**payload, "sha256": claimed}


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--artifact", type=Path, default=DEFAULT_ARTIFACT)
    parser.add_argument("--write-artifact", action="store_true")
    parser.add_argument("--recompute-artifact", action="store_true")
    args = parser.parse_args()
    if args.write_artifact and args.recompute_artifact:
        raise SystemExit("choose at most one artifact mode")
    if args.write_artifact:
        data = recompute()
        write_artifact(args.artifact, data)
        payload = load_artifact(args.artifact)
    elif args.recompute_artifact:
        stored = load_artifact(args.artifact)
        regenerated = artifact_payload(recompute())
        if stored != regenerated:
            raise SystemExit("EBCH128 outer: recomputed artifact differs from stored artifact")
        payload = stored
    else:
        payload = load_artifact(args.artifact)
    complete = payload["complete"]
    print("ebch128_outer_fullsplit_status,PASS")
    print(f"diagnostic_total_log2_upper,{complete['diagnostic_total_log2_upper']:.12f}")
    print(f"theorem_safe_threshold,{complete['exact_threshold']}")
    print(f"certified_bits,{complete['certified_bits_hundredths'] / 100:.2f}")
    print(f"artifact_sha256,{payload['sha256']}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
