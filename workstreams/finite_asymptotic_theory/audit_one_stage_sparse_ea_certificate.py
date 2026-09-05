#!/usr/bin/env python3
"""Audit the complete one-stage sparse-EA finite certificate chain."""

from __future__ import annotations

import hashlib
import json
from pathlib import Path


WORKSTREAM = Path(__file__).resolve().parent
PURE = WORKSTREAM / "pure_expander_accumulate"
OUTPUT = WORKSTREAM / "ONE_STAGE_SPARSE_EA_CERTIFICATE_MANIFEST.json"

FILES = {
    "sector_zero": PURE / "sector_zero_complement_orbits_K256_B512_r33_p1_159_rust20_outward.json",
    "sectors_one_two": PURE / "primal_schur_sectors1_2_K256_B512_r33_p1_159_rust20_outward.json",
    "sector_two_mass": PURE / "primal_schur_sector2_K256_B512_r33_p2_159_mass_rust20_outward.json",
    "defect_low": PURE / "dominant_character_likelihood_K256_B512_r33_w42_79_outward.json",
    "central": PURE / "dominant_character_deviation_K256_B512_r33_w80_432_outward.json",
    "defect_high": PURE / "dominant_character_likelihood_K256_B512_r33_w433_470_outward.json",
    "small_exact": PURE / "dominant_character_deviation_K4_B8_r3_all_shells_exact.json",
    "parallel_validation": PURE / "dominant_character_deviation_peach_workers8_validation.json",
    "caps": WORKSTREAM / "one_stage_sparse_ea_K256_B512_r33_spectrum_caps_outward.json",
    "transfer": WORKSTREAM / "one_stage_sparse_ea_K256_B512_r33_randomstepconv_M22_d109_outward.json",
}


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def semantic_shell_claim_sha256(payload: dict[str, object]) -> str:
    selected = {"shells": payload["shells"], "claim": payload["claim"]}
    encoded = json.dumps(
        selected, sort_keys=True, separators=(",", ":")
    ).encode("utf-8")
    return hashlib.sha256(encoded).hexdigest()


def load(name: str) -> dict[str, object]:
    return json.loads(FILES[name].read_text(encoding="utf-8"))


def check_embedded_sources(payload: dict[str, object], owner: Path) -> None:
    for row in payload.get("sources", []):  # type: ignore[union-attr]
        path = owner.parent / row["file"]
        if not path.exists():
            # The cap receipt names sources in its pure subdirectory.
            path = PURE / row["file"]
        if sha256(path) != row["sha256"]:
            raise AssertionError(f"stale embedded source hash: {path}")


def main() -> None:
    sector_zero = load("sector_zero")
    if sector_zero.get("status") != "OUTWARD_BINARY64_UPPER_BOUND":
        raise AssertionError("sector-zero receipt is not accepting")
    sector_zero_pairs = [
        (int(row["sector"]), int(row["level"]))
        for row in sector_zero["entries"]  # type: ignore[index]
    ]
    if sector_zero_pairs != [(0, level) for level in range(1, 160)]:
        raise AssertionError("sector-zero levels are incomplete")

    sectors = load("sectors_one_two")
    if sectors.get("status") != "OUTWARD_BINARY64_UPPER_BOUND":
        raise AssertionError("sector-one/two receipt is not accepting")
    sector_pairs = {
        (int(row["sector"]), int(row["level"]))
        for row in sectors["entries"]  # type: ignore[index]
    }
    expected_pairs = {
        *((1, level) for level in range(1, 160)),
        *((2, level) for level in range(2, 160)),
    }
    if sector_pairs != expected_pairs:
        raise AssertionError("sector-one/two levels are incomplete")

    sector_two = load("sector_two_mass")
    if sector_two.get("status") != "OUTWARD_BINARY64_UPPER_BOUND":
        raise AssertionError("sector-two mass receipt is not accepting")
    sector_two_rows = {
        int(row["level"]): float(row["scaled_diagonal_upper"])
        for row in sector_two["entries"]  # type: ignore[index]
    }
    combined_sector_two = {
        int(row["level"]): float(row["scaled_diagonal_upper"])
        for row in sectors["entries"]  # type: ignore[index]
        if int(row["sector"]) == 2
    }
    if sector_two_rows != combined_sector_two:
        raise AssertionError("sector-two receipts disagree")

    shell_sources = {
        "defect_low": range(42, 80),
        "central": range(80, 433),
        "defect_high": range(433, 471),
    }
    for name, shell_range in shell_sources.items():
        payload = load(name)
        if payload.get("status") != "OUTWARD_ARB_AND_BINARY64_CERTIFICATE":
            raise AssertionError(f"{name} receipt is not accepting")
        shells = [int(row["shell_weight"]) for row in payload["shells"]]  # type: ignore[index]
        if shells != list(shell_range):
            raise AssertionError(f"{name} shell coverage is incomplete")
    for name in ("defect_low", "defect_high"):
        if not load(name)["claim"]["all_requested_shells_pass_factor_512"]:  # type: ignore[index]
            raise AssertionError(f"{name} misses the factor-512 gate")

    small = load("small_exact")
    if small.get("status") != "EXACT_RATIONAL_SMALL_MODEL":
        raise AssertionError("small exact receipt has the wrong status")
    if not small["verified"]["all_shell_aggregates_dominated"]:  # type: ignore[index]
        raise AssertionError("small exact deviation audit failed")

    parallel = load("parallel_validation")
    if parallel.get("status") != "EXACT_SEMANTIC_FIELD_COMPARISON":
        raise AssertionError("parallel validation has the wrong status")
    if not parallel["comparison"]["shells_equal"]:  # type: ignore[index]
        raise AssertionError("parallel and serial shell rows differ")
    if not parallel["comparison"]["claim_equal"]:  # type: ignore[index]
        raise AssertionError("parallel and serial claims differ")
    central = load("central")
    if sha256(FILES["central"]) != parallel["sha256"]["current_local_serial_receipt"]:  # type: ignore[index]
        raise AssertionError("parallel validation names a stale local serial receipt")
    if semantic_shell_claim_sha256(central) != parallel["sha256"]["semantic_shells_and_claim"]:  # type: ignore[index]
        raise AssertionError("current local serial semantic fields differ from Peach validation")

    caps = load("caps")
    if caps.get("status") != "OUTWARD_CAP_CERTIFICATE":
        raise AssertionError("cap receipt is not accepting")
    if not all(caps["band_majorant_dominated_by_frozen_random"].values()):  # type: ignore[index]
        raise AssertionError("a cap band exceeds the frozen transfer envelope")
    if float(caps["claim"]["total_setup_margin_bits_lower"]) <= 40:  # type: ignore[index]
        raise AssertionError("setup event misses 40 bits")
    check_embedded_sources(caps, FILES["caps"])

    transfer = load("transfer")
    if transfer.get("status") != "OUTWARD_FINITE_DISTANCE_CERTIFICATE":
        raise AssertionError("finite transfer is not accepting")
    if not transfer["claim"]["comparison_to_2^-40"]:  # type: ignore[index]
        raise AssertionError("finite transfer misses 40 bits")
    if float(transfer["claim"]["overall_margin_bits_lower"]) <= 40:  # type: ignore[index]
        raise AssertionError("finite transfer margin is not above 40 bits")
    check_embedded_sources(transfer, FILES["transfer"])

    payload = {
        "schema": "one-stage-sparse-ea-certificate-manifest-v1",
        "status": "COMPLETE_CERTIFICATE_CHAIN_AUDIT",
        "claim": {
            "sector_zero_levels": [1, 159],
            "sector_one_levels": [1, 159],
            "sector_two_levels": [2, 159],
            "defect_shells": [[42, 79], [433, 470]],
            "central_shells": [80, 432],
            "cap_support": [42, 470],
            "bounded_rank_test_attempts": 16,
            "distance_cutoff": int(transfer["claim"]["minimum_distance_at_least"]),  # type: ignore[index]
            "relative_distance_lower": float(transfer["claim"]["relative_distance_lower"]),  # type: ignore[index]
            "overall_margin_bits_lower": float(transfer["claim"]["overall_margin_bits_lower"]),  # type: ignore[index]
        },
        "checks": [
            "complete outward sector-zero level range",
            "complete outward controlling sector-one and sector-two ranges",
            "independent sector-two receipt equality",
            "complete low, central, and high shell ranges",
            "factor-512 defect-shell gates",
            "exact small-model likelihood and covariance audit",
            "exact serial-versus-parallel shell and claim equality, tied to the current local receipt by a semantic hash",
            "cap-source hashes and frozen-band domination",
            "transfer-source hashes and end-to-end 40-bit gate",
        ],
        "artifacts": [
            {"role": name, "file": str(path.relative_to(WORKSTREAM)), "sha256": sha256(path)}
            for name, path in FILES.items()
        ],
    }
    OUTPUT.write_text(json.dumps(payload, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(payload["claim"], indent=2))
    print(f"output={OUTPUT}")


if __name__ == "__main__":
    main()
