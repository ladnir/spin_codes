#!/usr/bin/env python3
"""Audit spectra, coverage, and summary values for the rate-half Q1 sweep."""

from __future__ import annotations

import hashlib
import json
import math
import sys
from collections import defaultdict
from pathlib import Path


HERE = Path(__file__).resolve().parent
WORKSTREAM = HERE.parent
REPOSITORY = WORKSTREAM.parent.parent
sys.path.insert(0, str(WORKSTREAM))

from small_k_replay.evaluate_exact_spectra_q1_phase import (  # noqa: E402
    CONSTITUENTS,
    load_spectrum,
)


RECEIPT = HERE / "rate_half_family_k_margin_d100.json"
EXPECTED_MAXIMA = {
    "extended BCH [8,4,4] exact": -1.281736,
    "extended BCH [32,16,8] exact": 4.106765,
    "Philips shortened-XBCH [64,32,12] exact": 11.915619,
    "extended BCH [128,64,22] exact": 30.013770,
    "RM(1,3) [8,4,4] exact": -1.281736,
    "RM(2,5) [32,16,8] exact": 4.106765,
    "RM(3,7) [128,64,16] exact": 13.587380,
    "RM(4,9) [512,256,32] exact": 38.925445,
    "random full-rank [8,4] reused": -5.566247,
    "random full-rank [16,8] reused": -3.415112,
    "random full-rank [32,16] reused": 0.866304,
    "random full-rank [64,32] reused": 8.905753,
    "random full-rank [128,64] reused": 24.885346,
    "random full-rank [256,128] reused": 58.745221,
    "random full-rank [512,256] reused": 129.461348,
    "random full-rank [1024,512] reused": 274.879395,
}


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def main() -> None:
    for key in ("ebch8", "ebch32", "xbch64", "ebch128", "rm13", "rm25", "rm37", "rm49"):
        load_spectrum(CONSTITUENTS[key])

    manifest = json.loads(
        (REPOSITORY / "scripts" / "xbch64_32_philips_manifest.json").read_text()
    )
    xbch_path = CONSTITUENTS["xbch64"].spectrum_path
    assert sha256(xbch_path) == manifest["spectrum"]["sha256"]

    payload = json.loads(RECEIPT.read_text(encoding="utf-8"))
    assert payload["status"] == "BINARY64_DIAGNOSTIC"
    assert payload["coverage"]["bch_derived_block_bits"] == [8, 32, 64, 128]
    assert payload["coverage"]["rm_block_bits"] == [8, 32, 128, 512]
    assert payload["coverage"]["random_block_bits"] == [8, 16, 32, 64, 128, 256, 512, 1024]

    grouped = defaultdict(list)
    for row in payload["cases"]:
        assert row["memory_bits"] == row["message_exponent"] + 2
        assert row["distance_cutoff"] == math.ceil(0.10 * row["output_bits"])
        grouped[row["series"]].append(row)
    assert set(grouped) == set(EXPECTED_MAXIMA)
    for series, expected in EXPECTED_MAXIMA.items():
        observed = max(row["q1_margin_bits_diagnostic"] for row in grouped[series])
        assert abs(observed - expected) < 5e-7, (series, observed, expected)

    common = json.loads(
        (HERE / "matched_constituent_k_margin_d100.json").read_text(encoding="utf-8")
    )["cases"]
    expanded = {
        (row["series"], row["message_exponent"]): row["q1_margin_bits_diagnostic"]
        for row in payload["cases"]
    }
    for row in common:
        key = row["series"], row["message_exponent"]
        assert expanded[key] == row["q1_margin_bits_diagnostic"]

    print("PASS: spectra, manifest hash, family coverage, schedule, maxima, and baseline reproduction")


if __name__ == "__main__":
    main()
