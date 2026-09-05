#!/usr/bin/env python3
"""Audit the four endpoint-aware dimension-33 frontier certificates."""

from __future__ import annotations

import hashlib
import json
from pathlib import Path

from analyze_systematic_group_kernel import systematic_state_columns
from analyze_riffle_dp_g4_autonomous_map import (
    accumulate,
    apply_polynomial,
    build_apply,
    kernel_basis,
)
from audit_riffle_dp_2lap_g4_component_endpoint_dimension30 import (
    COMPONENT_DEGREES,
    audit_receipt,
    code_columns,
    observations,
    support_annihilator,
    support_degrees,
    support_dimension,
)
from probe_riffle_dp_g4_component_mixing import transpose_columns


ROOT = Path(__file__).resolve().parents[1]
EXPLORATIONS = ROOT / "explorations"
SOURCE = ROOT / "scripts" / "certify_riffle_dp_2lap_g4_component_24node.cpp"
RUNNER = ROOT / "scripts" / "run_riffle_dp_2lap_g4_component_dimension33.py"
BASE_RUNNER = ROOT / "scripts" / "run_riffle_dp_2lap_g4_component_dimension32.py"
BASE_AUDIT = EXPLORATIONS / "riffle_dp_2lap_g4_component_dimension32_audit.json"
PRIMARY_PROBE = EXPLORATIONS / "riffle_dp_2lap_g4_component_dimension33_primary_probe.json"
PRIMARY_RUN = EXPLORATIONS / "riffle_dp_2lap_g4_component_dimension33_primary_run.json"
INDEPENDENT_RUN = EXPLORATIONS / "riffle_dp_2lap_g4_component_dimension33_independent_run.json"
OUTPUT = EXPLORATIONS / "riffle_dp_2lap_g4_component_dimension33_audit.json"
FRONTIER = (
    (0x2E, (2, 4, 9, 18)),
    (0x35, (1, 4, 10, 18)),
    (0x4C, (4, 9, 20)),
    (0x53, (1, 2, 10, 20)),
)
DIMENSION = 33
EXPECTED_VECTORS = 12_216_115_184


def digest(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def executable(mode: str, window: int) -> Path:
    stem = "certify" if mode == "primary" else "verify"
    return ROOT / "scripts" / (
        f"{stem}_riffle_dp_2lap_g4_component_w{window}_d{DIMENSION}.exe"
    )


def validate_run(run: dict, mode: str, cases: int, probe_only: bool) -> None:
    if (
        run["schema"] != "riffle-dp-2lap-g4-component-dimension33-run-v1"
        or run["mode"] != mode
        or run["probe_only"] != probe_only
        or run["result"] != "PASS"
        or not run["completed"]
        or run["expected_case_count"] != cases
        or run["completed_case_count"] != cases
        or len(run["rows"]) != cases
    ):
        raise RuntimeError("dimension-33 audit: aggregate run changed")
    if run["runner_source_sha256"] != digest(RUNNER):
        raise RuntimeError("dimension-33 audit: runner hash mismatch")
    if run["base_runner_source_sha256"] != digest(BASE_RUNNER):
        raise RuntimeError("dimension-33 audit: base runner hash mismatch")
    if run["certificate_source_sha256"] != digest(SOURCE):
        raise RuntimeError("dimension-33 audit: certificate source hash mismatch")
    for key, expected_hash in run["executable_sha256_by_window_and_dimension"].items():
        window, dimension = key.split("_d")
        binary = executable(mode, int(window[1:]))
        if dimension != str(DIMENSION) or digest(binary) != expected_hash:
            raise RuntimeError("dimension-33 audit: executable hash mismatch")


def main() -> None:
    target32 = {
        mask
        for mask in range(1, 1 << len(COMPONENT_DEGREES))
        if support_dimension(mask) <= 32
    }
    target33 = {
        mask
        for mask in range(1, 1 << len(COMPONENT_DEGREES))
        if support_dimension(mask) <= 33
    }
    frontier_masks = {mask for mask, _ in FRONTIER}
    if len(target32) != 65 or len(target33) != 69:
        raise RuntimeError("dimension-33 audit: target support count changed")
    if target33 - target32 != frontier_masks:
        raise RuntimeError("dimension-33 audit: frontier changed")

    base = json.loads(BASE_AUDIT.read_text())
    if base["result"] != "PASS" or base["target_support_count"] != 65:
        raise RuntimeError("dimension-33 audit: base certificate changed")
    probe = json.loads(PRIMARY_PROBE.read_text())
    primary = json.loads(PRIMARY_RUN.read_text())
    independent = json.loads(INDEPENDENT_RUN.read_text())
    validate_run(probe, "primary", 16, True)
    validate_run(primary, "primary", 64, False)
    validate_run(independent, "independent", 64, False)

    full_runs = {"primary": primary, "independent": independent}
    row_maps = {
        mode: {
            (
                row["window_nodes"],
                int(row["component_support_mask_hex"], 16),
                row["packet_value"],
            ): row
            for row in run["rows"]
        }
        for mode, run in full_runs.items()
    }
    probe_rows = {
        (
            row["window_nodes"],
            int(row["component_support_mask_hex"], 16),
            row["packet_value"],
        ): row
        for row in probe["rows"]
    }

    apply_state = build_apply(systematic_state_columns())

    def step(state: int) -> int:
        return apply_state(accumulate(state))

    transpose_apply = build_apply(
        transpose_columns(tuple(step(1 << bit) for bit in range(64)))
    )
    state_cache = {
        (window, value): observations(step, value, window)
        for window, values in ((24, range(1, 16)), (12, (15,)))
        for value in values
    }

    certificate_rows = []
    full_count = 0
    probe_count = 0
    for mask, degrees in FRONTIER:
        if support_degrees(mask) != degrees or support_dimension(mask) != DIMENSION:
            raise RuntimeError("dimension-33 audit: support identity changed")
        annihilator = support_annihilator(mask)
        basis = kernel_basis(
            tuple(
                apply_polynomial(transpose_apply, annihilator, 1 << bit)
                for bit in range(64)
            )
        )
        if len(basis) != DIMENSION:
            raise RuntimeError("dimension-33 audit: kernel dimension changed")
        for window, values in ((24, range(1, 16)), (12, (15,))):
            for value in values:
                identity = (window, mask, value)
                columns = code_columns(state_cache[(window, value)], basis)
                receipts = {}
                for mode in ("primary", "independent"):
                    row = row_maps[mode][identity]
                    path = ROOT / row["receipt"]
                    if digest(path) != row["receipt_sha256"]:
                        raise RuntimeError("dimension-33 audit: receipt hash mismatch")
                    receipt = json.loads(path.read_text())
                    if (
                        int(receipt["component_support_mask_hex"], 16) != mask
                        or tuple(receipt["component_degrees"]) != degrees
                        or int(receipt["component_annihilator_hex"], 16)
                        != annihilator
                    ):
                        raise RuntimeError("dimension-33 audit: receipt identity mismatch")
                    audit_receipt(receipt, columns, mode)
                    receipts[mode] = receipt
                    full_count += 1
                if receipts["primary"]["information_sets"] == receipts["independent"]["information_sets"]:
                    raise RuntimeError("dimension-33 audit: verifier partitions coincide")
                certificate_rows.append(
                    {
                        "window_nodes": window,
                        "component_support_mask_hex": hex(mask),
                        "component_degrees": list(degrees),
                        "packet_value": value,
                        "primary_receipt_sha256": row_maps["primary"][identity]["receipt_sha256"],
                        "independent_receipt_sha256": row_maps["independent"][identity]["receipt_sha256"],
                        "partitions_are_distinct": True,
                        "passed": True,
                    }
                )
                if identity in probe_rows:
                    row = probe_rows[identity]
                    path = ROOT / row["receipt"]
                    if digest(path) != row["receipt_sha256"]:
                        raise RuntimeError("dimension-33 audit: probe receipt hash mismatch")
                    audit_receipt(json.loads(path.read_text()), columns, "primary")
                    probe_count += 1

    if full_count != 128 or probe_count != 16:
        raise RuntimeError("dimension-33 audit: receipt count changed")
    if primary["total_information_vectors_checked"] != EXPECTED_VECTORS:
        raise RuntimeError("dimension-33 audit: primary vector count changed")
    if independent["total_information_vectors_checked"] != EXPECTED_VECTORS:
        raise RuntimeError("dimension-33 audit: independent vector count changed")

    support_ledger = [
        {
            "component_support_mask_hex": hex(mask),
            "component_degrees": list(support_degrees(mask)),
            "dimension": support_dimension(mask),
            "status": (
                "CERTIFIED_BY_DIMENSION33_FRONTIER"
                if mask in frontier_masks
                else "CERTIFIED_BY_DIMENSION32_BASE"
            ),
        }
        for mask in sorted(target33)
    ]
    merged_masks = target33 | {0x60}
    payload = {
        "schema": "riffle-dp-2lap-g4-component-dimension33-audit-v1",
        "candidate": "Riffle DP-2Lap g=4",
        "evidence_label": "EXACT_EXHAUSTIVE_CERTIFICATE_WITH_INDEPENDENT_VERIFICATION",
        "audit_source_sha256": digest(Path(__file__).resolve()),
        "certificate_source_sha256": digest(SOURCE),
        "runner_source_sha256": digest(RUNNER),
        "base_runner_source_sha256": digest(BASE_RUNNER),
        "base_dimension32_audit_sha256": digest(BASE_AUDIT),
        "primary_probe_sha256": digest(PRIMARY_PROBE),
        "primary_run_sha256": digest(PRIMARY_RUN),
        "independent_run_sha256": digest(INDEPENDENT_RUN),
        "target_maximum_dimension": 33,
        "base_support_count": len(target32),
        "frontier_support_count": len(frontier_masks),
        "target_support_count": len(target33),
        "frontier_support_masks_hex": [hex(mask) for mask, _ in FRONTIER],
        "certificate_case_count_per_mode": len(certificate_rows),
        "audited_full_receipt_count": full_count,
        "audited_probe_receipt_count": probe_count,
        "information_vectors_checked_per_full_mode": EXPECTED_VECTORS,
        "information_vectors_checked_both_full_modes": 2 * EXPECTED_VECTORS,
        "certificate_rows": certificate_rows,
        "support_ledger": support_ledger,
        "merged_component_coverage": {
            "all_nonempty_component_supports": 127,
            "dimension_at_most_33_supports": len(target33),
            "additional_c18_c20_support_mask_hex": "0x60",
            "distinct_certified_supports": len(merged_masks),
            "remaining_supports": 127 - len(merged_masks),
        },
        "result": "PASS",
        "scope_limitation": (
            "The endpoint-aware certificate covers all 69 supports of total "
            "dimension at most 33. The separate C18+C20 result raises merged "
            "coverage to 70 supports. The other 57 nonempty supports remain open."
        ),
    }
    OUTPUT.write_text(json.dumps(payload, indent=2) + "\n")
    print(f"target_supports={len(target33)}")
    print(f"audited_full_receipts={full_count}")
    print(f"audited_probe_receipts={probe_count}")
    print(f"merged_certified_supports={len(merged_masks)}")
    print(f"output={OUTPUT}")
    print("status=PASS")


if __name__ == "__main__":
    main()
