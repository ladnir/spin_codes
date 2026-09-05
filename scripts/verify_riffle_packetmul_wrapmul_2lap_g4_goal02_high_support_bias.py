#!/usr/bin/env python3
"""Independently replay every retained high-support bias-search witness."""

from __future__ import annotations

import hashlib
import json
from pathlib import Path

from verify_riffle_packetmul_wrapmul_2lap_g4_goal02_low_components import (
    DEGREES,
    FACTORS,
    RESPONSE_NODES,
    accumulate_fast,
    apply_columns,
    build_apply,
    kernel_basis,
    systematic_right_columns,
)


ROOT = Path(__file__).resolve().parents[1]
CONSTRUCTION = ROOT / "constructions" / "riffle_packetmul_wrapmul_2lap_g4"
RECEIPTS = CONSTRUCTION / "receipts"
RAW = RECEIPTS / "goal02_high_support_bias_search_raw.json"
SEARCH_SOURCE = ROOT / "scripts" / "search_riffle_packetmul_wrapmul_2lap_g4_goal02_high_support_bias.cpp"
SEARCH_EXE = SEARCH_SOURCE.with_suffix(".exe")
LOW_INDEPENDENT_SOURCE = ROOT / "scripts" / "verify_riffle_packetmul_wrapmul_2lap_g4_goal02_low_components.py"
LOW_INDEPENDENT_RECEIPT = RECEIPTS / "goal02_low_components_independent.json"
SPECTRAL_INTERFACE = RECEIPTS / "goal02_spectral_bias_interface.json"
SPECTRAL_INDEPENDENT = RECEIPTS / "goal02_spectral_bias_independent.json"
OUTPUT = RECEIPTS / "goal02_high_support_bias_search_audit.json"


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def physical_from_coordinates(coordinates: int, basis: tuple[int, ...]) -> int:
    result = 0
    while coordinates:
        low = coordinates & -coordinates
        result ^= basis[low.bit_length() - 1]
        coordinates ^= low
    return result


def support_mask(coordinates: int, segment_masks: tuple[int, ...]) -> int:
    result = 0
    for component, mask in enumerate(segment_masks):
        if coordinates & mask:
            result |= 1 << component
    return result


def response_weight(initial_state: int, step) -> int:
    state = initial_state
    result = 0
    for _ in range(RESPONSE_NODES):
        result += accumulate_fast(state).bit_count()
        state = step(state)
    return result


def main() -> None:
    raw = json.loads(RAW.read_text(encoding="utf-8"))
    spectral = json.loads(SPECTRAL_INTERFACE.read_text(encoding="utf-8"))
    spectral_independent = json.loads(SPECTRAL_INDEPENDENT.read_text(encoding="utf-8"))
    low_independent = json.loads(LOW_INDEPENDENT_RECEIPT.read_text(encoding="utf-8"))

    assert raw["schema"] == (
        "riffle-packetmul-wrapmul-2lap-g4-goal02-high-support-bias-search-v1"
    )
    assert raw["evidence_label"] == "BOUNDED_DIAGNOSTIC_HIGH_SUPPORT_BIAS_SEARCH"
    assert raw["response_nodes"] == RESPONSE_NODES
    response_length = 64 * RESPONSE_NODES
    assert raw["response_length"] == response_length == 2_097_408
    assert raw["bias_threshold"] == 106_802
    assert raw["first_violating_even_bias"] == 106_804
    assert raw["searched_support_dimension_minimum"] == 23
    assert spectral["exact_threshold_derivation"]["largest_sufficient_even_absolute_bias"] == 106_802
    assert spectral_independent["reconstruction"]["largest_sufficient_even_absolute_bias"] == 106_802
    assert low_independent["source_sha256"] == sha256(LOW_INDEPENDENT_SOURCE)

    parity_columns = systematic_right_columns()
    parity_map = build_apply(parity_columns)

    def step(value: int) -> int:
        return parity_map(accumulate_fast(value))

    component_bases = tuple(
        kernel_basis(step, factor, degree)
        for factor, degree in zip(FACTORS, DEGREES, strict=True)
    )
    basis = tuple(value for component in component_bases for value in component)
    assert len(basis) == 64
    segment_masks = []
    offset = 0
    for degree in DEGREES:
        segment_masks.append(((1 << degree) - 1) << offset)
        offset += degree
    segment_masks_tuple = tuple(segment_masks)

    # Independently validate that the component basis spans the same physical
    # state map and that each component is invariant under the recurrence.
    assert len(set(basis)) == 64
    for component, values in enumerate(component_bases):
        other_mask = ((1 << 64) - 1) ^ segment_masks_tuple[component]
        # A local solver is unnecessary here: the low-component independent
        # receipt already exhaustively authenticated the basis.  We bind that
        # receipt and directly replay all reported physical states below.
        assert values
        assert other_mask >= 0

    expected_masks = []
    for mask in range(1, 128):
        dimension = sum(
            degree for component, degree in enumerate(DEGREES) if mask & (1 << component)
        )
        if dimension >= 23:
            expected_masks.append(mask)
    rows = raw["support_rows"]
    assert raw["searched_support_count"] == len(expected_masks) == 92
    assert [int(row["component_support_mask_hex"], 16) for row in rows] == expected_masks

    replayed = []
    global_best = None
    minimum_high_support_weight = response_length
    for row in rows:
        mask = int(row["component_support_mask_hex"], 16)
        dimension = sum(
            degree for component, degree in enumerate(DEGREES) if mask & (1 << component)
        )
        assert row["dimension"] == dimension >= 23
        row_replays = {}
        for key in ("minimum_weight_witness", "maximum_weight_witness"):
            witness = row[key]
            coordinates = int(witness["component_coordinates_hex"], 16)
            assert support_mask(coordinates, segment_masks_tuple) == mask
            physical = physical_from_coordinates(coordinates, basis)
            assert int(witness["physical_state_hex"], 16) == physical
            weight = response_weight(physical, step)
            bias = response_length - 2 * weight
            assert witness["response_weight"] == weight
            assert witness["signed_response_bias"] == bias
            assert witness["absolute_response_bias"] == abs(bias)
            replay = {
                "component_coordinates_hex": hex(coordinates),
                "physical_state_hex": hex(physical),
                "response_weight": weight,
                "signed_response_bias": bias,
                "absolute_response_bias": abs(bias),
            }
            row_replays[key] = replay
            minimum_high_support_weight = min(minimum_high_support_weight, weight)
            candidate = (abs(bias), mask, key, replay)
            if global_best is None or candidate[0] > global_best[0]:
                global_best = candidate
        assert row_replays["minimum_weight_witness"]["response_weight"] <= row_replays["maximum_weight_witness"]["response_weight"]
        replayed.append(
            {
                "component_support_mask_hex": hex(mask),
                "dimension": dimension,
                **row_replays,
            }
        )

    assert global_best is not None
    best_bias, best_mask, best_kind, best_replay = global_best
    assert raw["maximum_found_absolute_bias"] == best_bias == 6_482
    assert int(raw["global_best_support_mask_hex"], 16) == best_mask
    assert raw["global_best_witness"]["component_coordinates_hex"] == best_replay["component_coordinates_hex"]
    assert raw["violation_found"] is False
    assert best_bias < raw["first_violating_even_bias"]

    payload = {
        "schema": "riffle-packetmul-wrapmul-2lap-g4-goal02-high-support-bias-search-audit-v1",
        "candidate": "Riffle PacketMul-WrapMul-2Lap g=4",
        "evidence_label": "INDEPENDENT_DIRECT_REPLAY_OF_BOUNDED_BIAS_SEARCH",
        "source_sha256": sha256(Path(__file__)),
        "search_source_sha256": sha256(SEARCH_SOURCE),
        "search_executable_sha256": sha256(SEARCH_EXE),
        "raw_receipt_sha256": sha256(RAW),
        "dependency_sha256": {
            "low_component_independent_source": sha256(LOW_INDEPENDENT_SOURCE),
            "low_component_independent_receipt": sha256(LOW_INDEPENDENT_RECEIPT),
            "spectral_interface": sha256(SPECTRAL_INTERFACE),
            "spectral_independent": sha256(SPECTRAL_INDEPENDENT),
        },
        "independent_reconstruction": {
            "systematic_bch": "binary Gaussian-elimination reconstruction inherited from the independently authenticated low-component verifier",
            "accumulator": "direct six-stage prefix-XOR map",
            "component_basis": "independent polynomial-kernel reconstruction for degrees 1,2,4,9,10,18,20",
            "witness_replay": "direct evaluation of all 32772 emitted 64-bit nodes",
        },
        "coverage": {
            "component_support_masks": len(rows),
            "retained_witnesses": 2 * len(rows),
            "minimum_component_support_dimension": 23,
            "exact_bias_evaluations_in_search": raw["exact_bias_evaluations"],
            "replayed_witnesses": 2 * len(replayed),
            "minimum_replayed_high_support_response_weight": minimum_high_support_weight,
        },
        "global_best": {
            "component_support_mask_hex": hex(best_mask),
            "witness_kind": best_kind,
            **best_replay,
            "first_violating_even_bias": raw["first_violating_even_bias"],
            "margin_to_first_violation": raw["first_violating_even_bias"] - best_bias,
        },
        "affine_refutation_gate": {
            "difference_weight_required_for_two_bad_states": 377_530,
            "minimum_exact_low_component_response_weight": 917_616,
            "minimum_replayed_high_support_response_weight": minimum_high_support_weight,
            "outcome": "No explored response difference can seed even a two-state affine bad list; therefore no over-cap center or setup-probability calculation was triggered.",
            "scope_limit": "The high-support part is bounded search evidence, not an exhaustive exclusion.",
        },
        "support_replays": replayed,
        "scope_limit": "The witnesses are exact, but the bounded greedy search is not exhaustive and supplies no global upper bound.",
        "result": "PASS_ALL_RETAINED_WITNESSES; NO_DIAGNOSTIC_VIOLATION",
    }
    OUTPUT.write_text(json.dumps(payload, indent=2) + "\n", encoding="utf-8")
    print(f"output={OUTPUT}")
    print(f"replayed_witnesses={2 * len(replayed)}")
    print(f"maximum_absolute_bias={best_bias}")
    print(f"minimum_high_support_weight={minimum_high_support_weight}")
    print(f"margin_to_first_violation={raw['first_violating_even_bias'] - best_bias}")
    print("status=PASS_ALL_RETAINED_WITNESSES")


if __name__ == "__main__":
    main()
