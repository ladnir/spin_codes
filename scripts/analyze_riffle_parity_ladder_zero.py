#!/usr/bin/env python3
"""Count packet histograms for the zero-parity Riffle baseline.

The count is an ensemble expectation over the independent uniform bit
permutation attached to every BCH block.  The global packet permutation does
not change a packet-weight histogram and therefore does not enter this stage.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
from collections import defaultdict
from fractions import Fraction
from pathlib import Path


SCRIPT_DIRECTORY = Path(__file__).resolve().parent
REPOSITORY_ROOT = SCRIPT_DIRECTORY.parent
DEFAULT_SPECTRUM = SCRIPT_DIRECTORY / "ebch128_64_spectrum.csv"
DEFAULT_OUTPUT = (
    REPOSITORY_ROOT
    / "constructions"
    / "riffle_shiftalpha64_bchblockperm_parallelacc_g4"
    / "receipts"
    / "goal17_zero_parity_histogram_count.json"
)

DATA_BLOCKS = 1 << 14
PACKET_WIDTH = 4
PACKETS_PER_BCH_BLOCK = 32
BCH_LENGTH = PACKET_WIDTH * PACKETS_PER_BCH_BLOCK
BCH_DIMENSION = 64
BCH_MINIMUM_WEIGHT = 22


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as source:
        for chunk in iter(lambda: source.read(1 << 20), b""):
            digest.update(chunk)
    return digest.hexdigest().upper()


def load_spectrum(path: Path) -> dict[int, int]:
    with path.open(newline="", encoding="utf-8") as source:
        spectrum = {
            int(row["weight"]): int(row["count"])
            for row in csv.DictReader(source)
        }
    if sum(spectrum.values()) != 1 << BCH_DIMENSION:
        raise RuntimeError("EBCH spectrum has the wrong total mass")
    if spectrum.get(0) != 1 or spectrum.get(BCH_MINIMUM_WEIGHT) != 243_840:
        raise RuntimeError("EBCH spectrum failed its endpoint checks")
    for weight, count in spectrum.items():
        if spectrum.get(BCH_LENGTH - weight, 0) != count:
            raise RuntimeError("EBCH spectrum failed complement symmetry")
    return spectrum


def weak_compositions(total: int, parts: int):
    if parts == 1:
        yield (total,)
        return
    for first in range(total + 1):
        for suffix in weak_compositions(total - first, parts - 1):
            yield (first,) + suffix


def packetization_mass(histogram: tuple[int, ...]) -> int:
    """Count 128-bit strings with the prescribed four-bit packet histogram."""
    ways = math.factorial(PACKETS_PER_BCH_BLOCK)
    for count in histogram:
        ways //= math.factorial(count)
    for packet_weight, count in enumerate(histogram):
        ways *= math.comb(PACKET_WIDTH, packet_weight) ** count
    return ways


def weight_slice(weight: int) -> list[tuple[tuple[int, ...], int]]:
    rows = []
    for histogram in weak_compositions(PACKETS_PER_BCH_BLOCK, PACKET_WIDTH + 1):
        if sum(index * count for index, count in enumerate(histogram)) != weight:
            continue
        rows.append((histogram, packetization_mass(histogram)))
    if sum(mass for _, mass in rows) != math.comb(BCH_LENGTH, weight):
        raise RuntimeError(f"packetization mass mismatch at binary weight {weight}")
    return rows


def rational_payload(value: Fraction) -> dict[str, int | str]:
    return {
        "numerator": value.numerator,
        "denominator": value.denominator,
        "decimal": format(float(value), ".17g"),
    }


def minimum_shell_payload(spectrum: dict[int, int]) -> dict[str, object]:
    weight = BCH_MINIMUM_WEIGHT
    bch_words = spectrum[weight]
    labeled_outer_words = DATA_BLOCKS * bch_words
    slice_size = math.comb(BCH_LENGTH, weight)
    global_packet_positions = DATA_BLOCKS * PACKETS_PER_BCH_BLOCK
    rows = []
    support_mass: defaultdict[int, int] = defaultdict(int)
    expected_sum = Fraction(0)

    for local_histogram, mass in weight_slice(weight):
        active_support = sum(local_histogram[1:])
        support_mass[active_support] += mass
        expected_count = Fraction(labeled_outer_words * mass, slice_size)
        expected_sum += expected_count
        global_histogram = (
            global_packet_positions - PACKETS_PER_BCH_BLOCK + local_histogram[0],
            *local_histogram[1:],
        )
        rows.append(
            {
                "local_histogram_h0_through_h4": list(local_histogram),
                "global_histogram_h0_through_h4": list(global_histogram),
                "active_packet_support": active_support,
                "weight_slice_outcomes": mass,
                "histogram_probability": rational_payload(Fraction(mass, slice_size)),
                "expected_outer_words": rational_payload(expected_count),
            }
        )

    if expected_sum != labeled_outer_words:
        raise RuntimeError("minimum-shell expected counts do not recover shell size")

    most_likely = max(rows, key=lambda row: int(row["weight_slice_outcomes"]))
    support_rows = []
    support_probability_sum = Fraction(0)
    for support in sorted(support_mass):
        probability = Fraction(support_mass[support], slice_size)
        support_probability_sum += probability
        support_rows.append(
            {
                "active_packet_support": support,
                "weight_slice_outcomes": support_mass[support],
                "probability": rational_payload(probability),
                "expected_outer_words": rational_payload(
                    labeled_outer_words * probability
                ),
            }
        )
    if support_probability_sum != 1:
        raise RuntimeError("packet-support probabilities do not sum to one")

    return {
        "binary_weight": weight,
        "bch_words_per_data_position": bch_words,
        "data_positions": DATA_BLOCKS,
        "labeled_outer_words": labeled_outer_words,
        "labeled_outer_words_log2": math.log2(labeled_outer_words),
        "binary_weight_slice_size": slice_size,
        "distinct_packet_histograms": len(rows),
        "active_packet_support_range": [
            min(support_mass),
            max(support_mass),
        ],
        "most_likely_local_histogram": most_likely,
        "packet_support_distribution": support_rows,
        "packet_histograms": rows,
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--spectrum", type=Path, default=DEFAULT_SPECTRUM)
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    args = parser.parse_args()

    spectrum = load_spectrum(args.spectrum)
    occupied_one_all_weights = DATA_BLOCKS * ((1 << BCH_DIMENSION) - 1)
    packet_positions = DATA_BLOCKS * PACKETS_PER_BCH_BLOCK
    minimum_shell = minimum_shell_payload(spectrum)
    payload = {
        "schema": "riffle-parity-ladder-zero-histogram-count-v1",
        "evidence_label": "EXACT_ENSEMBLE_EXPECTATION",
        "stage": {
            "parity_blocks": 0,
            "outer_field_code": f"[{DATA_BLOCKS},{DATA_BLOCKS},1] identity code",
            "minimum_field_symbol_occupation": 1,
        },
        "parameters": {
            "data_blocks": DATA_BLOCKS,
            "bch": f"[{BCH_LENGTH},{BCH_DIMENSION},{BCH_MINIMUM_WEIGHT}] extended binary BCH",
            "packet_width": PACKET_WIDTH,
            "packets_per_bch_block": PACKETS_PER_BCH_BLOCK,
            "global_packet_positions": packet_positions,
            "global_binary_length": packet_positions * PACKET_WIDTH,
        },
        "complete_count": {
            "local_expected_histogram_polynomial": (
                "F(t)=sum_r A_r/binom(128,r) * "
                "sum_{h:|h|=32,||h||=r} P_r(h)t^h"
            ),
            "all_messages_expected_histogram_polynomial": (
                f"F(t)^{DATA_BLOCKS}"
            ),
            "nonzero_messages_expected_histogram_polynomial": (
                f"F(t)^{DATA_BLOCKS}-t0^{packet_positions}"
            ),
            "total_messages_symbolic": f"2^(64*{DATA_BLOCKS})",
            "nonzero_messages_symbolic": f"2^(64*{DATA_BLOCKS})-1",
            "occupation_one_outer_words": occupied_one_all_weights,
            "occupation_one_outer_words_log2": math.log2(occupied_one_all_weights),
        },
        "minimum_shell": minimum_shell,
        "accumulator_interface": {
            "conditional_count": "A_N(h,w)/T_N(h)",
            "expected_output_enumerator": (
                "sum_h C_0(h) A_N(h,w)/T_N(h)"
            ),
            "global_histogram_convention": (
                "The h0 coordinate includes all packets in inactive data blocks."
            ),
        },
        "validation": {
            "exact_bch_spectrum_mass": "PASS",
            "bch_complement_symmetry": "PASS",
            "minimum_slice_packetization_mass": "PASS",
            "minimum_shell_expected_count_mass": "PASS",
            "packet_support_probability_mass": "PASS",
        },
        "inputs": {
            "spectrum": str(args.spectrum.relative_to(REPOSITORY_ROOT)).replace("\\", "/"),
            "spectrum_sha256": sha256(args.spectrum),
        },
        "scope": (
            "This stage computes C_0(h), not the accumulator output sum. "
            "The full count is represented by an exact generating function. "
            "The complete minimum shell is also expanded coefficient by coefficient."
        ),
    }

    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print(json.dumps({
        "output": str(args.output),
        "minimum_shell_outer_words": minimum_shell["labeled_outer_words"],
        "minimum_shell_packet_histograms": minimum_shell["distinct_packet_histograms"],
        "status": "PASS",
    }, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
