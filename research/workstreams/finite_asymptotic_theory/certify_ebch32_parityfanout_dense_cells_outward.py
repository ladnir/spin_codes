#!/usr/bin/env python3
"""Arb outward verifier for saved EBCH32--PF31x33 dense-cover cells.

This verifier can check a dominant prefix while the diagnostic cover remains
incomplete.  A partial run is not a distance certificate.  Witness values
parsed from JSON are treated as exact IEEE-754 binary64 rationals.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import math
import os
from pathlib import Path

import flint
from flint import arb, ctx
import numpy as np

from analyze_golay_ba_rm2sub_joint import DEFAULT_SELECTION, load_envelope
from certify_ebch32_parityfanout_ba_setup import (
    B,
    LOWER_WEIGHT,
    UPPER_WEIGHT,
    conditioned_spectrum_upper,
)
from diagnose_ebch32_parityfanout_ba_three_band_cover_all_q import (
    BANDS,
    CELL_CONTRIBUTION_MARGIN,
    DELEGATED_EXTERNAL_TYPE,
    DISTANCE,
    EPOCHS,
    L,
    bounding_box_points,
)


WORKSTREAM = Path(__file__).resolve().parent
INPUT = WORKSTREAM / os.environ.get(
    "SPIN_EBCH_DENSE_COVER_OUTPUT",
    "ebch32_parityfanout31x33_ba3_B256_three_band_cover_all_q_d11.json",
)
OUTPUT = WORKSTREAM / "ebch32_parityfanout31x33_ba3_B256_dense_cells_outward_partial.json"
SETUP_VERIFIER = WORKSTREAM / "certify_ebch32_parityfanout_ba_setup.py"
SETUP_RECEIPT = WORKSTREAM / "ebch32_parityfanout31x33_ba3_B256_setup_outward.json"
SPECTRUM_SOURCE = WORKSTREAM / "ebch32_16_delta8_spectrum.csv"
COVER_GENERATOR = WORKSTREAM / "diagnose_ebch32_parityfanout_ba_three_band_cover_all_q.py"
PARTITION_AUDIT = WORKSTREAM / os.environ.get(
    "SPIN_EBCH_PARTITION_AUDIT_OUTPUT",
    "ebch32_parityfanout31x33_ba3_B256_dense_cover_partition_audit.json",
)
DELEGATED_OUTWARD = WORKSTREAM / "ebch32_parityfanout31x33_ba3_B256_split_low_rejected_point_cover_outward.json"
DELEGATED_VERIFIER = WORKSTREAM / "certify_ebch32_parityfanout_split_low_interval_cover_outward.py"


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def exact_float(value: float) -> arb:
    numerator, denominator = float(value).as_integer_ratio()
    return arb(numerator) / denominator


def definitely_positive(value: arb) -> bool:
    return bool(value.lower() > 0)


def definitely_negative(value: arb) -> bool:
    return bool(value.upper() < 0)


def arb_matrix_multiply(left, right):
    size = len(left)
    return [
        [sum((left[i][k] * right[k][j] for k in range(size)), arb(0)) for j in range(size)]
        for i in range(size)
    ]


def arb_matrix_power(matrix, exponent: int):
    size = len(matrix)
    result = [
        [arb(1) if i == j else arb(0) for j in range(size)]
        for i in range(size)
    ]
    factor = matrix
    power = exponent
    while power:
        if power & 1:
            result = arb_matrix_multiply(result, factor)
        power >>= 1
        if power:
            factor = arb_matrix_multiply(factor, factor)
    return result


class DenseOutwardVerifier:
    def __init__(self, precision_bits: int) -> None:
        ctx.prec = precision_bits
        self.precision_bits = precision_bits
        self.spectrum, self.tail_upper, self.good_lower = conditioned_spectrum_upper()
        self.envelope = load_envelope(DEFAULT_SELECTION, DISTANCE / (B * L))
        self.log2 = arb(2).log()
        self.band_cache: dict[tuple[str, str, int, int], arb] = {}
        self.inner_cache: dict[tuple[str, str], arb] = {}

    def band_log_moment(self, band, value_probability: float, order: float) -> arb:
        key = (float(value_probability).hex(), float(order).hex(), band[0], band[1])
        if key in self.band_cache:
            return self.band_cache[key]
        probability = exact_float(value_probability)
        one_minus = arb(1) - probability
        renyi_order = exact_float(order)
        if not definitely_positive(probability) or not definitely_positive(one_minus):
            raise ArithmeticError("invalid Bernoulli value witness")
        total = arb(0)
        for weight in range(band[0], band[1] + 1):
            multiplicity_upper = self.spectrum[weight]
            if multiplicity_upper == 0.0:
                continue
            multiplicity = exact_float(multiplicity_upper)
            reference = (
                arb(math.comb(B, weight))
                * probability**weight
                * one_minus ** (B - weight)
            )
            if not definitely_positive(reference):
                raise ArithmeticError("reference shell mass is not positive")
            term = (renyi_order * multiplicity.log() + (arb(1) - renyi_order) * reference.log()).exp()
            total += term
        if not definitely_positive(total):
            raise ArithmeticError("Renyi band moment is not positive")
        result = total.log()
        self.band_cache[key] = result
        return result

    def inner_log_moment(self, bit_probability: arb, surprisal: float) -> arb:
        key = (str(bit_probability), float(surprisal).hex())
        if key in self.inner_cache:
            return self.inner_cache[key]
        beta = bit_probability
        s = exact_float(surprisal)
        if not definitely_positive(beta) or not definitely_positive(arb(1) - beta):
            raise ArithmeticError("invalid actual-bit probability")
        z = (-s).exp()
        u = arb(1) - beta + beta * z
        v = beta + (arb(1) - beta) * z
        signed = arb(1) - beta - beta * z
        classes = [int(value) for value in self.envelope.classes]
        counts = [int(value) for value in self.envelope.counts]
        association = [
            [int(self.envelope.association[i, j]) for j in range(len(classes))]
            for i in range(len(classes))
        ]
        live = [u ** (self.envelope.step_bits - weight) * v**weight for weight in classes]
        fourier = [
            u ** (self.envelope.step_bits - weight) * signed**weight
            for weight in classes
        ]
        state_space = arb(self.envelope.state_space)
        zero_to_zero = sum(
            (arb(count) * value for count, value in zip(counts, fourier)), arb(0)
        ) / state_space
        zero_total = u**self.envelope.step_bits
        zero_to_deterministic = zero_total - zero_to_zero
        paired_total = sum(
            (
                fourier[i]
                * association[i][j]
                * live[j]
                for i in range(len(classes))
                for j in range(len(classes))
            ),
            arb(0),
        ) / state_space
        paired_nonzero = paired_total - zero_to_zero * live[0]
        uniform_numerator = sum(
            (arb(count) * value for count, value in zip(counts, live)), arb(0)
        ) - live[0]
        for name, value in (
            ("zero_to_zero", zero_to_zero),
            ("zero_to_deterministic", zero_to_deterministic),
            ("paired_nonzero", paired_nonzero),
            ("uniform_live_numerator", uniform_numerator),
        ):
            if not definitely_positive(value):
                raise ArithmeticError(f"{name} positivity is unresolved")
        deterministic = paired_nonzero / zero_to_deterministic
        uniform_live = uniform_numerator / self.envelope.live_states
        punctured = (
            arb(self.envelope.live_states)
            / (self.envelope.live_states - 1)
            * uniform_live
        )
        transfer = [
            [zero_to_zero, zero_to_deterministic, arb(0)],
            [deterministic / self.envelope.live_states, arb(0), deterministic],
            [punctured / self.envelope.live_states, arb(0), punctured],
        ]
        powered = arb_matrix_power(transfer, EPOCHS)
        row_sum = sum(powered[0], arb(0))
        if not definitely_positive(row_sum):
            raise ArithmeticError("powered transfer row sum is not positive")
        result = row_sum.log()
        self.inner_cache[key] = result
        return result

    def log_multinomial(self, counts: list[arb]) -> arb:
        total = arb(L + 1).lgamma()
        for count in counts:
            if bool(count.lower() < 0):
                raise ArithmeticError("negative type count")
            total -= (count + 1).lgamma()
        return total

    def verify_witness_at_vertices(self, witness, vertices) -> tuple[arb, list[str]]:
        probabilities = [exact_float(float(value)) for value in witness["reference_type_probabilities"]]
        values = [float(value) for value in witness["reference_value_probabilities"]]
        surprisal_float = float(witness["surprisal"])
        surprisal = exact_float(surprisal_float)
        order_float = float(witness["holder_order"])
        order = exact_float(order_float)
        if len(probabilities) != 4 or len(values) != 3:
            raise ArithmeticError("witness dimension mismatch")
        if not definitely_positive(order - 1):
            raise ArithmeticError("Renyi order is not above one")
        for probability in probabilities:
            if not definitely_positive(probability):
                raise ArithmeticError("reference type probability is not positive")
        probability_sum = sum(probabilities, arb(0))
        if not definitely_positive(probability_sum):
            raise ArithmeticError("reference type probability sum is not positive")
        probabilities = [probability / probability_sum for probability in probabilities]
        band_moments = [
            self.band_log_moment(band, value, order_float)
            for band, value in zip(BANDS, values)
        ]
        bit_probability_ball = sum(
            (probabilities[index + 1] * exact_float(values[index]) for index in range(3)),
            arb(0),
        )
        inner = self.inner_log_moment(bit_probability_ball, surprisal_float)
        conjugate = (order - 1) / order
        convex_coefficient = arb(1) - B * conjugate
        if not definitely_negative(convex_coefficient):
            raise ArithmeticError("convexity coefficient sign is unresolved")
        vertex_values = []
        raw_upper_strings = []
        for active in vertices:
            active_counts = [exact_float(float(value)) for value in active]
            occupation = sum(active_counts, arb(0))
            counts = [arb(L) - occupation, *active_counts]
            log_type_count = self.log_multinomial(counts)
            log_conditioning = log_type_count + sum(
                (count * probability.log() for count, probability in zip(counts, probabilities)),
                arb(0),
            )
            raw = inner + DISTANCE * surprisal - B * log_conditioning
            if not definitely_negative(raw):
                raise ArithmeticError("raw reference exponent is not proved negative")
            value = (
                log_type_count
                + sum(
                    (active_counts[index] * band_moments[index] for index in range(3)),
                    arb(0),
                )
                / order
                + conjugate * raw
            )
            vertex_values.append(value)
            raw_upper_strings.append(str(raw.upper()))
        maximum = vertex_values[0]
        for value in vertex_values[1:]:
            if bool(value.upper() > maximum.upper()):
                maximum = value
        return maximum, raw_upper_strings

    def verify_cell(self, row) -> tuple[dict[str, object], arb | None]:
        mode = row["optimizer_mode"]
        box_count = int(row["bounding_box_lattice_count"])
        if mode == "delegated-split-low-point":
            points = set(
                bounding_box_points(
                    np.asarray(row["vertices_active_counts"], dtype=np.float64)
                )
            )
            if points != {DELEGATED_EXTERNAL_TYPE} or box_count != 1:
                raise ArithmeticError("delegated cell is not the certified external point")
            return {
                "mode": mode,
                "depth": int(row["depth"]),
                "box_count": box_count,
                "delegated_external_type": list(DELEGATED_EXTERNAL_TYPE),
                "log2_contribution_upper": None,
                "log2_contribution_upper_display": -198.0,
                "cell_threshold_check": True,
                "excluded_from_this_sum": True,
            }, None
        if box_count == 0:
            return {
                "mode": mode,
                "box_count": 0,
                "log2_contribution_upper": None,
                "cell_threshold_check": True,
            }, None
        if mode == "enumerated-lattice-box":
            contribution = arb(0)
            point_rows = []
            for point in row["point_witnesses"]:
                maximum, raw_strings = self.verify_witness_at_vertices(
                    point, point["vertices_active_counts"][:1]
                )
                contribution += maximum.exp()
                point_rows.append({"raw_upper": raw_strings[0], "exponent_upper": str(maximum.upper())})
            log_contribution = contribution.log()
            raw_strings = []
        else:
            maximum, raw_strings = self.verify_witness_at_vertices(
                row, row["vertices_active_counts"]
            )
            log_contribution = arb(box_count).log() + maximum
            point_rows = []
        check = definitely_negative(
            log_contribution + CELL_CONTRIBUTION_MARGIN * self.log2
        ) or bool((log_contribution + CELL_CONTRIBUTION_MARGIN * self.log2).upper() == 0)
        if not check:
            raise ArithmeticError("cell contribution does not pass the configured threshold")
        return {
            "mode": mode,
            "depth": int(row["depth"]),
            "box_count": box_count,
            "log_contribution_upper": str(log_contribution.upper()),
            "log2_contribution_upper_display": float(log_contribution.upper()) / math.log(2.0),
            "cell_threshold_check": True,
            "raw_vertex_upper_bounds": raw_strings,
            "point_rows": point_rows,
        }, log_contribution


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--maximum-cells", type=int, default=20)
    parser.add_argument("--precision-bits", type=int, default=256)
    parser.add_argument("--output", type=Path, default=OUTPUT)
    args = parser.parse_args()
    if args.maximum_cells <= 0:
        raise ValueError("--maximum-cells must be positive")
    source = json.loads(INPUT.read_text(encoding="utf-8"))
    candidates = [
        row for row in source["tetrahedra"]
        if int(row["bounding_box_lattice_count"]) > 0
    ]
    candidates.sort(
        key=lambda row: float(row["bounding_box_log2_contribution_upper"]),
        reverse=True,
    )
    selected = candidates[: args.maximum_cells]
    verifier = DenseOutwardVerifier(args.precision_bits)
    rows = []
    aggregate = arb(0)
    for index, row in enumerate(selected, 1):
        result, log_contribution = verifier.verify_cell(row)
        if log_contribution is not None:
            aggregate += log_contribution.exp()
        result["diagnostic_rank"] = index
        result["diagnostic_log2_contribution"] = float(
            row["bounding_box_log2_contribution_upper"]
        )
        rows.append(result)
        if index <= 20 or index % 100 == 0 or index == len(selected):
            print(
                f"cell={index}/{len(selected)},mode={result['mode']},"
                f"outward_log2={result['log2_contribution_upper_display']:.9f}",
                flush=True,
            )
    aggregate_log = aggregate.log() if definitely_positive(aggregate) else None
    all_nonempty = len(selected) == len(candidates)
    delegated_selected = sum(
        row["mode"] == "delegated-split-low-point" for row in rows
    )
    payload = {
        "schema": "ebch32-parityfanout31x33-b256-dense-cells-arb-outward-partial-v1",
        "status": "PARTIAL_OUTWARD_CELL_VERIFICATION",
        "claim": {
            "verified_dominant_accepted_cells": len(rows),
            "source_accepted_cells": int(source["accepted_tetrahedra"]),
            "source_pending_cells": int(source["pending_tetrahedra"]),
            "source_rejected_cells": len(source["rejected_tetrahedra"]),
            "all_selected_cells_pass_configured_threshold": all(
                bool(row["cell_threshold_check"]) for row in rows
            ),
            "all_nonempty_accepted_cells_verified": all_nonempty,
            "delegated_cells_verified_geometrically": delegated_selected,
            "delegated_external_type_added_separately": list(DELEGATED_EXTERNAL_TYPE),
            "accepted_cells_log2_sum_upper": (
                None
                if aggregate_log is None
                else str((aggregate_log / verifier.log2).upper())
            ),
            "accepted_cells_margin_bits_lower_display": (
                None
                if aggregate_log is None
                else -float((aggregate_log / verifier.log2).upper())
            ),
            "full_dense_certificate": False,
        },
        "arithmetic": {
            "library": "python-flint Arb",
            "library_version": flint.__version__,
            "precision_bits": args.precision_bits,
            "witness_interpretation": "JSON binary64 values are converted to exact integer ratios",
        },
        "parameters": {
            "outer_bits": B,
            "outer_rows": L,
            "distance": DISTANCE,
            "epochs": EPOCHS,
            "bands": [list(band) for band in BANDS],
            "permitted_outer_weights": [LOWER_WEIGHT, UPPER_WEIGHT],
        },
        "outer_conditioning": {
            "tail_upper_hex": verifier.tail_upper.hex(),
            "good_probability_lower_hex": verifier.good_lower.hex(),
        },
        "cache": {
            "band_moments": len(verifier.band_cache),
            "inner_moments": len(verifier.inner_cache),
        },
        "rows": rows,
        "dependencies": [
            {"path": INPUT.name, "sha256": sha256(INPUT)},
            {"path": Path(__file__).name, "sha256": sha256(Path(__file__).resolve())},
            {"path": SETUP_VERIFIER.name, "sha256": sha256(SETUP_VERIFIER)},
            {"path": SETUP_RECEIPT.name, "sha256": sha256(SETUP_RECEIPT)},
            {"path": SPECTRUM_SOURCE.name, "sha256": sha256(SPECTRUM_SOURCE)},
            {"path": COVER_GENERATOR.name, "sha256": sha256(COVER_GENERATOR)},
            {"path": PARTITION_AUDIT.name, "sha256": sha256(PARTITION_AUDIT)},
            {"path": DELEGATED_OUTWARD.name, "sha256": sha256(DELEGATED_OUTWARD)},
            {"path": DELEGATED_VERIFIER.name, "sha256": sha256(DELEGATED_VERIFIER)},
            {"path": DEFAULT_SELECTION.as_posix(), "sha256": sha256(DEFAULT_SELECTION)},
        ],
        "limitations": [
            "Only the selected dominant accepted cells are verified.",
            "The diagnostic cover still has pending cells.",
            "Delegated cells are geometry-checked here but excluded from this sum; the split-low receipt must be added once.",
            "The aggregate dense sum and the sparse-plus-dense total are not verified here.",
        ],
    }
    args.output.write_text(json.dumps(payload, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(payload["claim"], indent=2))
    print(f"output={args.output}")


if __name__ == "__main__":
    main()
