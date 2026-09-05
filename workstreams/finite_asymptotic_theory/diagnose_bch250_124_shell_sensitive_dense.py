#!/usr/bin/env python3
"""Optimize a shell-sensitive dense-occupation bound for wrapped BCH250.

Fix ``Q`` active outer rows.  A Bernoulli value law with parameter ``theta``
majorizes the expected wrapped-row counting measure by ``M(theta)``.  An
independent Bernoulli candidate mask with parameter ``r`` is conditioned to
contain exactly ``Q`` positions in every transposed region.  RM2Sub therefore
sees Bernoulli input bits with probability ``r*theta``.  The resulting bound is

    C(L,Q) M(theta)^Q Bin(L,Q;r)^(-B)
    z^(-D) c(r,theta,z) lambda(r,theta,z)^E.

The program optimizes ``r``, ``theta``, and ``z``.  It uses nearest binary64
arithmetic and emits a diagnostic, not an outward certificate.
"""

from __future__ import annotations

import argparse
import json
import math
from pathlib import Path
import sys

import numpy as np
from scipy.optimize import minimize


WORKSTREAM = Path(__file__).resolve().parent
sys.path.insert(0, str(WORKSTREAM))

from analyze_golay_ba_rm2sub_joint import (  # noqa: E402
    DEFAULT_SELECTION,
    load_envelope,
    logistic,
    logit,
)


B = 250
L = 8576
N = B * L
EPOCHS = N // 128
DISTANCE = 235_840
DEFAULT_OLD_COVER = (
    WORKSTREAM
    / "bch250_124_parityfanout31x33_l256_finite_interval_cover_q32_8576_d11_diagnostic.json"
)


def log_choose(total: int, selected: int) -> float:
    return (
        math.lgamma(total + 1)
        - math.lgamma(selected + 1)
        - math.lgamma(total - selected + 1)
    )


def log_binomial_mass(
    total: int, selected: int, probability: float
) -> float:
    result = log_choose(total, selected)
    if selected:
        result += selected * math.log(probability)
    if selected < total:
        result += (total - selected) * math.log1p(-probability)
    return result


def load_spectrum(path: Path) -> dict[int, float]:
    payload = json.loads(path.read_text(encoding="utf-8"))
    result: dict[int, float] = {}
    for row in payload["spectrum"]:
        weight = int(row["weight"])
        if weight == 0:
            continue
        value = row.get("log2_expected_multiplicity")
        if value is not None:
            result[weight] = float(value) * math.log(2.0)
    if not result:
        raise ValueError("spectrum has no positive nonzero shell")
    return result


def log_majorant(spectrum: dict[int, float], theta: float) -> tuple[float, int]:
    best = -math.inf
    best_weight = -1
    for weight, log_multiplicity in spectrum.items():
        value = log_multiplicity - log_choose(B, weight)
        value -= weight * math.log(theta)
        value -= (B - weight) * math.log1p(-theta)
        if value > best:
            best = value
            best_weight = weight
    return best, best_weight


def collatz_data(transfer: np.ndarray) -> tuple[float, float, list[float]]:
    if np.min(transfer) < -1.0e-10:
        raise ArithmeticError("transfer has a material negative entry")
    transfer = np.maximum(transfer, 0.0)
    eigenvalues, eigenvectors = np.linalg.eig(transfer)
    index = int(np.argmax(np.abs(eigenvalues)))
    vector = np.abs(np.real(eigenvectors[:, index]))
    vector = np.maximum(vector, 1.0e-300)
    vector /= float(np.max(vector))
    row_ratios = transfer @ vector / vector
    radius = max(float(abs(eigenvalues[index])), float(np.max(row_ratios)))
    radius = math.nextafter(radius, math.inf)
    prefactor = float(vector[0] * np.max(1.0 / vector))
    return radius, prefactor, [float(value) for value in vector]


def evaluate(
    point: np.ndarray,
    *,
    envelope,
    spectrum: dict[int, float],
    occupation: int,
) -> tuple[float, dict[str, object]]:
    r = logistic(float(point[0]))
    theta = logistic(float(point[1]))
    log_surprisal = min(2.0, max(-12.0, float(point[2])))
    surprisal = math.exp(log_surprisal)
    actual_bit_probability = r * theta
    if not 0.0 < actual_bit_probability < 0.5:
        return math.inf, {}

    transfer, _ = envelope.transfer(
        candidate_probability=2.0 * actual_bit_probability,
        surprisal=surprisal,
    )
    radius, prefactor, vector = collatz_data(transfer)
    majorant, maximizing_weight = log_majorant(spectrum, theta)
    value = log_choose(L, occupation)
    value += occupation * majorant
    value -= B * log_binomial_mass(L, occupation, r)
    value += DISTANCE * surprisal
    value += EPOCHS * math.log(radius)
    value += math.log(prefactor)
    return value, {
        "candidate_probability": r,
        "value_probability": theta,
        "actual_bit_probability": actual_bit_probability,
        "surprisal": surprisal,
        "z": math.exp(-surprisal),
        "majorant_log2": majorant / math.log(2.0),
        "majorant_maximizing_weight": maximizing_weight,
        "radius": radius,
        "matrix_prefactor": prefactor,
        "collatz_vector": vector,
    }


def old_witnesses(path: Path) -> list[dict[str, object]]:
    payload = json.loads(path.read_text(encoding="utf-8"))
    return list(payload["intervals"])


def old_start(rows: list[dict[str, object]], occupation: int) -> tuple[float, float]:
    for row in rows:
        if int(row["lower_occupation"]) <= occupation <= int(
            row["upper_occupation"]
        ):
            r = float(row["candidate_probability"])
            z = float(row["z"])
            return r, -math.log(z)
    raise ValueError("old interval cover does not contain occupation")


def optimize_occupation(
    envelope,
    spectrum: dict[int, float],
    old_rows: list[dict[str, object]],
    occupation: int,
) -> dict[str, object]:
    old_r, old_surprisal = old_start(old_rows, occupation)
    alpha = occupation / L
    starts = []
    for r in (old_r, alpha, min(0.999, max(0.001, 1.1 * alpha))):
        r = min(1.0 - 1.0e-12, max(1.0e-12, r))
        for theta in (0.35, 0.45, 0.5, 0.55, 0.65):
            if r * theta >= 0.5:
                theta = min(theta, 0.499999 / r)
            for scale in (0.7, 1.0, 1.4):
                surprisal = max(1.0e-8, old_surprisal * scale)
                starts.append(
                    np.asarray(
                        [logit(r), logit(theta), math.log(surprisal)],
                        dtype=np.float64,
                    )
                )

    best = None
    for start in starts:
        result = minimize(
            lambda point: evaluate(
                point,
                envelope=envelope,
                spectrum=spectrum,
                occupation=occupation,
            )[0],
            start,
            method="Nelder-Mead",
            options={"maxiter": 1000, "xatol": 1.0e-9, "fatol": 1.0e-7},
        )
        if best is None or float(result.fun) < float(best.fun):
            best = result
    assert best is not None
    value, details = evaluate(
        best.x,
        envelope=envelope,
        spectrum=spectrum,
        occupation=occupation,
    )
    return {
        "occupation": occupation,
        "log2_upper": value / math.log(2.0),
        "margin_bits": -value / math.log(2.0),
        "optimizer_success": bool(best.success),
        "optimizer_message": str(best.message),
        **details,
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--spectrum", type=Path, required=True)
    parser.add_argument("--selection", type=Path, default=DEFAULT_SELECTION)
    parser.add_argument("--old-cover", type=Path, default=DEFAULT_OLD_COVER)
    parser.add_argument("--fanout-layers", type=int, required=True)
    parser.add_argument(
        "--occupations",
        type=int,
        nargs="+",
        default=(32, 48, 64, 96, 128, 192, 256, 384, 512, 768, 1024,
                 1536, 2048, 3072, 4096, 6144, 8192, 8448, L),
    )
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    if any(not 32 <= value <= L for value in args.occupations):
        parser.error(f"occupations must lie in [32,{L}]")

    envelope = load_envelope(args.selection, DISTANCE / N)
    spectrum = load_spectrum(args.spectrum)
    old_rows = old_witnesses(args.old_cover)
    rows = []
    for occupation in args.occupations:
        row = optimize_occupation(envelope, spectrum, old_rows, occupation)
        rows.append(row)
        print(
            f"occupation={occupation},margin_bits={row['margin_bits']:.9f},"
            f"theta={row['value_probability']:.9f}",
            flush=True,
        )

    result = {
        "schema": "bch250-124-shell-sensitive-dense-diagnostic-v1",
        "status": "BINARY64_SAMPLED_DIAGNOSTIC",
        "parameters": {
            "outer_bits": B,
            "outer_dimension": 124,
            "outer_rows": L,
            "output_bits": N,
            "distance": DISTANCE,
            "inner_epochs": EPOCHS,
            "fanout_layers": args.fanout_layers,
            "spectrum": str(args.spectrum),
            "selection": str(args.selection),
            "occupations": args.occupations,
        },
        "worst_sampled": min(rows, key=lambda row: row["margin_bits"]),
        "all_sampled_below_2^-40": all(row["margin_bits"] > 40.0 for row in rows),
        "rows": rows,
        "limitations": [
            "The occupation set is sampled and is not an interval cover.",
            "The spectrum and optimization use nearest binary64 arithmetic.",
            "An outward verifier must replay fixed witnesses and shell bounds.",
        ],
    }
    args.output.write_text(
        json.dumps(result, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    print(f"worst_margin_bits={result['worst_sampled']['margin_bits']:.9f}")
    print(f"output={args.output}")


if __name__ == "__main__":
    main()
