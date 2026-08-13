#!/usr/bin/env python3
"""Outward probe for exact graph-spectrum/puncture-averaged total-weight outer.

For a data outer word of pre-puncture weight T, 128 graph holes replace 128
randomized data coordinates.  Conditional on the 42-coordinate band-zero
occupancies, Maclaurin bounds the puncture MGF by the 128th power of their
arithmetic mean; ``(1+x)^128 <= exp(128x)`` then separates blocks.  Independent
uniform coordinate bijections make a block's band-zero occupancy at weight w
hypergeometric(128,w,42).  The graph word uses its exact 24-dimensional
subcode spectrum.  Thus, for a frozen 0<q<1, the final-weight coefficient is
bounded by

  q^-W [sum_w A_w q^w E exp((q^-1-1) J/5376)]^16384
       [2^-24 sum_g G_g q^g].

The numerical minimizer only selects q.  The reported interval reevaluates
that frozen binary64 q as an exact rational with guarded outward arithmetic.
"""

from __future__ import annotations

import argparse
import json
import math
from decimal import Decimal
from fractions import Fraction
from pathlib import Path

import numpy as np
from scipy.optimize import minimize_scalar
from scipy.special import gammaln, logsumexp

from certify_packet8_hard_face_drive_inner import exact_float
from certify_packet_group_triangle_ledger import _ln_sum_exp
import certify_packet_group_triangle_ledger as certifier
from certify_three_band_exact_length import load_graph_spectrum
from outward_log2 import LN2, Interval, log2_fraction
from packet_group_outer_profile import K, N, total_weight_outer
from punctured_ebch_outer import full_spectrum


BLOCKS = K // 64
BAND_ZERO = 42
GRAPH_DIMENSION = 24
PUNCTURES = 128
TILES = 256
LANES = 64
PUNCTURE_DENOMINATOR = TILES * LANES * BAND_ZERO // PUNCTURES  # 5376


def ln_fraction(value: Fraction) -> Interval:
    return log2_fraction(value) * LN2


def hypergeometric_rows(weight: int):
    low = max(0, BAND_ZERO - (128 - weight))
    high = min(BAND_ZERO, weight)
    denominator = math.comb(128, BAND_ZERO)
    return tuple(
        (
            selected,
            Fraction(
                math.comb(weight, selected)
                * math.comb(128 - weight, BAND_ZERO - selected),
                denominator,
            ),
        )
        for selected in range(low, high + 1)
    )


def outward_outer(profile: list[int], q: Fraction) -> Interval:
    physical_weight = sum(index * count for index, count in enumerate(profile))
    log_q = ln_fraction(q)
    lam = Fraction(1, PUNCTURE_DENOMINATOR) * (Fraction(1, 1) / q - 1)
    lam_interval = Interval.exact(lam.numerator) / Interval.exact(lam.denominator)
    block_terms = []
    for weight, count in enumerate(full_spectrum()):
        if not count:
            continue
        occupancy = _ln_sum_exp(
            ln_fraction(probability) + lam_interval.times_int(selected)
            for selected, probability in hypergeometric_rows(weight)
        )
        block_terms.append(ln_fraction(Fraction(count)) + log_q.times_int(weight) + occupancy)
    block = _ln_sum_exp(block_terms)

    graph_terms = []
    for weight, count in enumerate(load_graph_spectrum()):
        if count:
            graph_terms.append(
                ln_fraction(Fraction(count, 1 << GRAPH_DIMENSION))
                + log_q.times_int(weight)
            )
    graph = _ln_sum_exp(graph_terms)
    natural = block.times_int(BLOCKS) + graph - log_q.times_int(physical_weight)
    return natural / LN2


def optimize_q(profile: list[int]) -> float:
    physical_weight = sum(index * count for index, count in enumerate(profile))
    spectrum = full_spectrum()
    weights = np.asarray([weight for weight, count in enumerate(spectrum) if count])
    logs = np.log(np.asarray([spectrum[int(weight)] for weight in weights], dtype=np.float64))
    graph_spectrum = load_graph_spectrum()
    graph_weights = np.asarray([weight for weight, count in enumerate(graph_spectrum) if count])
    graph_logs = np.log(
        np.asarray([graph_spectrum[int(weight)] for weight in graph_weights], dtype=np.float64)
    ) - GRAPH_DIMENSION * math.log(2.0)
    hypergeom = {}
    denominator = gammaln(129) - gammaln(BAND_ZERO + 1) - gammaln(128 - BAND_ZERO + 1)
    for weight in weights:
        rows = hypergeometric_rows(int(weight))
        selected = np.asarray([row[0] for row in rows], dtype=np.float64)
        probabilities = np.asarray([float(row[1]) for row in rows])
        hypergeom[int(weight)] = (selected, np.log(probabilities))

    def objective(log_q: float) -> float:
        q = math.exp(log_q)
        lam = (1.0 / q - 1.0) / PUNCTURE_DENOMINATOR
        occupancy = np.asarray(
            [
                logsumexp(hypergeom[int(weight)][1] + hypergeom[int(weight)][0] * lam)
                for weight in weights
            ]
        )
        block = logsumexp(logs + weights * log_q + occupancy)
        graph = logsumexp(graph_logs + graph_weights * log_q)
        return BLOCKS * block + graph - physical_weight * log_q

    result = minimize_scalar(
        objective, bounds=(-4.0, -0.01), method="bounded", options={"xatol": 1e-13}
    )
    return float(result.x)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--profile", required=True, help="comma-separated g=4 profile")
    parser.add_argument("--inner-log2", type=Decimal)
    parser.add_argument(
        "--witness-artifact",
        type=Path,
        help="recompute the selected combined witness inner interval outward",
    )
    parser.add_argument(
        "--witness-reference",
        help="canonical artifact.json:index or unique witness name",
    )
    parser.add_argument("--iterations", type=int, default=40)
    parser.add_argument("--target-log2", type=Decimal, default=Decimal("-111.41506501642425"))
    args = parser.parse_args()
    profile = [int(value) for value in args.profile.split(",")]
    if len(profile) != 5 or sum(profile) != N // 4:
        raise SystemExit("graph-puncture outer: profile must have five counts summing to N/4")
    certifier.GROUP_BITS = 4
    certifier.BLOCK_ATOMS = 16
    certifier._OUTWARD_NORMALIZATION_CACHE.clear()
    log_q = optimize_q(profile)
    q = exact_float(math.exp(log_q))
    outer = outward_outer(profile, q)
    canonical_constant, canonical_charges, canonical_report = (
        certifier._exact_graph_puncture_total_weight_outer(q)
    )
    canonical_outer = canonical_constant
    for count, charge in zip(profile, canonical_charges):
        canonical_outer = canonical_outer - charge.times_int(count)
    if canonical_outer.hi < outer.lo or outer.hi < canonical_outer.lo:
        raise RuntimeError(
            "standalone and canonical graph-puncture outer intervals are disjoint"
        )
    physical_weight = sum(index * count for index, count in enumerate(profile))
    adversarial, _result, _interval = total_weight_outer(physical_weight)
    report = {
        "status": "OUTWARD_GRAPH_SPECTRUM_PUNCTURE_AVERAGED_OUTER_PROBE",
        "profile": profile,
        "physical_weight": physical_weight,
        "q_exact": f"{q.numerator}/{q.denominator}",
        "q_binary64": float(q),
        "outer_log2_interval": [str(outer.lo), str(outer.hi)],
        "adversarial_shift_outer_log2": adversarial,
        "certified_improvement_lower_bits": str(Decimal(str(adversarial)) - outer.hi),
        "lemma": (
            "distinct-tile Maclaurin puncture averaging, independent uniform "
            "coordinate bijections, exact graph24 subcode spectrum"
        ),
        "graph_spectrum_sha256": canonical_report["graph_spectrum_sha256"],
    }
    if args.witness_artifact is not None:
        if not args.witness_reference:
            raise SystemExit("graph-puncture outer: --witness-reference is required")
        if args.iterations <= 0:
            raise SystemExit("graph-puncture outer: --iterations must be positive")
        witnesses = certifier.load_witnesses([args.witness_artifact])
        if args.witness_reference not in witnesses:
            raise SystemExit(
                f"graph-puncture outer: witness {args.witness_reference!r} not found"
            )
        row, artifact = witnesses[args.witness_reference]
        row = dict(row)
        details = dict(row.get("outer_details", {}))
        details.update(
            {
                "pole_exact": f"{q.numerator}/{q.denominator}",
                "log_pole": math.log(float(q)),
            }
        )
        row["outer_type"] = "exact_graph_puncture_total_weight"
        row["outer_details"] = details
        hardened = certifier.harden_witness(
            args.witness_reference,
            row,
            artifact,
            certifier.split_cap_table(),
            args.iterations,
        )
        combined = certifier.evaluate_vertex(tuple(profile), hardened)
        hardened_report = hardened["report"]
        inner_constant_raw = hardened_report["inner_constant_log2_interval"]
        inner_constant = Interval(
            Decimal(inner_constant_raw[0]), Decimal(inner_constant_raw[1])
        )
        inner = inner_constant
        for count, raw_charge in zip(
            profile, hardened_report["inner_charge_log2_intervals"]
        ):
            if count:
                if raw_charge is None:
                    raise ValueError("selected witness is support-ineligible at its profile")
                inner = inner - Interval(
                    Decimal(raw_charge[0]), Decimal(raw_charge[1])
                ).times_int(count)
        inner = inner - certifier.outward_normalization(4, profile)
        recomposed = canonical_outer + inner
        if recomposed.hi < combined.lo or combined.hi < recomposed.lo:
            raise RuntimeError("recomposed outer+inner interval misses hardened combined interval")
        report.update(
            {
                "witness_artifact": str(args.witness_artifact),
                "witness_artifact_sha256": certifier._file_digest(
                    args.witness_artifact.resolve()
                ),
                "witness_reference": args.witness_reference,
                "inner_log2_interval": [str(inner.lo), str(inner.hi)],
                "combined_log2_interval": [str(combined.lo), str(combined.hi)],
                "target_log2": str(args.target_log2),
                "target_margin_lower_bits": str(args.target_log2 - combined.hi),
                "closes_target": combined.hi <= args.target_log2,
                "hardening_report": hardened_report,
            }
        )
    if args.inner_log2 is not None:
        if args.witness_artifact is not None:
            raise SystemExit(
                "graph-puncture outer: choose --inner-log2 or --witness-artifact, not both"
            )
        combined = outer + Interval.exact(args.inner_log2)
        report["inner_log2_input"] = str(args.inner_log2)
        report["combined_log2_interval"] = [str(combined.lo), str(combined.hi)]
        report["target_log2"] = str(args.target_log2)
        report["target_margin_lower_bits"] = str(args.target_log2 - combined.hi)
        report["closes_target"] = combined.hi <= args.target_log2
    print(json.dumps(report, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
