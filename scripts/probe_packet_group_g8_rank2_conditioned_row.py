#!/usr/bin/env python3
"""Layout-unbound rank-2 algebra prototype; invalid for the frozen layout.

DO NOT RUN THIS AS A FROZEN-CONSTRUCTION DIAGNOSTIC.  In the sloped layout,
a block lies in one common tile of all three bands only in lane l=0.  There is
no second aligned complete row.  The parity and rank-one calculations below
describe an abstract aligned-row model, not the committed construction.
Moreover, the puncture/graph path would require a proved 2-factor of the
three-partite tile hypergraph containing every punctured block.  No such
certificate exists.

Condition two BCH rows ``x`` and ``y`` in each 64-row tile.  At one physical
coordinate, their joint type is one of ``00,01,10,11``.  The packet moment
depends only on the joint Hamming weight, so the required local log factors
are

    ell[b,u] = p_b log E_w R_{w+u}(t)^(1/p_b),  u=0,1,2,

where ``w`` is Binomial(62,1/2).  Thus ``00`` uses ``u=0``, ``01`` and ``10``
use ``u=1``, and ``11`` uses ``u=2``.

The existing exact split spectra support a sound parity collapse.  Write
``d=x+y``.  On ``d=0`` coordinates, replace the two equal types by
``max(ell[b,0],ell[b,2])``.  On ``d=1`` coordinates, use ``ell[b,1]``.
The sum then depends only on the one codeword ``d``; the other codeword gives
an exact factor 2^64.  Existing band-(0,1), punctured band-(0,1), band-(1,2),
and graph spectra therefore apply unchanged.  This is exactly the sound
``conditioned_rows=2`` specialization of the canonical evaluator.

A second sound reuse path applies the entrywise rank-one envelope

    log G[x,y] <= a + d*x + a + d*y,

where ``a=max(ell[0]/2,ell[2]/2-d,(ell[1]-d)/2)`` in each band.  The two
codeword sums then factor into two existing one-row pair-Cauchy enumerators.
At a graph hole, one enumerator is punctured, one is normal, and
``exp(a_0+d_0*g)`` is the graph-bit factor.

Avoiding the parity maximum requires second-order split spectra for ordered
pairs of projected BCH codewords.  The current one-codeword split spectra do
not determine those tables.  This probe reports the exact missing artifact
schema and makes no full-joint rank-2 claim.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import math
import os
import tempfile
from pathlib import Path
from typing import Sequence

import numpy as np
from scipy.optimize import minimize
from scipy.special import logsumexp

from certify_three_band_exact_length import load_graph_spectrum
from packet_group_outer_profile import K, atom_count
from probe_packet_group_conditioned_row_outer import (
    conditioned_norms,
    evaluate_conditioned_outer,
    load_split_spectrum,
    pair_cauchy_log_enumerator,
)


ROOT = Path(__file__).resolve().parents[1]
SPECTRUM01 = ROOT / "out/ebch85_band01_split_spectrum.csv"
PUNCTURED01 = ROOT / "out/ebch84_punctured_band01_split_spectrum.csv"
SPECTRUM12 = ROOT / "out/ebch86_band12_split_spectrum.csv"
GRAPH_SPECTRUM = ROOT / "scripts/ebch128_graph24_spectrum.csv"
CANONICAL_EVALUATOR = ROOT / "scripts/probe_packet_group_conditioned_row_outer.py"
BAND_SIZES = (42, 43, 43)
JOINT_TYPES = ("00", "01", "10", "11")
TYPE_SHIFT = {"00": 0, "01": 1, "10": 1, "11": 2}
PINNED_SHA256 = {
    SPECTRUM01: "10315518bf5aa72199a48bb201e9f11d159c61bd9a83aa6f0fc1948887cf1f06",
    PUNCTURED01: "78fd4730ce818fabf923c237f88115e7c9ad4ecc6588f347c8e653ff93839ce3",
    SPECTRUM12: "d7bf4816ab602a03bfdb526cabfa437ab412e5894e6ea64987ae736752f448e1",
    GRAPH_SPECTRUM: "79a3f8f34280996d46ccd076c515ef1ccdc3b1f80f389ddaf2e7a24286c4f2b3",
    CANONICAL_EVALUATOR: "c6ac971e857aa2f043ad0c61175cae151eab3d54c14442d100627760cd2d3774",
}


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1 << 20), b""):
            digest.update(chunk)
    return digest.hexdigest()


def atomic_json(path: Path, payload: dict) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    descriptor, temporary_name = tempfile.mkstemp(
        prefix=path.name + ".", suffix=".tmp", dir=path.parent
    )
    try:
        with os.fdopen(descriptor, "w", encoding="utf-8", newline="\n") as handle:
            json.dump(payload, handle, sort_keys=True, separators=(",", ":"))
            handle.write("\n")
            handle.flush()
            os.fsync(handle.fileno())
        os.replace(temporary_name, path)
    except BaseException:
        try:
            os.unlink(temporary_name)
        except FileNotFoundError:
            pass
        raise


def validate_band_coefficients(values: Sequence[float]) -> tuple[float, float, float]:
    if len(values) != 3:
        raise ValueError("rank-2 probe requires three band coefficients")
    coefficients = tuple(float(value) for value in values)
    if any(not math.isfinite(value) or not 0.0 < value <= 1.0 for value in coefficients):
        raise ValueError("band coefficients must be finite and lie in (0,1]")
    for left, right in ((0, 1), (0, 2), (1, 2)):
        if coefficients[left] + coefficients[right] < 1.0:
            raise ValueError("band coefficients violate the linear-BL dimension condition")
    return coefficients


def rank2_joint_log_factors(
    log_variables: np.ndarray,
    band_coefficients: Sequence[float],
    group_bits: int = 8,
) -> list[dict[str, float]]:
    """Return the four joint-type log factors in each band."""

    coefficients = validate_band_coefficients(band_coefficients)
    if log_variables.shape != (group_bits + 1,) or log_variables[0] != 0.0:
        raise ValueError("rank-2 probe received malformed packet log variables")
    if not np.all(np.isfinite(log_variables)):
        raise ValueError("rank-2 packet log variables must be finite")
    result = []
    for coefficient in coefficients:
        shifts = conditioned_norms(
            log_variables, coefficient, conditioned_rows=2, group_bits=group_bits
        )
        result.append({joint: shifts[TYPE_SHIFT[joint]] for joint in JOINT_TYPES})
    return result


def parity_collapse(joint_factors: Sequence[dict[str, float]]) -> dict:
    if len(joint_factors) != 3:
        raise ValueError("rank-2 parity collapse requires three band maps")
    equal = []
    different = []
    hole0 = []
    hole1 = []
    for band in joint_factors:
        if set(band) != set(JOINT_TYPES):
            raise ValueError("rank-2 joint factor map has the wrong types")
        if band["01"] != band["10"]:
            raise ValueError("packet moment must be symmetric between types 01 and 10")
        equal.append(max(band["00"], band["11"]))
        different.append(band["01"])
        # The graph bit replaces the first conditioned bit at the puncture;
        # maximize over the second conditioned bit.
        hole0.append(max(band["00"], band["01"]))
        hole1.append(max(band["10"], band["11"]))
    return {
        "equal_parity_log_factors": equal,
        "different_parity_log_factors": different,
        "parity_log_ratios": [right - left for left, right in zip(equal, different)],
        "puncture_hole_log_factor_given_graph_0": hole0,
        "puncture_hole_log_factor_given_graph_1": hole1,
    }


def rank_one_envelope(
    joint_factors: Sequence[dict[str, float]], d_values: Sequence[float]
) -> dict:
    """Construct G_xy <= exp(a+d*x) exp(a+d*y) in every band."""

    if len(joint_factors) != 3 or len(d_values) != 3:
        raise ValueError("rank-one envelope requires three bands and three d values")
    rows = []
    for band, (joint, raw_d) in enumerate(zip(joint_factors, d_values)):
        d = float(raw_d)
        if not math.isfinite(d):
            raise ValueError("rank-one d values must be finite")
        l0, l1, l2 = joint["00"], joint["01"], joint["11"]
        candidates = (l0 / 2.0, l2 / 2.0 - d, (l1 - d) / 2.0)
        a = max(candidates)
        slacks = {
            "00": 2.0 * a - l0,
            "01": 2.0 * a + d - l1,
            "10": 2.0 * a + d - l1,
            "11": 2.0 * a + 2.0 * d - l2,
        }
        if min(slacks.values()) < -1e-12:
            raise AssertionError("rank-one envelope failed an entrywise inequality")
        rows.append(
            {
                "band": band,
                "d": d,
                "a": a,
                "active_lower_bounds_for_a": {
                    "l0_over_2": candidates[0],
                    "l2_over_2_minus_d": candidates[1],
                    "l1_minus_d_over_2": candidates[2],
                },
                "entrywise_log_slacks": slacks,
            }
        )
    return {
        "bands": rows,
        "a": [row["a"] for row in rows],
        "d": [row["d"] for row in rows],
        "status": "SOUND_ENTRYWISE_RANK_ONE_ENVELOPE",
    }


def evaluate_rank_one_factorized(
    profile: list[int],
    log_variables: np.ndarray,
    band_coefficients: Sequence[float],
    theta: float,
    spectra: tuple[np.ndarray, np.ndarray, np.ndarray],
    d_values: Sequence[float],
) -> tuple[float, dict]:
    """Evaluate the sound factorized rank-2 envelope for frozen d values."""

    coefficients = validate_band_coefficients(band_coefficients)
    if len(profile) != 9 or sum(profile) != atom_count(8):
        raise ValueError("rank-one rank-2 probe received a malformed g=8 profile")
    if not 0.0 <= theta <= 1.0:
        raise ValueError("rank-one rank-2 probe received an invalid Cauchy split")
    spectrum01, punctured01, spectrum12 = spectra
    joint = rank2_joint_log_factors(log_variables, coefficients)
    envelope = rank_one_envelope(joint, d_values)
    a = envelope["a"]
    d = envelope["d"]
    normal_enumerator = pair_cauchy_log_enumerator(
        spectrum01, spectrum12, tuple(d), theta
    )
    punctured_enumerator = pair_cauchy_log_enumerator(
        punctured01, spectrum12, tuple(d), theta
    )
    normal_row_log = sum(size * value for size, value in zip(BAND_SIZES, a)) + normal_enumerator
    punctured_row_log = (
        (BAND_SIZES[0] - 1) * a[0]
        + BAND_SIZES[1] * a[1]
        + BAND_SIZES[2] * a[2]
        + punctured_enumerator
    )
    free_rows_log = 62 * 64 * math.log(2.0)
    normal_tile_log = free_rows_log + 2.0 * normal_row_log
    hole_common_log = free_rows_log + normal_row_log + punctured_row_log
    graph_terms = [
        math.log(count)
        - 24 * math.log(2.0)
        + (128 - weight) * (hole_common_log + a[0])
        + weight * (hole_common_log + a[0] + d[0])
        for weight, count in enumerate(load_graph_spectrum())
        if count
    ]
    graph_average_log = float(logsumexp(graph_terms))
    natural = (
        128 * normal_tile_log
        + graph_average_log
        - float(np.asarray(profile, dtype=np.float64) @ log_variables)
    )
    return natural / math.log(2.0), {
        "rank_one_envelope": envelope,
        "normal_one_row_pair_cauchy_enumerator_log": normal_enumerator,
        "punctured_one_row_pair_cauchy_enumerator_log": punctured_enumerator,
        "normal_one_row_log": normal_row_log,
        "punctured_one_row_log": punctured_row_log,
        "free_62_rows_log": free_rows_log,
        "normal_tile_log": normal_tile_log,
        "hole_common_before_graph_factor_log": hole_common_log,
        "graph_bit_0_log_factor": a[0],
        "graph_bit_1_log_factor": a[0] + d[0],
        "graph_average_log": graph_average_log,
        "factorization": (
            "normal tile = 62*64*ln2 + 2*(sum_b |B_b| a_b + one-row enumerator); "
            "hole = one normal-row enumerator + one punctured-row enumerator + graph factor"
        ),
        "status": "SOUND_FIXED_D_RANK2_RANK_ONE_ENVELOPE",
    }


def optimize_rank_one_d(
    profile: list[int],
    log_variables: np.ndarray,
    band_coefficients: Sequence[float],
    theta: float,
    spectra: tuple[np.ndarray, np.ndarray, np.ndarray],
) -> tuple[list[float], dict]:
    """Bounded binary64 discovery optimization of the three envelope slopes."""

    joint = rank2_joint_log_factors(log_variables, band_coefficients)
    balanced = np.asarray(
        [(band["11"] - band["00"]) / 2.0 for band in joint], dtype=np.float64
    )

    def objective(point: np.ndarray) -> float:
        value, _details = evaluate_rank_one_factorized(
            profile,
            log_variables,
            band_coefficients,
            theta,
            spectra,
            point,
        )
        return value

    candidates = []
    for start in (np.zeros(3), balanced):
        result = minimize(
            objective,
            start,
            method="Powell",
            bounds=[(-160.0, 160.0)] * 3,
            options={"maxiter": 1000, "xtol": 1e-9, "ftol": 1e-11},
        )
        candidates.append((objective(result.x), result))
    _value, result = min(candidates, key=lambda row: row[0])
    return result.x.tolist(), {
        "success": bool(result.success),
        "message": str(result.message),
        "iterations": int(result.nit),
        "evaluations": int(result.nfev),
        "status": "DISCOVERY_BINARY64_FIXED_FORMULA_NO_OUTWARD_CLAIM",
    }


def missing_second_order_certificate() -> dict:
    """Describe the exact artifacts needed to avoid parity maximization."""

    normal_fields = [
        "left_n01",
        "left_n10",
        "left_n11",
        "right_n01",
        "right_n10",
        "right_n11",
        "count",
    ]
    return {
        "status": "MISSING_REQUIRED_SECOND_ORDER_SPLIT_CERTIFICATE",
        "reason": (
            "One-codeword split weights do not determine joint types of an "
            "ordered pair (x,y), equivalently the bandwise weights of x,y,x+y."
        ),
        "normal_tables": [
            {
                "name": "ebch_rank2_band01_joint_split",
                "band_sizes": [42, 43],
                "fields": normal_fields,
                "row_semantics": (
                    "count ordered projected-codeword pairs (x,y) by the three "
                    "non-00 joint-type counts in each band; n00 is the band size "
                    "minus n01+n10+n11"
                ),
                "required_mass": "2^128",
            },
            {
                "name": "ebch_rank2_band12_joint_split",
                "band_sizes": [43, 43],
                "fields": normal_fields,
                "row_semantics": (
                    "same ordered-pair joint-type counts for bands one and two"
                ),
                "required_mass": "2^128",
            },
        ],
        "punctured_transform": {
            "source": "ebch_rank2_band01_joint_split",
            "puncture_choices": 42,
            "fields": [
                "deleted_type",
                "left_remaining_n01",
                "left_remaining_n10",
                "left_remaining_n11",
                "right_n01",
                "right_n10",
                "right_n11",
                "pair_count",
            ],
            "deleted_type_domain": list(JOINT_TYPES),
            "required_mass": "42*2^128",
            "reason_deleted_type_is_required": (
                "after the graph bit replaces the first deleted bit, the hole "
                "factor still depends on the second deleted bit"
            ),
        },
        "required_certificate_checks": [
            "exact integer nonnegative counts",
            "declared total masses",
            "one all-zero ordered pair",
            "symmetry under swapping x and y (n01 <-> n10)",
            "punctured table equals the exact 42-coordinate deletion transform",
            "bandwise marginals agree with the committed one-codeword split spectra",
            "content digest and generator/projection-code closure",
        ],
        "reusable_existing_inputs": {
            "graph_weight_spectrum": (
                "still sufficient after the punctured joint table produces one "
                "hole factor for graph bit 0 and one for graph bit 1"
            ),
            "pair_cauchy_split": (
                "still combines exact rank-2 band01 and band12 enumerators"
            ),
        },
    }


def puncture_joint_cell(
    left_counts: dict[str, int], right_counts: dict[str, int], count: int
) -> list[dict]:
    """Exact deletion transform for one unpunctured joint-spectrum cell."""

    if set(left_counts) != set(JOINT_TYPES) or set(right_counts) != set(JOINT_TYPES):
        raise ValueError("joint cell must contain all four types in both bands")
    if sum(left_counts.values()) != 42 or sum(right_counts.values()) != 43:
        raise ValueError("joint cell has incorrect band sizes")
    if count < 0 or any(value < 0 for value in (*left_counts.values(), *right_counts.values())):
        raise ValueError("joint cell counts must be nonnegative")
    rows = []
    for deleted_type in JOINT_TYPES:
        multiplicity = left_counts[deleted_type]
        if not multiplicity or not count:
            continue
        remaining = dict(left_counts)
        remaining[deleted_type] -= 1
        rows.append(
            {
                "deleted_type": deleted_type,
                "left_remaining": remaining,
                "right": dict(right_counts),
                "pair_count": count * multiplicity,
            }
        )
    return rows


def load_inputs() -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    for path, expected in PINNED_SHA256.items():
        observed = sha256_file(path)
        if observed != expected:
            raise ValueError(
                f"rank-2 frozen input digest mismatch for {path}: "
                f"expected {expected}, observed {observed}"
            )
    return (
        load_split_spectrum(SPECTRUM01),
        load_split_spectrum(PUNCTURED01, count_field="pair_count", divisor=42),
        load_split_spectrum(SPECTRUM12),
    )


def evaluate_probe(
    profile: list[int],
    log_variables: np.ndarray,
    band_coefficients: Sequence[float],
    theta: float,
    spectra: tuple[np.ndarray, np.ndarray, np.ndarray],
    rank_one_d: Sequence[float] = (0.0, 0.0, 0.0),
    optimize_d: bool = False,
) -> dict:
    coefficients = validate_band_coefficients(band_coefficients)
    spectrum01, punctured01, spectrum12 = spectra
    outer_log2, canonical = evaluate_conditioned_outer(
        profile,
        log_variables,
        coefficients[1],
        theta,
        spectrum01,
        punctured01,
        spectrum12,
        conditioned_rows=2,
        band_coefficients=coefficients,
        group_bits=8,
    )
    joint = rank2_joint_log_factors(log_variables, coefficients)
    collapsed = parity_collapse(joint)

    canonical_joint = canonical["conditioned_log_factors"]
    expected_joint = [
        [band["00"], band["01"], band["11"]] for band in joint
    ]
    if not np.array_equal(np.asarray(canonical_joint), np.asarray(expected_joint)):
        raise AssertionError("rank-2 joint factors disagree with canonical evaluator")
    expected_parity = [
        [collapsed["equal_parity_log_factors"][band], collapsed["different_parity_log_factors"][band]]
        for band in range(3)
    ]
    if not np.array_equal(
        np.asarray(canonical["conditioned_parity_log_factors"]),
        np.asarray(expected_parity),
    ):
        raise AssertionError("rank-2 parity collapse disagrees with canonical evaluator")
    punctured_common = canonical["hole0_tile_log"] - collapsed[
        "puncture_hole_log_factor_given_graph_0"
    ][0]
    if canonical["hole1_tile_log"] != (
        punctured_common
        + collapsed["puncture_hole_log_factor_given_graph_1"][0]
    ):
        raise AssertionError("rank-2 graph-hole collapse disagrees with canonical evaluator")

    optimizer = None
    chosen_d = list(rank_one_d)
    if optimize_d:
        chosen_d, optimizer = optimize_rank_one_d(
            profile, log_variables, coefficients, theta, spectra
        )
    rank_one_log2, rank_one_details = evaluate_rank_one_factorized(
        profile,
        log_variables,
        coefficients,
        theta,
        spectra,
        chosen_d,
    )

    return {
        "schema": "permute-conv.packet-group-g8-rank2-conditioned-row-probe.v1",
        "status": "LAYOUT_UNBOUND_INVALID_FOR_FROZEN_SLOPED_CONSTRUCTION_DO_NOT_RUN",
        "parity_collapsed_outer_log2": outer_log2,
        "rank_one_factorized_outer_log2": rank_one_log2,
        "best_sound_rank2_outer_log2": min(outer_log2, rank_one_log2),
        "rank_one_saving_vs_parity_collapse_bits": outer_log2 - rank_one_log2,
        "rank_one_d_optimizer": optimizer,
        "profile": profile,
        "log_variables": log_variables.tolist(),
        "band_coefficients": list(coefficients),
        "pair_cauchy_theta": theta,
        "joint_type_shift": TYPE_SHIFT,
        "required_packet_moment_variables": {
            "base": "R_w(t), w=0,...,64",
            "rank2_local": "ell[b,u], b=0,1,2 and u=0,1,2",
            "joint_mapping": {"00": "u=0", "01": "u=1", "10": "u=1", "11": "u=2"},
        },
        "joint_log_factors": joint,
        "parity_collapse": collapsed,
        "rank_one_factorized": rank_one_details,
        "canonical_details": canonical,
        "missing_full_joint_certificate": missing_second_order_certificate(),
        "source_bindings": {
            "spectrum01": {"path": str(SPECTRUM01), "sha256": sha256_file(SPECTRUM01)},
            "punctured01": {"path": str(PUNCTURED01), "sha256": sha256_file(PUNCTURED01)},
            "spectrum12": {"path": str(SPECTRUM12), "sha256": sha256_file(SPECTRUM12)},
            "graph_spectrum": {"path": str(GRAPH_SPECTRUM), "sha256": sha256_file(GRAPH_SPECTRUM)},
            "canonical_evaluator": {
                "path": str(CANONICAL_EVALUATOR),
                "sha256": sha256_file(CANONICAL_EVALUATOR),
            },
        },
        "scope_limit": (
            "Abstract aligned-row algebra only. The frozen sloped layout has no "
            "second aligned row, and the puncture path lacks a required hypergraph "
            "2-factor certificate. This output is invalid for the construction."
        ),
    }


def static_self_test() -> None:
    if TYPE_SHIFT != {"00": 0, "01": 1, "10": 1, "11": 2}:
        raise AssertionError("rank-2 joint type map changed")
    spectra = load_inputs()
    profile = [atom_count(8)] + [0] * 8
    zero = evaluate_probe(
        profile,
        np.zeros(9, dtype=np.float64),
        (0.4, 0.6, 0.6),
        0.5,
        spectra,
        rank_one_d=(0.0, 0.0, 0.0),
    )
    for key in ("parity_collapsed_outer_log2", "rank_one_factorized_outer_log2"):
        if not math.isclose(zero[key], K, rel_tol=0.0, abs_tol=1e-6):
            raise AssertionError(f"rank-2 zero-tilt normalization failed for {key}")

    tilted = evaluate_probe(
        profile,
        np.asarray([0.0, -0.2, 0.1, -0.4, 0.3, -0.1, -0.5, 0.2, -0.3]),
        (0.43, 0.57, 0.61),
        0.37,
        spectra,
        rank_one_d=(0.2, -0.1, 0.3),
    )
    if tilted["joint_log_factors"][0]["01"] != tilted["joint_log_factors"][0]["10"]:
        raise AssertionError("rank-2 01/10 packet symmetry failed")
    if not math.isfinite(tilted["rank_one_factorized_outer_log2"]):
        raise AssertionError("rank-2 rank-one factorized value is not finite")
    for band in tilted["rank_one_factorized"]["rank_one_envelope"]["bands"]:
        if min(band["entrywise_log_slacks"].values()) < -1e-12:
            raise AssertionError("rank-2 rank-one static envelope has negative slack")

    # Directly check the entrywise factorization on a synthetic 2x2 matrix.
    synthetic_joint = [
        {"00": 0.2, "01": 0.7, "10": 0.7, "11": 1.0},
        {"00": -0.1, "01": 0.4, "10": 0.4, "11": 0.9},
        {"00": 0.0, "01": -0.2, "10": -0.2, "11": 0.3},
    ]
    synthetic_envelope = rank_one_envelope(synthetic_joint, (0.3, -0.2, 0.1))
    for band, row in enumerate(synthetic_envelope["bands"]):
        for left in (0, 1):
            for right in (0, 1):
                joint_type = f"{left}{right}"
                if synthetic_joint[band][joint_type] > (
                    2.0 * row["a"] + row["d"] * (left + right) + 1e-12
                ):
                    raise AssertionError("synthetic rank-one entrywise bound failed")

    left = {"00": 10, "01": 11, "10": 9, "11": 12}
    right = {"00": 13, "01": 8, "10": 10, "11": 12}
    punctured = puncture_joint_cell(left, right, count=7)
    if sum(row["pair_count"] for row in punctured) != 42 * 7:
        raise AssertionError("rank-2 punctured joint transform lost mass")
    for row in punctured:
        if sum(row["left_remaining"].values()) != 41:
            raise AssertionError("rank-2 punctured joint transform has wrong length")

    contract = missing_second_order_certificate()
    if contract["normal_tables"][0]["required_mass"] != "2^128":
        raise AssertionError("rank-2 missing-table contract has wrong mass")
    print("g8 rank-2 conditioned-row static self-test: PASS")
    print(
        f"zero_tilt_parity_log2={zero['parity_collapsed_outer_log2']:.9f} "
        f"rank_one_log2={zero['rank_one_factorized_outer_log2']:.9f}"
    )
    print("joint_types=00:u0,01:u1,10:u1,11:u2")
    print("existing_tables=SOUND_PARITY_COLLAPSE_AND_RANK_ONE_FACTORIZATION")
    print("rank_one_factorization=SOUND_WITH_EXISTING_ONE_ROW_SPLIT_AND_GRAPH_TABLES")
    print("full_joint_status=MISSING_SECOND_ORDER_BAND01_AND_BAND12_SPLIT_CERTIFICATE")
    print("frozen_layout_status=INVALID_NO_SECOND_ALIGNED_ROW_DO_NOT_RUN")


def parse_csv(text: str, cast) -> list:
    return [cast(token.strip()) for token in text.split(",")]


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--self-test", action="store_true")
    parser.add_argument("--profile")
    parser.add_argument("--log-variables")
    parser.add_argument("--band-coefficients", default="0.4,0.6,0.6")
    parser.add_argument("--rank-one-d", default="0,0,0")
    parser.add_argument("--optimize-d", action="store_true")
    parser.add_argument("--theta", type=float, default=0.5)
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    if args.self_test:
        static_self_test()
        return
    parser.error(
        "disabled: rank-2 aligned-row model is invalid for the frozen sloped layout"
    )
    if args.profile is None or args.log_variables is None:
        parser.error("--profile and --log-variables are required outside --self-test")
    profile = parse_csv(args.profile, int)
    variables = np.asarray(parse_csv(args.log_variables, float), dtype=np.float64)
    coefficients = parse_csv(args.band_coefficients, float)
    rank_one_d = parse_csv(args.rank_one_d, float)
    if len(profile) != 9 or sum(profile) != atom_count(8):
        parser.error("--profile must contain nine counts summing to the g=8 atom count")
    if variables.shape != (9,):
        parser.error("--log-variables must contain nine values")
    if len(rank_one_d) != 3:
        parser.error("--rank-one-d must contain three values")
    try:
        report = evaluate_probe(
            profile,
            variables,
            coefficients,
            args.theta,
            load_inputs(),
            rank_one_d=rank_one_d,
            optimize_d=args.optimize_d,
        )
    except ValueError as error:
        parser.error(str(error))
    rendered = json.dumps(report, indent=2, sort_keys=True)
    print(rendered)
    if args.output is not None:
        atomic_json(args.output, report)


if __name__ == "__main__":
    main()
