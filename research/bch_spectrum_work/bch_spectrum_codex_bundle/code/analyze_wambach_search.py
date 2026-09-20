#!/usr/bin/env python3
"""Summarize the calibrated systematic-basis BCH search.

Rigorous fields use only disjoint affine-orbit closures.  The Chao1 fields are
explicitly heuristic because row-combination discoveries are not IID samples.
"""

from __future__ import annotations

import json
import math
import csv
from collections import Counter
from pathlib import Path

from bch_quotient import binary_poly_divmod, generator_polynomial


ROOT = Path(__file__).resolve().parents[1]
GENERATED = ROOT / "generated"


def load_shell(filename: str, weight: int) -> dict:
    payload = json.loads((GENERATED / filename).read_text(encoding="utf-8"))
    return payload["shells"][str(weight)]


def expected_full_rank(n: int, k: int, weight: int) -> float:
    return (2**k - 1) * math.comb(n, weight) / (2**n - 1)


def logsumexp(values: list[float]) -> float:
    maximum = max(values)
    return maximum + math.log(sum(math.exp(value - maximum) for value in values))


def logaddexp(left: float, right: float) -> float:
    maximum = max(left, right)
    return maximum + math.log(
        math.exp(left - maximum) + math.exp(right - maximum)
    )


def incidence_summary(shell: dict) -> dict[str, float | int | dict[str, int]]:
    frequencies = Counter(
        int(record["search_hits_in_orbit"]) for record in shell["orbits"]
    )
    observed = len(shell["orbits"])
    singletons = frequencies[1]
    doubletons = frequencies[2]
    chao1 = (
        observed + singletons * singletons / (2 * doubletons)
        if doubletons
        else math.inf
    )
    rigorous_mass = int(shell["rigorous_lower_bound"])
    heuristic_mass = rigorous_mass * chao1 / observed
    return {
        "observed_orbits": observed,
        "singletons": singletons,
        "doubletons": doubletons,
        "hit_multiplicity_histogram": {
            str(key): value for key, value in sorted(frequencies.items())
        },
        "rigorous_orbit_closure_lower_bound": rigorous_mass,
        "chao1_orbit_count_heuristic": chao1,
        "chao1_shell_mass_heuristic": heuristic_mass,
    }


def load_tsv_orbits(filename: str) -> list[dict[str, int]]:
    with (GENERATED / filename).open(newline="", encoding="ascii") as stream:
        records = []
        for row in csv.DictReader(stream, delimiter="\t"):
            word = sum(int(row[f"limb{index}"], 16) << (64 * index) for index in range(4))
            records.append(
                {
                    "word": word,
                    "orbit_size": int(row["orbit_size"]),
                    "stabilizer_size": int(row["stabilizer_size"]),
                    "search_hits": int(row["search_hits"]),
                }
            )
    return records


def tsv_incidence_summary(records: list[dict[str, int]]) -> dict:
    frequencies = Counter(record["search_hits"] for record in records)
    observed = len(records)
    singletons = frequencies[1]
    doubletons = frequencies[2]
    chao1 = observed + singletons * singletons / (2 * doubletons)
    rigorous_mass = sum(record["orbit_size"] for record in records)
    return {
        "observed_orbits": observed,
        "singletons": singletons,
        "doubletons": doubletons,
        "rigorous_orbit_closure_lower_bound": rigorous_mass,
        "chao1_orbit_count_heuristic": chao1,
        "chao1_shell_mass_heuristic": rigorous_mass * chao1 / observed,
    }


def main() -> None:
    exact_calibration = [
        ("B8", 4, 14, "wambach_b8_orbits_w4.json"),
        ("B32", 8, 620, "wambach_b32_orbits_w8.json"),
        ("B128", 22, 243840, "wambach_b128_orbits_w22.json"),
    ]
    calibration = []
    for name, weight, exact_count, filename in exact_calibration:
        shell = load_shell(filename, weight)
        recovered = int(shell["rigorous_lower_bound"])
        calibration.append(
            {
                "code": name,
                "minimum_weight": weight,
                "exact_minimum_shell": exact_count,
                "recovered_by_radius4_affine_closure": recovered,
                "recall": recovered / exact_count,
            }
        )

    p38_r4 = load_shell("wambach_p_orbits_w38.json", 38)
    p38_r5 = load_shell("wambach_p_orbits_r5_w38.json", 38)
    q40_r4 = load_shell("wambach_q_orbits_w40.json", 40)
    q40_r5 = load_shell("wambach_q_orbits_r5_w40.json", 40)
    p40_r5 = load_shell("wambach_p_orbits_r5_w40.json", 40)

    p38 = incidence_summary(p38_r5)
    q40 = incidence_summary(q40_r5)
    p40 = incidence_summary(p40_r5)

    c38_lower = 31 * int(p38_r5["rigorous_lower_bound"]) // 255
    c38_heuristic = 31 * float(p38["chao1_shell_mass_heuristic"]) / 255

    q_generator = generator_polynomial(39)
    p40_outside_q = 0
    for record in p40_r5["orbits"]:
        punctured = int(record["canonical_hex"], 16) & ((1 << 255) - 1)
        if binary_poly_divmod(punctured, q_generator)[1] != 0:
            p40_outside_q += int(record["orbit_size"])
    c40_lower = int(q40_r5["rigorous_lower_bound"]) + 31 * p40_outside_q // 255
    c40_heuristic = float(q40["chao1_shell_mass_heuristic"]) + 31 * (
        float(p40["chao1_shell_mass_heuristic"])
        - float(q40["chao1_shell_mass_heuristic"])
    ) / 255

    random_c38 = expected_full_rank(256, 128, 38)
    random_c40 = expected_full_rank(256, 128, 40)
    diagnostic = json.loads(
        (GENERATED / "random_inner_oa15_majorant_diagnostic.json").read_text(
            encoding="utf-8"
        )
    )
    application = diagnostic["proofish_low_prefix_reduction"]
    caps = {
        int(record["weight"]): float(record["sufficient_log2_multiplicity_cap"])
        for record in application["sufficient_common_multiplier_shell_caps"]
    }


    def additional_shell(weight: int, radius: int = 5) -> dict:
        p_records = load_tsv_orbits(
            f"wambach_p_orbits_r{radius}_w{weight}.tsv"
        )
        q_records = load_tsv_orbits(
            f"wambach_q_orbits_r{radius}_w{weight}.tsv"
        )
        p_summary = tsv_incidence_summary(p_records)
        q_summary = tsv_incidence_summary(q_records)
        p_outside_q = sum(
            record["orbit_size"]
            for record in p_records
            if binary_poly_divmod(record["word"] & ((1 << 255) - 1), q_generator)[1]
            != 0
        )
        c_lower = q_summary["rigorous_orbit_closure_lower_bound"] + 31 * p_outside_q // 255
        c_heuristic = q_summary["chao1_shell_mass_heuristic"] + 31 * (
            p_summary["chao1_shell_mass_heuristic"]
            - q_summary["chao1_shell_mass_heuristic"]
        ) / 255
        random_c = expected_full_rank(256, 128, weight)
        cap = 2**caps[weight]
        return {
            "weight": weight,
            "search_radius": radius,
            "p_search": p_summary,
            "q_search": q_summary,
            "p_mass_proved_outside_q": p_outside_q,
            "c_rigorous_union_lower_bound": c_lower,
            "c_chao1_heuristic": c_heuristic,
            "c_random_full_rank_expectation": random_c,
            "c_heuristic_to_random_ratio": c_heuristic / random_c,
            "application_cap_log2": caps[weight],
            "application_cap_to_heuristic_ratio": cap / c_heuristic,
        }

    additional_shells = [
        additional_shell(42),
        additional_shell(44),
        additional_shell(46),
        additional_shell(48),
        additional_shell(50, radius=4),
    ]
    heuristic_counts = {
        38: c38_heuristic,
        40: c40_heuristic,
        **{
            int(record["weight"]): float(record["c_chao1_heuristic"])
            for record in additional_shells
        },
    }
    coefficients = diagnostic["log_shell_coefficients_natural"]
    low_log = logsumexp(
        [
            math.log(count)
            + logaddexp(coefficients[weight], coefficients[256 - weight])
            for weight, count in heuristic_counts.items()
        ]
    )
    high_margin = application["weights_52_and_above_margin_bits"]
    high_log = -high_margin * math.log(2.0)
    combined_log = logaddexp(low_log, high_log)
    target_log = -40.0 * math.log(2.0)
    low_multiplier = (
        math.exp(target_log) - math.exp(high_log)
    ) / math.exp(low_log)

    report = {
        "classification": {
            "orbit_closure_counts": "rigorous lower bounds",
            "radius4_calibration": "exact comparison with published spectra",
            "chao1": "heuristic only; discoveries are not IID",
            "completeness_claim_for_256": False,
        },
        "radius4_minimum_shell_calibration": calibration,
        "p_weight38": {
            "radius4_orbits": len(p38_r4["orbits"]),
            "radius4_lower_bound": p38_r4["rigorous_lower_bound"],
            "radius5": p38,
            "random_full_rank_expectation": expected_full_rank(256, 131, 38),
        },
        "q_weight40": {
            "radius4_orbits": len(q40_r4["orbits"]),
            "radius4_lower_bound": q40_r4["rigorous_lower_bound"],
            "radius5": q40,
            "random_full_rank_expectation": expected_full_rank(256, 123, 40),
        },
        "p_weight40": {
            "radius5": p40,
            "radius5_mass_proved_outside_q": p40_outside_q,
            "random_full_rank_expectation": expected_full_rank(256, 131, 40),
        },
        "intermediate_c": {
            "weight38_rigorous_lower_bound": c38_lower,
            "weight38_chao1_heuristic": c38_heuristic,
            "weight38_random_full_rank_expectation": random_c38,
            "weight38_heuristic_to_random_ratio": c38_heuristic / random_c38,
            "weight38_application_cap_log2": caps[38],
            "weight38_application_cap_to_heuristic_ratio": 2**caps[38]
            / c38_heuristic,
            "weight40_rigorous_union_lower_bound": c40_lower,
            "weight40_chao1_heuristic": c40_heuristic,
            "weight40_random_full_rank_expectation": random_c40,
            "weight40_heuristic_to_random_ratio": c40_heuristic / random_c40,
            "weight40_application_cap_log2": caps[40],
            "weight40_application_cap_to_heuristic_ratio": 2**caps[40]
            / c40_heuristic,
            "joint_low_shell_multiplier_before_40_bit_failure": application[
                "joint_low_prefix_multiplier_before_40_bit_failure"
            ],
        },
        "additional_intermediate_c_shells": additional_shells,
        "complete_functional_heuristic": {
            "classification": (
                "Chao1 low-shell heuristic combined with exact combinatorial "
                "upper bounds for every shell of weight at least 52; binary64 "
                "RandomStepConv transfer; not a certificate"
            ),
            "low_shells": sorted(heuristic_counts),
            "low_shell_only_margin_bits": -low_log / math.log(2.0),
            "weights_52_and_above_rigorous_margin_bits": high_margin,
            "combined_margin_bits": -combined_log / math.log(2.0),
            "joint_multiplier_on_heuristic_low_shells_before_40_bits": low_multiplier,
            "joint_multiplier_log2": math.log2(low_multiplier),
        },
    }

    output = GENERATED / "wambach_search_evidence.json"
    output.write_text(json.dumps(report, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(report["intermediate_c"], indent=2))
    print(json.dumps(report["additional_intermediate_c_shells"], indent=2))
    print(json.dumps(report["complete_functional_heuristic"], indent=2))
    print(f"wrote {output}")


if __name__ == "__main__":
    main()
