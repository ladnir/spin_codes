"""Run all reproducibility checks shipped with the bundle."""

import json
from fractions import Fraction
from pathlib import Path

from anchors import L71_HALF, U187_DUAL_HALF, U187_LOW_38_60
from spectrum_exact import (
    full_from_symmetric_half,
    macwilliams,
    dual_min_distance,
)
from affine_wambach import (
    P_WORD_OCTAL,
    Q_WORD_OCTAL,
    punctured_weight,
    verify_seed,
)
from bch_quotient import verify_quotient_algebra
from verify_scaled_rational_solution import verify as verify_scaled_lp
from verify_locator_reduction import verify_locator_reduction


def first_nonzero(A, start=1):
    for i in range(start, len(A)):
        if A[i]:
            return i, A[i]
    return None


def main() -> None:
    L = full_from_symmetric_half(L71_HALF)
    assert sum(L) == 1 << 71
    Ld = macwilliams(L, 71)
    assert sum(Ld) == 1 << 185
    assert dual_min_distance(Ld) == 16

    expected_first_dual = {
        16: 927792,
        18: 522240,
        20: 253338624,
        22: 28394188800,
        24: 2816730489600,
        26: 232205133324288,
        28: 16178387961871360,
        30: 962446607116646400,
    }
    for w, a in expected_first_dual.items():
        assert Ld[w] == a

    Ud = full_from_symmetric_half(U187_DUAL_HALF)
    assert sum(Ud) == 1 << 69
    U = macwilliams(Ud, 69)
    assert sum(U) == 1 << 187

    for w, expected in U187_LOW_38_60.items():
        assert U[w] == expected, (w, U[w], expected)

    assert punctured_weight(P_WORD_OCTAL) == 37
    assert punctured_weight(Q_WORD_OCTAL) == 39

    # This is the slowest check but still only 2 * 4 * 65280 affine maps.
    verify_seed(P_WORD_OCTAL, 13, 37)
    verify_seed(Q_WORD_OCTAL, 7, 39)
    quotient = verify_quotient_algebra()
    assert quotient["p_dimension"] == 131
    assert quotient["q_dimension"] == 123
    assert quotient["quotient_dimension"] == 8
    assert quotient["nonzero_shift_orbit"] == 255
    assert quotient["quotient_factor"] == 0x12D

    model_path = Path(__file__).resolve().parents[1] / "generated" / "coupled_lp_exact.json"
    model = json.loads(model_path.read_text())
    constraints = {item["name"]: item for item in model["constraints"]}
    divisibility = {
        item["name"]: item for item in model["divisibility_constraints"]
    }
    assert constraints["Wambach_h_38_lower"]["rhs"] == "256"
    assert constraints["Wambach_q_40_lower"]["rhs"] == "65280"
    assert constraints["orbit_search_h_38_lower"]["rhs"] == "29440"
    assert constraints["orbit_search_q_40_lower"]["rhs"] == "12990720"
    assert constraints["orbit_search_h_40_lower"]["rhs"] == "163072"
    assert constraints["orbit_search_q_42_lower"]["rhs"] == "9204480"
    assert constraints["orbit_search_h_42_lower"]["rhs"] == "1067264"
    assert constraints["orbit_search_q_44_lower"]["rhs"] == "81991680"
    assert constraints["orbit_search_h_44_lower"]["rhs"] == "5977088"
    assert constraints["orbit_search_q_46_lower"]["rhs"] == "381757440"
    assert constraints["orbit_search_h_46_lower"]["rhs"] == "27937536"
    assert constraints["orbit_search_q_48_lower"]["rhs"] == "1980725760"
    assert constraints["orbit_search_h_48_lower"]["rhs"] == "114343168"
    assert constraints["orbit_search_q_50_lower"]["rhs"] == "679042560"
    assert constraints["orbit_search_h_50_lower"]["rhs"] == "28884992"
    assert divisibility["h_38_multiple_128"]["modulus"] == "128"
    assert divisibility["P_shell_2_design_38"]["linear_form"] == {
        "q_38": "703",
        "h_38": "179265",
    }

    johnson_path = model_path.with_name("johnson_n256_w52_d38.json")
    johnson = json.loads(johnson_path.read_text())
    assert johnson["classification"] == "exact rational Johnson-scheme Delsarte LP"
    assert johnson["qsopt_status"] == "OPTIMAL"
    assert johnson["integer_shell_upper"] == 2979928035058718446462766169
    assert all(johnson["independent_exact_primal_checks"].values())
    assert all(johnson["independent_exact_dual_checks"].values())
    johnson54 = json.loads(
        model_path.with_name("johnson_n256_w54_d38.json").read_text()
    )
    assert johnson54["qsopt_status"] == "OPTIMAL"
    assert johnson54["integer_shell_upper"] == 47955805775385501045909339217
    assert all(johnson54["independent_exact_primal_checks"].values())
    assert all(johnson54["independent_exact_dual_checks"].values())
    johnson38 = json.loads(
        model_path.with_name("johnson_n256_w38_d38.json").read_text()
    )
    assert johnson38["qsopt_status"] == "OPTIMAL"
    assert johnson38["integer_shell_upper"] == 971390652389758748
    assert all(johnson38["independent_exact_primal_checks"].values())
    assert all(johnson38["independent_exact_dual_checks"].values())

    residual38 = json.loads(
        model_path.with_name("johnson_residual_n254_w36_d38.json").read_text()
    )
    assert residual38["qsopt_status"] == "OPTIMAL"
    assert residual38["integer_shell_upper"] == 22789092036927452
    assert all(residual38["independent_exact_primal_checks"].values())
    assert all(residual38["independent_exact_dual_checks"].values())

    generated = model_path.parent
    h38_certificate = verify_scaled_lp(
        model_path,
        generated / "max_h_38_scaled_rational.lp",
        generated / "max_h_38_scaled_rational.sol",
    )
    assert h38_certificate["derived_A38_C_cap"] == 773397229452928
    expected_shell_caps = {
        40: 4210950726984874,
        42: 38281539335776556,
        44: 1144487702470434098,
        46: 7169795968863348791,
        48: 56162197953802356662,
        50: 645473488233155997799,
    }
    for weight, expected in expected_shell_caps.items():
        certificate = verify_scaled_lp(
            model_path,
            generated / f"max_c_{weight}_scaled_rational.lp",
            generated / f"max_c_{weight}_scaled_rational.sol",
        )
        assert certificate["physical_lattice_cap"] == expected
    cap_summary = json.loads(
        (generated / "exact_bch_lp_caps.json").read_text()
    )
    assert abs(cap_summary["combined_margin_bits"] - 32.70601117191865) < 1e-12
    assert abs(cap_summary["weight_38_term_margin_bits"] - 32.78173761656781) < 1e-12
    verify_locator_reduction()

    locator_bounds = json.loads(
        (generated / "locator_incidence_bounds.json").read_text()
    )
    assert locator_bounds["exact_lp_locator_cap"] == 3703262943449
    assert (
        locator_bounds["degenerate_family_g_prime_zero"]["locator_count_upper"]
        == 87550900
    )
    assert abs(locator_bounds["remaining_gap_bits"] - 7.658719250509947) < 1e-12
    assert (
        locator_bounds["smooth_24_point_target"][
            "smooth_factor_budget_after_degenerate_cap"
        ]
        > 6.069
    )

    locator_jacobian = json.loads(
        (generated / "locator_incidence_jacobian.json").read_text()
    )
    assert locator_jacobian["source_affine_orbits"] == 138
    assert locator_jacobian["normalized_locator_polynomials"] == 5244
    assert locator_jacobian["agreement_count_histogram"] == {"37": 5244}
    assert locator_jacobian["interpolant_degree_histogram"] == {
        "16": 7,
        "17": 31,
        "18": 5206,
    }
    assert locator_jacobian["jacobian_rank_histogram"] == {"6": 20976}

    locator_monte_carlo = json.loads(
        (generated / "locator_extension_monte_carlo.json").read_text()
    )
    assert locator_monte_carlo["samples"] == 5000000
    assert locator_monte_carlo["agreement_count_upper"] == 37
    sixth_moment = locator_monte_carlo["binomial_moment_comparison"][6]
    assert sixth_moment["order"] == 6
    assert 0.8 < sixth_moment["ratio"] < 1.3
    assert sixth_moment["normal_95_upper_ratio"] < 1.5
    degree_strata = {
        row["degree"]: row
        for row in locator_monte_carlo["interpolant_degree"]["strata"]
    }
    assert sum(row["samples"] for row in degree_strata.values()) == 5000000
    assert degree_strata[17]["global_reference_contribution_ratio"] < 0.01
    assert degree_strata[18]["normal_95_upper_conditional_ratio"] < 1.5
    ridge_strata = {
        row["in_ridge_image"]: row
        for row in locator_monte_carlo[
            "affine_ridge_image_by_linear_coefficient"
        ]["strata"]
    }
    assert sum(row["samples"] for row in ridge_strata.values()) == 5000000
    assert ridge_strata[True]["normal_95_upper_conditional_ratio"] < 1.5
    assert abs(
        sum(row["global_reference_contribution_ratio"] for row in ridge_strata.values())
        - sixth_moment["ratio"]
    ) < 1e-12

    locator_ladder = json.loads(
        (generated / "locator_size_ladder.json").read_text()
    )
    q8, q32, q128, q256 = locator_ladder["families"]
    assert q8["normalized_locator_list_size"] == 1
    assert q32["normalized_locator_list_size"] == 5
    assert q32["factorial_moment_comparison"][4]["ratio_exact"] == (
        "1048576/525915"
    )
    assert q128["normalized_locator_list_size"] == 330
    assert q128["factorial_moment_comparison"][6]["ratio"] < 1.0
    assert q256["factorial_moment_comparison"][6]["ratio"] < 1.3
    subfield_contributions = sum(
        stratum["global_reference_contribution_ratio"]
        for stratum in q256["middle_subfield"]["strata"]
    )
    assert abs(subfield_contributions - sixth_moment["ratio"]) < 1e-12

    monomial_spectrum = json.loads(
        (generated / "monomial_146_spectrum.json").read_text()
    )
    assert monomial_spectrum["differential_uniformity"] == 6
    assert monomial_spectrum["maximum_absolute_walsh_coefficient"] == 64
    assert monomial_spectrum["walsh_coefficient_histogram"] == {
        "-32": 4080,
        "-16": 13260,
        "0": 26775,
        "16": 17340,
        "32": 3570,
        "64": 255,
    }

    threshold = json.loads(
        (generated / "random_inner_threshold_outward.json").read_text()
    )
    assert threshold["classification"].startswith("directed outward certificate")
    assert threshold["a38_sufficient_cap"]["certified_integer_cap"] == (
        "3827351840403"
    )
    r6_threshold = Fraction(
        int(threshold["r6_factor_threshold_lower"]["numerator"]),
        int(threshold["r6_factor_threshold_lower"]["denominator"]),
    )
    assert Fraction(6) < r6_threshold < Fraction(6099, 1000)
    assert threshold["constant_checks"]["6"]["certified_below_integer_cap"]
    assert (
        threshold["smooth_r6_factor_threshold_after_degenerate_cap"]["decimal"]
        > 6.069
    )

    low_degree_slices = json.loads(
        (generated / "low_degree_146_slices.json").read_text()
    )
    assert [
        row["maximum_number_of_roots"] for row in low_degree_slices["slices"]
    ] == [4, 7, 10]
    assert (
        low_degree_slices["exceptional_walsh_ridge"][
            "signed_total_over_all_256_linear_coefficients"
        ]
        == -64
    )

    moment_lp = json.loads((generated / "rs_coset_moment_lp.json").read_text())
    assert all(moment_lp["exact_certificate_checks"].values())
    assert 3.3e7 < moment_lp["incidence_ratio_R6_upper"]["decimal"] < 3.4e7
    assert list(map(int, moment_lp["nonzero_extremizer"])) == [
        *range(18),
        37,
    ]

    symmetric_incidence = json.loads(
        (generated / "symmetric_incidence_identity.json").read_text()
    )
    assert symmetric_incidence["known_locator_subsets_checked"] == 37
    assert symmetric_incidence["random_sets_checked"] == 1000
    assert symmetric_incidence["random_admissible_sets_seen"] == 0

    binary_family = json.loads(
        (generated / "binary_coefficient_146_family.json").read_text()
    )
    assert binary_family["family_size"] == 1 << 18
    assert binary_family["polynomials_with_at_least_24_agreements"] == 82
    assert binary_family["locator_polynomials_with_37_agreements"] == 10
    assert binary_family["R6_contribution"]["decimal"] < 4e-9

    frobenius_incidence = json.loads(
        (generated / "frobenius_invariant_incidence.json").read_text()
    )
    assert frobenius_incidence["invariant_24_sets"] == 5365
    assert frobenius_incidence["admissible_24_sets"] == 112
    assert frobenius_incidence["new_affine_orbits"] == 4
    assert frobenius_incidence["new_affine_orbit_sizes"] == [65280] * 4
    assert frobenius_incidence["new_rigorous_shell_consequences"] == {
        "additional_A38_P_words": 261120,
        "improved_A38_P_lower": 7768320,
        "improved_h38_lower": 30464,
        "improved_A38_C_lower": 944384,
    }
    assert frobenius_incidence["R6_contribution"]["decimal"] < 2e-17

    frobenius4_incidence = json.loads(
        (generated / "frobenius4_invariant_incidence.json").read_text()
    )
    assert frobenius4_incidence["invariant_24_sets"] == 267516561
    assert frobenius4_incidence["admissible_24_sets"] == 67900
    assert frobenius4_incidence["distinct_admissible_polynomials"] == 40712
    assert frobenius4_incidence["locator_polynomials_with_37_agreements"] == 34
    assert frobenius4_incidence["distinct_affine_orbits_represented_by_locators"] == 15
    assert frobenius4_incidence["new_affine_orbits"] == 10
    assert frobenius4_incidence["new_affine_orbit_sizes"] == [65280] * 10
    assert frobenius4_incidence["new_rigorous_shell_consequences"] == {
        "additional_A38_P_words": 652800,
        "improved_A38_P_lower": 8421120,
        "improved_h38_lower": 33024,
        "improved_A38_C_lower": 1023744,
    }
    assert 1.03 < frobenius4_incidence["conditional_GF4_rank_model_ratio"]["decimal"] < 1.05
    assert frobenius4_incidence["GF4_coefficient_family_I24_contribution"] == 121384172168
    assert frobenius4_incidence["GF4_coefficient_family_R6_contribution"]["decimal"] < 1.2e-8
    assert frobenius4_incidence["GF4_endpoint_R6_contribution"][
        "fraction_of_GF4_family_I24"
    ] > 0.997
    assert frobenius4_incidence["invariant_set_R6_contribution"]["decimal"] < 7e-15

    frobenius16_sample = json.loads(
        (generated / "frobenius16_locator_extension_monte_carlo.json").read_text()
    )
    assert frobenius16_sample["samples"] == 10000000
    assert frobenius16_sample["new_affine_orbits"] == 9
    assert frobenius16_sample["new_affine_orbit_sizes"] == [65280] * 9
    assert frobenius16_sample["new_GF16_endpoint_locators_from_affine_closure"] == 64
    assert frobenius16_sample["certified_GF16_endpoint_inventory_after_closure"] == {
        "affine_orbits_containing_GF16_locators": 30,
        "normalized_GF16_locators": 190,
    }
    assert frobenius16_sample["weighted_invariant_extension_estimator"][
        "normal_95_upper_ratio"
    ] < 1.02
    assert frobenius16_sample["new_rigorous_shell_consequences"] == {
        "additional_A38_P_words": 587520,
        "improved_A38_P_lower": 9008640,
        "improved_h38_lower": 35328,
        "improved_A38_C_lower": 1095168,
    }

    print("PASS: EBCH(256,71) exact spectrum / MacWilliams checks")
    print("PASS: d(EBCH(256,71)^perp) = 16")
    print("PASS: EBCH(256,187)^perp exact spectrum / MacWilliams checks")
    print("PASS: EBCH(256,187) coefficients w=38..60")
    print("PASS: Wambach P/Q seed weights and affine stabilizers")
    print("PASS: BCH chain and transitive 255-element quotient shift orbit")
    print("PASS: exact LP export includes seed bounds and shell lattices")
    print("PASS: calibrated affine-orbit bounds are present in the exact LP")
    print("PASS: exact Johnson LP gives A_52(C) <= 2979928035058718446462766169")
    print("PASS: exact Johnson LP gives A_54(C) <= 47955805775385501045909339217")
    print("PASS: exact diagnostic Johnson LP gives A_38(C) <= 971390652389758748")
    print("PASS: exact residual Johnson LP gives A_38(C) <= 128630323189940096")
    print("PASS: seven exact BCH-sandwich LP primal-dual certificates")
    print("  A_38(C) <= 773397229452928")
    for weight, cap in expected_shell_caps.items():
        print(f"  A_{weight}(C) <= {cap}")
    print("PASS: Wambach minimum word satisfies the normalized locator reduction")
    print("PASS: exact locator incidence thresholds and degenerate-family bound")
    print("PASS: diagnostic locator agreement and Jacobian-rank receipt")
    print("PASS: five-million-sample locator-incidence Monte Carlo receipt")
    print("PASS: exact/sampled locator size ladder and GF(16) stratification")
    print("PASS: exact differential and Walsh spectra of x^146")
    print("PASS: directed RandomStepConv threshold certifies R_6 <= 6")
    print("PASS: exact degree-1/2/3 perturbation slices and ridge cancellation")
    print("PASS: exact primal-dual Reed-Solomon coset moment bound")
    print("PASS: symmetric-function form of the 24-point incidence equations")
    print("PASS: exact binary/GF(4) and sampled GF(16) Frobenius strata")
    print("PASS: twenty-three new certified affine weight-38 orbits")
    print()
    print("Rigorous seed consequences:")
    print("  A_38(P) >= 65280")
    print("  H_38 >= 256")
    print("  A_38(C) >= 7936 and A_38(C) is a multiple of 3968")
    print("  A_40(Q) >= 65280, hence A_40(C) >= 65280")
    print()
    print("Rigorous calibrated orbit-search consequences:")
    print("  A_38(P) >= 7507200")
    print("  H_38 >= 29440 and A_38(C) >= 912640")
    print("  Frobenius enumeration improves A_38(P) >= 7768320")
    print("  H_38 >= 30464 and A_38(C) >= 944384")
    print("  GF(4) enumeration improves A_38(P) >= 8421120")
    print("  H_38 >= 33024 and A_38(C) >= 1023744")
    print("  GF(16) witnesses improve A_38(P) >= 9008640")
    print("  H_38 >= 35328 and A_38(C) >= 1095168")
    print("  A_40(Q) >= 12990720")
    print("  H_40 >= 163072 and A_40(C) >= 18045952")
    print("  H_42 >= 1067264 and A_42(C) >= 42289664")
    print("  H_44 >= 5977088 and A_44(C) >= 267281408")
    print("  H_46 >= 27937536 and A_46(C) >= 1247821056")
    print("  H_48 >= 114343168 and A_48(C) >= 5525363968")
    print("  H_50 >= 28884992 and A_50(C) >= 1574477312 (radius four)")


if __name__ == "__main__":
    main()
