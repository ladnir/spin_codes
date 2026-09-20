"""Exact sample budgets for a direct BCH weight-38 endpoint test.

This analyzes a prospective IID experiment. The retained fixed-seed MT19937
run is diagnostic data, not a certified source of independent random bases.
All decisions use rational arithmetic; floats only label the JSON output.
"""

from __future__ import annotations

import hashlib
import json
import math
from fractions import Fraction
from math import comb
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
GENERATED = ROOT / "generated"


def ceil_fraction(x: Fraction) -> int:
    return -(-x.numerator // x.denominator)


def log_unit_interval(x: Fraction, terms: int = 64) -> tuple[Fraction, Fraction]:
    """Enclose log(x) for 1 <= x <= 2 by the positive atanh series."""
    assert 1 <= x <= 2 and terms > 0
    z = (x - 1) / (x + 1)
    lower = 2 * sum((z ** (2*j + 1) / (2*j + 1) for j in range(terms)), Fraction())
    remainder = 2 * z ** (2*terms + 1) / ((2*terms + 1) * (1-z*z))
    return lower, lower + remainder


def log_integer(n: int) -> tuple[Fraction, Fraction]:
    assert n >= 1
    k = n.bit_length() - 1
    a, b = log_unit_interval(Fraction(2))
    c, d = log_unit_interval(Fraction(n, 1 << k))
    return k*a + c, k*b + d


def neg_log_one_minus(p: Fraction, terms: int = 4) -> tuple[Fraction, Fraction]:
    assert 0 < p < 1 and terms > 0
    lower = sum((p**j / j for j in range(1, terms+1)), Fraction())
    return lower, lower + p**(terms+1) / ((terms+1) * (1-p))


def ratio_record(x: Fraction) -> dict:
    return {"numerator": str(x.numerator), "denominator": str(x.denominator),
            "approximate": float(x)}


def check_small_field_identity(ladder: dict) -> list[dict]:
    checks = []
    for row in ladder["families"]:
        if row["field_size"] not in (8, 32):
            continue
        q, t, d = row["field_size"], row["base_points"], row["locator_degree"]
        histogram = {int(k): v for k, v in row["extra_agreement_histogram"].items()}
        bases, locators = comb(q-1, t), row["normalized_locator_list_size"]
        assert sum(histogram.values()) == bases
        assert histogram.get(d-t, 0) == locators * comb(d, t)
        checks.append({"q": q, "endpoint_bases": histogram[d-t],
                       "locators": locators, "passed": True})
    assert len(checks) == 2
    return checks


def main() -> None:
    threshold_path = GENERATED / "random_inner_threshold_outward.json"
    sample_path = GENERATED / "locator_extension_monte_carlo.json"
    ladder_path = GENERATED / "locator_size_ladder.json"
    threshold = json.loads(threshold_path.read_text())
    sample = json.loads(sample_path.read_text())
    ladder = json.loads(ladder_path.read_text())
    assert tuple(sample[k] for k in (
        "field_size", "field_modulus", "target_exponent",
        "base_points_per_interpolant", "agreement_count_upper")) == (256, 333, 146, 18, 37)
    hist = {int(k): v for k, v in sample["extra_agreement_histogram"].items()}
    n = sample["samples"]
    assert all(0 <= k <= 19 and isinstance(v, int) and v >= 0 for k, v in hist.items())
    assert sum(hist.values()) == n and n > 0
    hits = hist.get(19, 0)
    assert hits == 0, "This analysis implements the zero-endpoint acceptance rule only."

    cap = int(threshold["a38_sufficient_cap"]["certified_integer_cap"])
    # A38(C)=3968*L37/19 and gcd(3968,19)=1 force L37 in 19*Z.
    assert math.gcd(3968, 19) == 1
    first_bad_a38 = 3968 * (cap // 3968 + 1)
    first_bad_list = 19 * (cap // 3968 + 1)
    endpoint_mass_per_locator = Fraction(comb(37, 18), comb(255, 18))
    p_bad = first_bad_list * endpoint_mass_per_locator
    decay_lo, decay_hi = neg_log_one_minus(p_bad)
    plans = []
    for label, inverse_alpha in (("95 percent", 20), ("99 percent", 100),
                                 ("2^-40 error", 1 << 40), ("2^-80 error", 1 << 80)):
        log_lo, log_hi = log_integer(inverse_alpha)
        # Certify both sufficiency of n_min and failure of n_min-1.
        n_min = ceil_fraction(log_hi / decay_lo)
        assert n_min * decay_lo >= log_hi
        assert (n_min - 1) * decay_hi < log_lo
        assert n * decay_hi < log_lo
        plans.append({"label": label, "alpha": f"1/{inverse_alpha}",
                      "minimum_fresh_samples_if_zero_hits": n_min,
                      "minimality_checked_exactly": True,
                      "existing_run_insufficient_even_under_IID": True})

    log20_lo, log20_hi = log_integer(20)
    # Since (1-p)^n <= exp(-np), p <= log(20)/n is a conservative
    # 95% zero-hit upper confidence rule under the IID model.
    p95_upper = log20_hi / n
    list95_upper = p95_upper / endpoint_mass_per_locator
    a38_95_upper = 3968 * (list95_upper // 19)
    assert a38_95_upper > cap

    checks = check_small_field_identity(ladder)
    # Small exact Bernoulli cases independently check the series direction
    # and integer cutoff logic against direct rational powers.
    for p, denom in ((Fraction(1, 2), 20), (Fraction(1, 4), 100)):
        a, b = neg_log_one_minus(p, terms=128)
        c, d = log_integer(denom)
        m = ceil_fraction(d/a)
        assert (m-1)*b < c
        assert (1-p)**m <= Fraction(1, denom) < (1-p)**(m-1)

    result = {
        "classification": "exact counting identity and prospective IID test budgets; no executed certificate",
        "scope": "fixed BCH [256,128,38], current RandomStepConv occupation-one shell cap",
        "source_sha256": {p.name: hashlib.sha256(p.read_bytes()).hexdigest()
                          for p in (threshold_path, sample_path, ladder_path, Path(__file__))},
        "identity": "Pr[N_B=37] = L37 * C(37,18) / C(255,18)",
        "a38_cap": cap,
        "first_bad_a38_on_lattice": first_bad_a38,
        "first_bad_L37_on_lattice": first_bad_list,
        "endpoint_mass_per_locator": ratio_record(endpoint_mass_per_locator),
        "first_bad_endpoint_probability": ratio_record(p_bad),
        "negative_log_no_hit_probability_per_sample": {
            "lower": ratio_record(decay_lo), "upper": ratio_record(decay_hi)},
        "acceptance_rule": "At a predeclared sample count, accept the cap only if no endpoint was observed; otherwise inconclusive.",
        "plans": plans,
        "existing_diagnostic": {
            "samples": n, "endpoint_hits": hits, "seed": sample["seed"],
            "randomness": "fixed-seed std::mt19937_64; IID guarantee is not established",
            "zero_hit_probability_at_first_bad_shell_approximate": math.exp(-n*float(decay_lo)),
            "hypothetical_IID_conservative_95pct_a38_upper": int(a38_95_upper),
            "certifies_current_cap": False},
        "checks": {"small_field_double_counting": checks,
                   "small_exact_bernoulli_cutoffs": True,
                   "all_sample_count_minimality_inequalities": True},
        "execution_requirements": [
            "Freeze the code, field representation, cap, sample count, and acceptance rule before fresh sampling.",
            "Model and audit the random-bit source; a fixed MT19937 seed does not establish the IID premise.",
            "Use exact uniform subset sampling and exact finite-field evaluation.",
            "Count every endpoint including already known locators; do not filter endpoint hits.",
            "Do not combine invariant-base GF16 samples with uniform-base samples.",
            "Do not retry until zero hits without accounting for the total false-accept probability.",
            "Keep statistical false acceptance separate from the SPIN setup failure bound.",
            "Run no concurrent benchmark."
        ]
    }
    output = GENERATED / "endpoint_sampling_analysis.json"
    output.write_text(json.dumps(result, indent=2) + "\n", encoding="utf-8")
    print(json.dumps({"first_bad_endpoint_probability": float(p_bad), "plans": plans,
                      "existing_diagnostic": result["existing_diagnostic"],
                      "checks": result["checks"], "receipt": str(output)}, indent=2))


if __name__ == "__main__":
    main()
