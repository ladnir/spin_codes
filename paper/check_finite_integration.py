"""Check manuscript numbers and selected-map transcription against pinned inputs.

This is an integration check, not a new numerical replay of the distance proof.
Run from any directory with Python 3.10 or newer.
"""
from collections import Counter
from fractions import Fraction
from pathlib import Path
import hashlib
import json
import math
import re
import sys

ROOT = Path(__file__).resolve().parents[1]
PAPER = ROOT / "paper"
BRIDGE = ROOT / "workstreams/bch_rm2sub_bridge"
sys.set_int_max_str_digits(0)


def read(path):
    return json.loads(path.read_text(encoding="utf-8"))


def require(condition, message):
    if not condition:
        raise AssertionError(message)


def find_union(value):
    for key in ("union_upper", "covered_union_upper", "full_union_upper"):
        if key in value:
            row = value[key]
            return Fraction(int(row["numerator"]), int(row["denominator"]))
    if "coverage" in value:
        return find_union(value["coverage"])
    raise AssertionError("Missing exact union")


finite = (PAPER / "finite_certificates.tex").read_text()
appendix = (PAPER / "finite_appendix.tex").read_text()
implementation = (PAPER / "implementation.tex").read_text()
rows = re.findall(r"^(16|18|20|22|24) & (.*?) & (.*?) & ([0-9.]+) \\\\", finite, re.M)
require(len(rows) == 5, "Expected five margin table rows")
margin_checks = []
for m_text, l_text, h_text, display in rows:
    m = int(m_text)
    filename = (f"t128_s19_m{m}_full_split_coverage_v1.json" if m < 20
                else f"t128_s19_m{m}_ladder_full_v1.json")
    ledger = read(BRIDGE / "generated" / filename)
    spec = ledger["instance"]
    require(int(l_text.replace(r"\,", "")) == spec["rows"] == 2**(m-7), "Wrong row count")
    require(int(h_text.replace(r"\,", "")) == spec["cutoff"] == 2**(m+1)//10, "Wrong cutoff")
    bound = find_union(ledger)
    require(bound < Fraction(1, 2**40), "Insufficient exact margin")
    margin = math.log2(bound.denominator) - math.log2(bound.numerator)
    require(abs(float(display) - margin) < 0.00000051, "Misrounded table margin")
    match = re.search(r"\(" + str(m) + r",([0-9.]+)\)", finite)
    require(match is not None and abs(float(match[1])-margin) < 1e-9, "Wrong figure point")
    margin_checks.append({"m": m, "margin": margin})
    if m >= 20:
        for key in ("q1_upper", "sparse_rest_upper", "dense_upper"):
            value = ledger[key]
            component = Fraction(int(value["numerator"]), int(value["denominator"]))
            bits = math.log2(component.denominator)-math.log2(component.numerator)
            require(f"{bits:.6f}" in appendix, "Wrong appendix component margin")
            if key == "q1_upper":
                require(0 <= bits-margin < 1e-6, "Q1 dominance claim does not hold")

selection_path = BRIDGE / "generated/larger_state_inputs_v1/t128_s19_selection.json"
selection = read(selection_path)
require(hashlib.sha256(selection_path.read_bytes()).hexdigest()
        == "7b58cef9156db676bbecbe19179057abc774ee765cd6251c7a86a5ac34eab236",
        "Wrong selected map")
generators = [int(v, 16) for v in selection["selected"]["A_generator_words_hex"]]
transcribed = re.findall(r"^(\d+) & \\texttt\{([0-9a-f]+)\}", appendix, re.M)
require([int(i) for i, _ in transcribed] == list(range(19)), "Wrong generator indices")
require([int(v, 16) for _, v in transcribed] == generators, "Wrong generator word")
require(all((a & b).bit_count() % 2 == 0 for a in generators for b in generators), "CA != 0")
pivots = {}
for word in generators:
    while word:
        pivot = word.bit_length()-1
        if pivot not in pivots:
            pivots[pivot] = word
            break
        word ^= pivots[pivot]
require(len(pivots) == 19, "Wrong image rank")
spectrum = Counter({0: 1})
word = 0
for i in range(1, 2**19):
    word ^= generators[(i & -i).bit_length()-1]
    spectrum[word.bit_count()] += 1
expected = {r["weight"]: r["count"] for r in read(
    BRIDGE / "generated/larger_state_inputs_v1/t128_s19_a_spectrum.json")["spectrum"]}
require(dict(spectrum) == expected, "Wrong image spectrum")
require(spectrum == {0: 1, 48: 5040, 56: 110848, 64: 292510,
                     72: 110848, 80: 5040, 128: 1}, "Wrong printed spectrum")
kernel = read(BRIDGE / "generated/larger_state_inputs_v1/t128_s19_b_kernel_spectrum.json")
for row in kernel["by_total_weight"]:
    j = row["total_weight"]
    total = sum(count * sum((-1)**v * math.comb(w, v) * math.comb(128-w, j-v)
            for v in range(max(0, j-128+w), min(j, w)+1))
            for w, count in spectrum.items())
    require(total % 2**19 == 0 and total//2**19 == row["kernel_words"], "Wrong kernel spectrum")
require(min(r["total_weight"] for r in kernel["by_total_weight"]
            if r["total_weight"] and r["kernel_words"]) == 6, "Wrong kernel distance")

manifest = read(ROOT / "workstreams/bare_bch_rm2sub/generated/MANIFEST.json")
for key in ("p_generator", "q_generator"):
    require(manifest["bch"][key][2:] in appendix, "Wrong BCH polynomial")
require(manifest["t128_s19"]["modulus_low"] == "0x27", "Wrong state field modulus")

perf = read(ROOT / "workstreams/bare_bch_rm2sub/PERFORMANCE.json")
timing_checks = []
for name, values in (("t128_s19", ["0.560", "2.322", "11.259"]),
                     ("t64_s20", ["0.625", "2.577", "12.059"])):
    for m, display in zip((16, 18, 20), values):
        row = next(r for r in perf["summary"] if r["configuration"] == name and r["m"] == m)
        require(f'{row["median_ms"]:.3f}' == display and display in implementation, "Wrong timing")
        timing_checks.append({"configuration": name, "m": m, "ms": display})
for name in ("structured_spin", "structured_appendix", "finite_certificates", "implementation"):
    prose = (PAPER / f"{name}.tex").read_text()
    require("under selection" not in prose and r"\papertodo{" not in prose, "Stale finite placeholder")
    require("Toeplitz" not in prose, "Discarded design in manuscript")
selected = next(r for r in perf["summary"]
                if r["configuration"] == "t128_s19" and r["m"] == 20)
for key, digits in (("setup_median_ms", 2), ("run_min_median_ms", 3),
                    ("run_max_median_ms", 3)):
    require(f'{selected[key]:.{digits}f}' in implementation, "Wrong implementation cost or variation")
require(f'{selected["retained_setup_bytes"]/2**20:.2f}' in implementation,
        "Wrong retained setup memory")
require(selected["workspace_bytes"] == 40*2**20, "Wrong workspace")
main = (PAPER / "main.tex").read_text()
require(main.index(r"\input{structured_spin}") < main.index(r"\input{finite_certificates}")
        < main.index(r"\input{scaling_complexity}"), "Wrong section order")
print(json.dumps({"status": "MANUSCRIPT_INTEGRATION_CHECK_PASSED",
    "full_interval_replay": False, "margin_checks": margin_checks,
    "timing_checks": timing_checks, "generator_words_checked": 19,
    "image_states_enumerated": 2**19, "kernel_spectrum_checked": True}, indent=2))
