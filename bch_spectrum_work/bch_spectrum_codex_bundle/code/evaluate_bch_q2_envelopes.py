"""Q=1,2 diagnostics for fixed BCH256 spectrum upper envelopes.

Snapshots the existing RandomStepConv support-pair evaluator before importing
it. All arithmetic is nearest binary64. A numerical pass is not a certificate.
"""
from __future__ import annotations

import argparse
import hashlib
import itertools
import json
import math
from pathlib import Path
import shutil
import sys
import time

sys.dont_write_bytecode = True
import numpy as np
from scipy.special import logsumexp

ROOT = Path(__file__).resolve().parents[1]
GEN = ROOT / "generated"
UPSTREAM = Path("C:/Users/peter/.codex/worktrees/3061/permute_conv/workstreams/finite_asymptotic_theory")
MODULES = ("evaluate_ebch128_randomstepconv_g1.py",
           "evaluate_ebch128_randomstepconv_q1_exact.py",
           "evaluate_single_random_constituent_q2.py")


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def envelopes():
    low = json.loads((GEN / "exact_bch_lp_caps.json").read_text())
    threshold = json.loads((GEN / "random_inner_threshold_outward.json").read_text())
    upper = {0: 0, 256: 1}
    low_caps = {r["weight"]: r["exact_Aw_C_cap"] for r in low["shells"]}
    for w in range(38, 129, 2):
        cap = min(1 << 128, math.comb(256, w-18)//math.comb(w, w-18))
        if w in low_caps:
            cap = min(cap, low_caps[w])
        if w in (52, 54):
            j = json.loads((GEN / f"johnson_n256_w{w}_d38.json").read_text())
            cap = min(cap, j["integer_shell_upper"])
        upper[w] = upper[256-w] = cap
    endpoint = upper.copy()
    endpoint[38] = endpoint[218] = int(threshold["a38_sufficient_cap"]["certified_integer_cap"])
    proposed = endpoint.copy()
    proposed[40] = proposed[216] = 50_000_000_000_000
    proposed[42] = proposed[214] = 5_000_000_000_000_000
    return {"LP_only": upper, "accepted_endpoint38": endpoint,
            "proposed_40_42": proposed}


def toy_check(base, uniform, pairs):
    zero, one = base.step_matrices(0.91, 4)
    region = uniform(base.log_entries(zero), base.log_entries(one), 5, 2)
    direct = []
    for weight in range(3):
        value = np.zeros((2, 2))
        for subset in itertools.combinations(range(5), weight):
            product = np.eye(2)
            for i in range(5):
                product = product @ (one if i in subset else zero)
            value += product / math.comb(5, weight)
        direct.append(value)
    assert np.max(np.abs(np.exp(region)-np.array(direct))) < 2e-14
    actual = np.exp(pairs(region, 4))
    expected = np.zeros_like(actual)
    for a in range(5):
        for b in range(5):
            for first in itertools.combinations(range(4), a):
                for second in itertools.combinations(range(4), b):
                    product = np.eye(2)
                    for i in range(4):
                        product = product @ direct[int(i in first)+int(i in second)]
                    expected[a,b] += product / (math.comb(4,a)*math.comb(4,b))
    error = float(np.max(np.abs(actual-expected)))
    assert error < 3e-14
    return {"all_25_weight_pairs_checked": True, "maximum_absolute_error": error}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, default=GEN / "bch256_q2_envelope_diagnostic.json")
    parser.add_argument("--tilts", type=float, nargs="+", default=[-9,-8.5,-8,-7.5,-7,-6.5,-6,-5.5,-5])
    args = parser.parse_args()
    assert not args.output.exists(), "Refusing to overwrite a diagnostic run"
    snapshot = args.output.with_suffix("").with_name(args.output.stem+"_sources")
    snapshot.mkdir(parents=True, exist_ok=False)
    for name in MODULES:
        shutil.copyfile(UPSTREAM / name, snapshot / name)
    hashes = {name: digest(snapshot/name) for name in MODULES}
    sys.path.insert(0, str(snapshot))
    import evaluate_ebch128_randomstepconv_g1 as base
    from evaluate_ebch128_randomstepconv_q1_exact import uniform_coefficients as uniform
    from evaluate_single_random_constituent_q2 import pair_support_coefficients as pairs
    checks = toy_check(base, uniform, pairs)
    rows = 8192
    cutoff = ((1 << 21)+9)//10
    best = np.full((257,257), np.inf)
    witness = np.full_like(best, np.nan)
    start = time.monotonic()
    for u in args.tilts:
        s = math.exp(u)
        zero, one = base.step_matrices(math.exp(-s), 22)
        region = uniform(base.log_entries(zero), base.log_entries(one), rows, 2)
        coefficients = pairs(region, 256)
        moments = np.logaddexp(coefficients[...,0,0], coefficients[...,0,1])
        value = np.minimum(0.0, moments+cutoff*s)
        selected = value < best
        witness[selected] = u
        best = np.minimum(best, value)
        print(f"tilt={u} elapsed={time.monotonic()-start:.1f}s", flush=True)
    result = {}
    spectra = envelopes()
    for name, spectrum in spectra.items():
        weights = sorted(w for w,v in spectrum.items() if w and v)
        logs = np.array([math.log(spectrum[w]) for w in weights])
        subtotal = math.log(math.comb(rows,2))+logs[:,None]+logs[None,:]+best[np.ix_(weights,weights)]
        aggregate = float(logsumexp(subtotal))
        top = np.argsort(subtotal.ravel())[-10:][::-1]
        dominant = []
        for i in top:
            a,b = np.unravel_index(i, subtotal.shape)
            dominant.append({"weights": [weights[a],weights[b]],
                             "term_margin_bits": -float(subtotal[a,b])/math.log(2),
                             "tilt": float(witness[weights[a],weights[b]])})
        q1 = float(logsumexp(math.log(rows)+logs+best[weights,0]))
        result[name] = {
            "q2_margin_bits": -aggregate/math.log(2), "q1_margin_bits": -q1/math.log(2),
            "q1_all_one_term_margin_bits": -(math.log(rows)+best[256,0])/math.log(2),
            "dominant_q2_pairs": dominant}
    cache = args.output.with_suffix(".npz")
    np.savez_compressed(cache, best_conditional_log_bound=best, best_tilt=witness)
    assert hashes == {name: digest(snapshot/name) for name in MODULES}
    receipt = {
        "classification": "nearest-binary64 diagnostic; coefficientwise spectrum envelopes, not exact spectra",
        "parameters": {"message_bits": 1 << 20, "output_bits": 1 << 21, "outer_rows": rows,
                       "memory_bits": 22, "cutoff": cutoff, "tilts": args.tilts},
        "probability_space": "Fixed C repeated across rows; independent row and region permutations and RandomStepConv inner setup",
        "counting_law": "For distinct active row positions, ordered codeword weight pairs have multiplicity A_a*A_b; the supports are independent after independent row permutations, even when local messages coincide.",
        "envelope_conditions": {"LP_only": "exact existing individual upper caps",
            "accepted_endpoint38": "also assumes the statistically accepted A38 cap, with its stated RNG qualification",
            "proposed_40_42": "also assumes the still-unproved A40<=5e13 and A42<=5e15"},
        "results": result, "toy_check": checks,
        "limitations": ["No outward rounding", "Only Q=1 and Q=2", "RandomStepConv, not RM2Sub",
                        "Pointwise upper envelopes need not themselves be realizable spectra",
                        "A finite tilt grid gives valid diagnostic witnesses, not guaranteed optima"],
        "source_snapshot": str(snapshot), "source_sha256": hashes,
        "code_sha256": digest(Path(__file__)), "cache_sha256": digest(cache),
        "spectrum_envelopes": spectra, "seconds": time.monotonic()-start}
    args.output.write_text(json.dumps(receipt,indent=2)+"\n",encoding="utf-8")
    print(json.dumps({name: {k:v for k,v in r.items() if k!="dominant_q2_pairs"} for name,r in result.items()},indent=2))


if __name__ == "__main__":
    main()
