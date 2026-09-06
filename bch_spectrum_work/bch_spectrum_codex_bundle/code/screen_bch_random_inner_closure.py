"""Sequential binary64 Q1 screen using deterministic BCH caps only.

No parameter change is silently substituted for the existing M22 theorem.
The saved witnesses can seed a separate interval certificate.
"""
from __future__ import annotations
import json
import math
import sys
from pathlib import Path
import numpy as np
from scipy.special import logsumexp

sys.dont_write_bytecode = True
ROOT = Path(__file__).resolve().parents[1]
GEN = ROOT / 'generated'
sys.path.insert(0, str(GEN / 'bch256_q2_envelope_diagnostic_sources'))
import evaluate_ebch128_randomstepconv_g1 as base
from evaluate_ebch128_randomstepconv_q1_exact import uniform_coefficients
from evaluate_bch_q2_envelopes import envelopes
from run_higher_endpoint_preflight import sha, write_new


def region(zero, one, length):
    # Truncated matrix-polynomial powering; coefficient one counts positions.
    result = np.full((2, 2, 2), -np.inf)
    result[0] = base.log_identity()
    power = np.stack((base.log_entries(zero), base.log_entries(one)))
    def mul(a, b):
        return np.stack((base.log_matmul(a[0], b[0]),
                         np.logaddexp(base.log_matmul(a[0], b[1]),
                                      base.log_matmul(a[1], b[0]))))
    exponent = length
    while exponent:
        if exponent & 1:
            result = mul(result, power)
        exponent >>= 1
        if exponent:
            power = mul(power, power)
    result[1] -= math.log(length)
    return result


def main():
    output = GEN / 'bch256_random_inner_closure_screen.json'
    assert not output.exists()
    caps = envelopes()['LP_only']
    weights = sorted(w for w, c in caps.items() if c)
    cases = []
    for memory in (22, 23, 24):
        best = {cutoff: np.full(257, np.inf) for cutoff in (209716, 230687)}
        witness = {cutoff: np.zeros(257) for cutoff in best}
        for i in range(36):
            u = -9.0 + i / 10
            s = math.exp(u)
            zero, one = base.step_matrices(math.exp(-s), memory)
            rz, ra = region(zero, one, 8192)
            values = uniform_coefficients(rz, ra, 256, 256)
            moments = np.logaddexp(values[:, 0, 0], values[:, 0, 1])
            for cutoff in best:
                candidate = math.log(8192) + np.minimum(0., moments + cutoff * s)
                improve = candidate < best[cutoff]
                witness[cutoff][improve] = s
                best[cutoff] = np.minimum(candidate, best[cutoff])
        for cutoff in best:
            terms = [math.log(caps[w]) + best[cutoff][w] for w in weights]
            row = dict(memory_bits=memory, distance_cutoff=cutoff,
                       Q1_margin_bits_diagnostic=-float(logsumexp(terms))/math.log(2),
                       coefficient_rows=[dict(weight=w, cap=caps[w],
                           s=float(witness[cutoff][w]),
                           log2_coefficient=float(best[cutoff][w])/math.log(2)) for w in weights])
            cases.append(row)
            print(json.dumps({k: row[k] for k in ('memory_bits', 'distance_cutoff', 'Q1_margin_bits_diagnostic')}), flush=True)
    sources = [Path(__file__), ROOT/'code/evaluate_bch_q2_envelopes.py',
        GEN/'exact_bch_lp_caps.json', GEN/'johnson_n256_w52_d38.json', GEN/'johnson_n256_w54_d38.json',
        GEN/'bch256_q2_envelope_diagnostic_sources/evaluate_ebch128_randomstepconv_g1.py',
        GEN/'bch256_q2_envelope_diagnostic_sources/evaluate_ebch128_randomstepconv_q1_exact.py']
    write_new(output, dict(classification='Binary64 diagnostic, not a proof',
        statistical_caps_used=False, cases=cases,
        source_sha256={str(p.relative_to(ROOT)): sha(p) for p in sources}))


if __name__ == '__main__':
    main()
