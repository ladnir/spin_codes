"""Exact lazy-refresh state law and a binary64 Q1 diagnostic, NOT a certificate.

Each round samples nonzero u, then v uniformly from u-perp, including v=0.
The linear map is I+u v^T. All messages share the sampled maps.
"""
from collections import Counter
from fractions import Fraction
import hashlib
import json
import math
from pathlib import Path
import sys

import numpy as np

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
sys.path.insert(0, str(ROOT/'workstreams/rate_quarter_bch'))
import evaluate_fixed_inner as fixed
import smaller_outer
q1 = fixed.q1


def dot(a, b):
    return (a & b).bit_count() & 1


def update(q, u, v):
    return q ^ (u if dot(q, v) else 0)


def projected_v(u, word):
    """Uniform word -> uniform u-perp; u must be nonzero."""
    if not u:
        raise ValueError('u must be nonzero')
    return word ^ ((u & -u) if dot(u, word) else 0)


def transition_counts(s, q):
    return Counter(update(q, u, v) for u in range(1, 1 << s)
                   for v in range(1 << s) if not dot(u, v))


def cancellation_weights(columns):
    """For input e_p, lazy cancellation occurs only at q=B e_p.

    Return wt(e_p+A B e_p), using the retained A=B^T map.
    """
    assert len(set(columns)) == len(columns) and all(columns)
    return [sum(dot(a, q) ^ (p == j) for j, a in enumerate(columns))
            for p, q in enumerate(columns)]


def epoch_logs(spectrum, columns, lambdas, rounds):
    """Z/D/L weighted-measure envelope; rounds=None means exact refresh.

    L is capped by kappa/M on nonzero states. Lazy surviving mass is sent to
    D, not L. Lazy cancellation uses exact single-input syndrome preimages.
    """
    t = len(columns)
    s = (sum(spectrum.values())+1).bit_length()-1
    den = (1 << s)-1
    assert sum(spectrum.values()) == den
    if rounds is not None and rounds < 1:
        raise ValueError('positive round count required')
    lam = np.asarray(lambdas)
    weights = np.array(list(spectrum))
    counts = np.array(list(spectrum.values()), dtype=float)
    f0 = -lam[:, None]*weights
    f1 = np.logaddexp(np.log(weights/t)-lam[:, None]*(weights-1),
                     np.log(np.maximum(1-weights/t, 1e-300))-lam[:, None]*(weights+1))
    # A weight-t word contributes no second summand.
    if t in spectrum:
        f1[:, list(spectrum).index(t)] = -lam*(t-1)
    m0 = np.logaddexp.reduce(f0+np.log(counts/den), axis=1)
    m1 = np.logaddexp.reduce(f1+np.log(counts/den), axis=1)
    d0, d1 = np.max(f0, axis=1), np.max(f1, axis=1)
    kappa = math.log1p(1/(den-1))
    log_m = math.log(den)
    epsilon = 0. if rounds is None else 2.**-rounds
    lazy = -np.inf if not epsilon else math.log(epsilon)
    fresh = math.log1p(-epsilon)
    cw = np.array(cancellation_weights(columns))
    cd = -lam*min(cw)-math.log(t)
    cl = kappa-log_m+np.logaddexp.reduce(-lam[:, None]*cw, axis=1)-math.log(t)
    zero = np.full((len(lam), 3, 3), -np.inf)
    one = np.full_like(zero, -np.inf)
    zero[:, 0, 0] = 0
    one[:, 0, 1] = -lam
    for matrix, d, m in ((zero, d0, m0), (one, d1, m1)):
        matrix[:, 1, 1] = lazy+d
        matrix[:, 1, 2] = fresh+d
        matrix[:, 2, 1] = lazy+kappa+m
        matrix[:, 2, 2] = fresh+kappa+m
    one[:, 1, 0] = np.logaddexp(lazy+cd, fresh+d1-log_m)
    one[:, 2, 0] = np.logaddexp(lazy+cl, fresh+kappa+m1-log_m)
    return zero, one


def main():
    spectrum, _, sources = fixed.load_inner()
    path = ROOT/'workstreams/bare_bch_rm2sub/generated/manifest.json'
    columns = json.loads(path.read_text())['t128_s19']['columns']
    outer = {w:n for w,n in enumerate(smaller_outer.spectrum()) if w and n}
    tilts = np.arange(-180, 1, dtype=float)/10
    rows, outer_length = 32768, 128
    results = []
    for rounds in (None, 1, 2, 3, 4, 6, 8):
        moments = q1.coefficient_logs(*q1.region_logs(
            *epoch_logs(spectrum, columns, np.exp(tilts), rounds), rows//128), outer_length)
        for delta in (Fraction(33,200), Fraction(19,100)):
            cutoff = outer_length*rows*delta.numerator//delta.denominator
            values = np.minimum(0., moments+cutoff*np.exp(tilts)[:,None])
            witnesses = np.argmin(values, axis=0)
            best = values[witnesses, np.arange(outer_length+1)]
            terms = [math.log(rows*n)+best[w] for w,n in outer.items()]
            dominant = list(outer)[int(np.argmax(terms))]
            result = dict(rounds=rounds, distance_target=str(delta),
                          q1_margin_bits=-float(np.logaddexp.reduce(terms))/math.log(2),
                          dominant_outer_weight=dominant,
                          dominant_log_surprisal=float(tilts[witnesses[dominant]]),
                          witness_at_grid_edge=bool(witnesses[dominant] in (0,len(tilts)-1)))
            results.append(result)
            print(json.dumps(result), flush=True)
    sources += [path, Path(__file__), Path(smaller_outer.__file__), Path(q1.__file__)]
    payload = dict(status='BINARY64_Q1_ONLY_NOT_A_CERTIFICATE', outer=[128,32,32],
                   message_bits=1<<20, t=128, s=19, q1_log_surprisals=tilts.tolist(),
                   cancellation_weight_counts=dict(sorted(Counter(cancellation_weights(columns)).items())),
                   source_sha256={p.relative_to(ROOT).as_posix():hashlib.sha256(p.read_bytes()).hexdigest()
                                  for p in sources}, results=results)
    (HERE/'TRANSVECTION_Q1.json').write_text(json.dumps(payload, indent=2)+'\n', encoding='utf-8')


if __name__ == '__main__':
    main()
