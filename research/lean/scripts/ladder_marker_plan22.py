"""Fine-grid marker-range planner for the s22 ladder tail (untrusted producer).

Writes .lake/turbo/ladder/g{e}/marker_plan.json: greedy marker ranges from q* to L.
Each range uses one (tilt, alpha, t); the marker term is log-convex in q
(TailMarkerRange22), so a range is accepted when both endpoints meet the target.
The float estimate here only proposes ranges; Lean checks every endpoint exactly.

Grids (much finer than ladder_tail22.stage_plan, which misses narrow optima at
L >= 2^17): tilt lambda*d = 10 * 2^(k/8), alpha in 0.05 steps, t in 2^(1/16) steps.
Usage: python scripts/ladder_marker_plan22.py e qstar [kmax]
"""
import json
import sys
from fractions import Fraction as Q
from math import lgamma, log, log2

import numpy as np

sys.path.insert(0, str(__import__('pathlib').Path(__file__).resolve().parent))
import ladder_tail22 as T  # noqa: E402


def lch(n, k):
    return (lgamma(n + 1) - lgamma(k + 1) - lgamma(n - k + 1)) / log(2)


def main():
    e, qstar = int(sys.argv[1]), int(sys.argv[2])
    kmax = int(sys.argv[3]) if len(sys.argv) > 3 else 8 * 22
    L, NP, D = T.geom(e)
    cfg = json.loads((T.ROOT / f'.lake/turbo/ladder/g{e}/config.json').read_text())
    target = -(cfg['q1_bits'] + 20)
    # The float estimate runs ~5-6 bits ahead of the integer producer; keep a margin.
    target -= float(__import__('os').environ.get('TURBO_MARKER_MARGIN', 8))
    ALS = [Q(k, 20) for k in range(8, 21)]          # 0.40 .. 1.00
    # TURBO_MARKER_ALPHA1=1: alpha = 1 only. The prefactor check raises the v-trace
    # mantissa (65535^(8q), unscaled) to the power ea; at q ~ 1.5e5..5e5 any ea > 1
    # exceeds the kernel numeral limit.
    if __import__('os').environ.get('TURBO_MARKER_ALPHA1') == '1':
        ALS = [Q(1)]
    # t values are snapped to small rationals first, so the planned value is the checked one
    tq_all = [Q(2.0 ** (x / 16)).limit_denominator(4096) for x in range(-192, 16 * 56)]
    # TURBO_MARKER_POW2T=1: t restricted to exact powers of two. The checker raises the
    # numerator/denominator of t to the power 64q; at q ~ 5e5 a non-dyadic t needs
    # numerals far beyond the kernel limit, while 2^k folds into the tracked exponent.
    if __import__('os').environ.get('TURBO_MARKER_POW2T') == '1':
        tq_all = [Q(2) ** k for k in range(-12, 56)]
    ts = np.array([float(x) for x in tq_all])
    lgt = np.log2(ts)
    c1 = 256 - 8 * log2(65535)
    sets = []
    for k in range(0, kmax + 1):
        lam = T.lam_of(round(10 * 2 ** (k / 8)), D)
        f = T.fam(lam)
        lz = log2(f['z_num'] / 2 ** 80)
        for al in ALS:
            rho = T.marker_rho(f, float(al), list(ts))
            lr = np.array([rho[t] for t in ts])
            sets.append((lam, al, lz, 64 * NP * lr))
    print(f'G{e}: {len(sets)} marker sets, target {target:.2f}', flush=True)

    def best(q0, q1):
        """Min over sets/t of max(term(q0), term(q1))."""
        out = (1e18, None)
        a0 = lch(L, q0) * -63
        a1 = lch(L, q1) * -63
        for lam, al, lz, lr64 in sets:
            fa = float(al)
            v0 = a0 + fa * (q0 * c1 - D * lz) - 64 * q0 * lgt + lr64
            v1 = a1 + fa * (q1 * c1 - D * lz) - 64 * q1 * lgt + lr64
            v = np.maximum(v0, v1)
            i = int(np.argmin(v))
            if v[i] < out[0]:
                out = (float(v[i]), (lam, al, tq_all[i]))
        return out

    pieces, q0 = [], qstar
    while q0 <= L:
        lo, hi, chosen = q0, L, None
        while lo <= hi:
            mid = (lo + hi) // 2
            v, arg = best(q0, mid)
            if v <= target:
                chosen, lo = (mid, v, arg), mid + 1
            else:
                hi = mid - 1
        assert chosen is not None, ('marker gap at', q0, best(q0, q0))
        q1, v, (lam, al, t) = chosen
        tq = t
        pieces.append(dict(kind='marker', tilt=f'{lam.numerator}/{lam.denominator}',
                           alpha=f'{al.numerator}/{al.denominator}', p=None,
                           t=f'{tq.numerator}/{tq.denominator}', qlo=q0, qhi=q1, range=True))
        print(f'  marker {q0}..{q1}: tilt {float(lam):.6g} alpha {al} t {tq} worst {v:.1f}', flush=True)
        q0 = q1 + 1
    path = T.ROOT / f'.lake/turbo/ladder/g{e}/marker_plan_q{qstar}.json'
    path.write_text(json.dumps(pieces, indent=1))
    print('MARKERDONE', len(pieces), flush=True)


if __name__ == '__main__':
    main()
