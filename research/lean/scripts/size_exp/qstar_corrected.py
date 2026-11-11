#!/usr/bin/env python3
"""Size experiment (A): marker crossover q* with the Lean-exact normalizer.

TailMarkerFracAB22 / pairedRegionDP22_marker_pow64_le bound the tail term
    C(L,q) X^a sum_k (DP_q/C(L,q))^64_{0k} <= X^a C(L,q)^-63 t^-64q (M_t^{64n})_{0.}
(final check `tailMarkerFracFinalCheck22`: ... <= a C^63 tn0^(64q)), but the planners
(ladder_tail22.marker_term, ladder_marker_plan22.best) use C(L,q)^-62, i.e. they are
pessimistic by log2 C(L,q) bits.  This script recomputes q* both ways over the
families already cached in .lake/turbo/q1screen/hi (read-only; nothing written there).

Usage: python scripts/size_exp/qstar_corrected.py e [e ...]
"""
import json
import sys
from fractions import Fraction as Q
from math import comb, lgamma, log, log2
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / 'scripts'))
import emit_turbo_k18_tail_exact22 as X  # noqa: E402

lch = lambda n, k: (lgamma(n + 1) - lgamma(k + 1) - lgamma(n - k + 1)) / log(2)
C1 = 256 - 8 * log2(65535)
FDIR = ROOT / '.lake/turbo/q1screen/hi'
BC = Path(__file__).resolve().parent / '_bcache'


def fam_B(path):
    """17 float 11x11 B matrices (alpha independent), cached as .npy next to this script."""
    BC.mkdir(exist_ok=True)
    cp = BC / (path.stem + '.npz')
    if cp.exists():
        z = np.load(cp)
        return z['B'], float(z['lz'])
    f = json.loads(path.read_text())
    for k, v in list(f.items()):
        if isinstance(v, str) and v[:1] in '[{':
            pass
    pi = X.phys_int(f)
    den = float(16 ** 16) * 2.0 ** 384
    B = np.array([[[X.row_sel(pi, a, i, k) for k in range(11)] for i in range(11)] for a in range(17)],
                 dtype=float) / den
    lz = log2(int(f['z_num']) / 2 ** 80)
    np.savez(cp, B=B, lz=lz)
    return B, lz


def paired(B, al):
    P = np.zeros((33, 11, 11))
    for a in range(17):
        for b in range(17):
            P[a + b] += comb(16, a) * comb(16, b) * np.power(B[a] @ B[b], al)
    return P


def log2rho(P, ts):
    lt = np.log2(ts)
    top = np.where(ts > 1, 32, 0)
    g = np.arange(33)
    W = np.exp2((g[None, :] - top[:, None]) * lt[:, None])        # (T,33)
    M = np.einsum('tg,gij->tij', W, P)
    r = np.abs(np.linalg.eigvals(M)).max(axis=1)
    return top * lt + np.log2(np.maximum(r, 1e-300))


def load_sets():
    files = sorted(FDIR.glob('*.json'))
    ALS = [Q(k, 20) for k in range(8, 21)]
    tq = [Q(2.0 ** (x / 16)).limit_denominator(4096) for x in range(-192, 16 * 56)]
    ts = np.array([float(x) for x in tq])
    sets = []
    for p in files:
        try:
            B, lz = fam_B(p)
        except Exception as ex:  # noqa
            print('skip', p.name, ex)
            continue
        for al in ALS:
            sets.append((p.stem, al, lz, log2rho(paired(B, float(al)), ts)))
    return sets, tq, np.log2(ts)


def main():
    es = [int(x) for x in sys.argv[1:]] or [22]
    sets, tq, lgt = load_sets()
    print(f'{len(sets)} (family, alpha) sets', flush=True)
    for e in es:
        cfg = json.loads((ROOT / f'.lake/turbo/ladder/g{e}/config.json').read_text())
        L, D = cfg['L'], cfg['D']
        NP = L // 32
        target = -(cfg['q1_bits'] + 20)

        def best(q0, q1, cC):
            out = (1e18, None)
            a0, a1 = cC * lch(L, q0), cC * lch(L, q1)
            for name, al, lz, lr in sets:
                fa = float(al)
                v = np.maximum(a0 + fa * (q0 * C1 - D * lz) - 64 * q0 * lgt + 64 * NP * lr,
                               a1 + fa * (q1 * C1 - D * lz) - 64 * q1 * lgt + 64 * NP * lr)
                i = int(np.argmin(v))
                if v[i] < out[0]:
                    out = (float(v[i]), (name, al, tq[i]))
            return out

        qs = sorted({int(round(2 ** (x / 8))) for x in range(8, 8 * (L.bit_length()))} | {L})
        qs = [q for q in qs if 3 <= q <= L]
        print(f'G{e}: L={L} target {target:.2f}', flush=True)
        for cC, tag in ((-62, 'planner C^-62'), (-63, 'Lean C^-63')):
            vals = {q: best(q, q, cC)[0] for q in qs}
            bad = [q for q in qs if vals[q] > target]
            qstar = (max(bad) + 1) if bad else 3
            # refine qstar downward past the last failing grid point
            if bad:
                lo = max(bad)
                hi = min(q for q in qs if q > lo)
                while hi - lo > 1:
                    m = (lo + hi) // 2
                    if best(m, m, cC)[0] <= target:
                        hi = m
                    else:
                        lo = m
                qstar = hi
            print(f'  [{tag}] pointwise q* = {qstar}; sample:', ' '.join(
                f'{q}:{vals[q]:.0f}' for q in qs if q in (3, 8, 16, 32, 64, 128, 256, 362, 512, 1024, 2048)), flush=True)
            if bad:
                print(f'    failing grid q: {bad[:12]}{" ..." if len(bad) > 12 else ""}', flush=True)
        # greedy marker ranges from the corrected q* (Lean normalizer); count + worst
        q0, rng = qstar, []
        while q0 <= L:
            lo, hi, ch = q0, L, None
            while lo <= hi:
                mid = (lo + hi) // 2
                v, arg = best(q0, mid, -63)
                if v <= target:
                    ch, lo = (mid, v, arg), mid + 1
                else:
                    hi = mid - 1
            if ch is None:
                print(f'    range gap at q0={q0}: {best(q0, q0, -63)}', flush=True)
                break
            rng.append((q0, ch[0], round(ch[1], 1), ch[2][0], str(ch[2][1]), str(ch[2][2])))
            q0 = ch[0] + 1
        print(f'  [Lean C^-63] {len(rng)} marker ranges from q*={qstar}:', flush=True)
        for r in rng:
            print('    ', r, flush=True)


if __name__ == '__main__':
    main()
