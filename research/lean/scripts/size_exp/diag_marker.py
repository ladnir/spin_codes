#!/usr/bin/env python3
"""Size experiment (A): why is the marker (Chernoff) bound ~75 bits/region looser
than the exact coefficient at G22?  Untrusted diagnostics only; nothing emitted.

Computes, for the G22 tilt 26214425/2^30, alpha 3/5, n = NP = 1024 pairs/region:
  * structure of the 11-state paired polynomial P(t) (degree <= 32);
  * per q: Chernoff t*, log2 of exact [P^n]_q entries vs t*^-q [P(t*)^n] entries;
  * the tilted distribution of g under P~ = t*^g P[g]/rho(t*) (row 0): where the
    mass sits (bimodality = non-convex rate = what the marker loses).
"""
import sys
from fractions import Fraction as Q
from math import log2
from pathlib import Path

import numpy as np
from flint import fmpz_poly

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / 'scripts'))
import ladder_tail22 as LT  # noqa: E402

E = 22
L, NP, D = LT.geom(E)
SP = 400


def rho(P, t):
    M = sum((t ** g) * P[g] for g in range(33))
    return max(abs(np.linalg.eigvals(M)))


def chernoff_t(P, q, n):
    # minimize -q log t + n log rho(t) over log-t grid then golden refine
    best = None
    for x in np.linspace(-12, 2, 2801):
        t = 2.0 ** x
        v = -q * x + n * log2(rho(P, t))
        if best is None or v < best[0]:
            best = (v, t)
    return best


def tilted_power(P, t, n, gmax):
    """Fixed-point tilted power (P~)^n truncated at gmax, P~[g] = t^g P[g]/rho(t)."""
    r = rho(P, t)
    Pt = [P[g] * (t ** g) / r for g in range(33)]
    M = [[fmpz_poly([int(np.ceil(Pt[g][i, j] * 2.0 ** 60)) << (SP - 60) for g in range(33)])
          for j in range(11)] for i in range(11)]
    for _ in range(n.bit_length() - 1):
        out = [[None] * 11 for _ in range(11)]
        for i in range(11):
            for j in range(11):
                acc = fmpz_poly([])
                for k in range(11):
                    acc += fmpz_poly((M[i][k] * M[k][j]).coeffs()[:gmax + 1])
                out[i][j] = fmpz_poly([-((-int(x)) >> SP) for x in acc.coeffs()])
        M = out
    A = np.full((gmax + 1, 11, 11), -np.inf)
    for i in range(11):
        for j in range(11):
            c = M[i][j].coeffs()
            for g, x in enumerate(c[:gmax + 1]):
                x = int(x)
                if x > 0:
                    A[g, i, j] = log2(x) - SP if x.bit_length() < 1000 else (x.bit_length() - SP)
    return A


def main():
    lam = Q(26214425, 1 << 30)
    al = 0.6
    f = LT.fam(lam)
    P = LT.paired_poly(f, al)
    nz = (P > 0)
    print('support (g ranges) per entry (i->k):')
    for i in range(11):
        row = []
        for k in range(11):
            gs = np.nonzero(nz[:, i, k])[0]
            row.append(f'{gs.min()}-{gs.max()}' if len(gs) else '.')
        print(f'  {i:2d}: ' + ' '.join(f'{s:>6}' for s in row))
    qs = [int(x) for x in sys.argv[1:]] or [352, 546, 760]
    for q in qs:
        v, t = chernoff_t(P, q, NP)
        print(f'\nq={q}: Chernoff per-region log2 t^-q rho^n = {v:.1f} at t = 2^{log2(t):.3f}')
        gmax = min(32 * NP, 3 * q + 64)
        A = tilted_power(P, t, NP, gmax)
        # in tilted measure Chernoff value is 0 bits; A[q] entries give the per-entry loss
        print('  log2 tilted [P~^n]_q row 0:', ' '.join(f'{A[q, 0, k]:7.1f}' for k in range(11)))
        print('  log2 tilted [P~^n]_q max entry:', f'{A[q].max():.1f}')
        Tq = np.exp2(A[q])
        print(f'  per-region loss of marker vs exact: -log2 rho(T_q) = {-log2(max(abs(np.linalg.eigvals(Tq)))):.2f} bits')
        # where is row-0 mass (sum over k) as function of g
        row = np.logaddexp2.reduce(A[:, 0, :], axis=1)
        top = np.argsort(-row)[:5]
        print('  row0 mass peaks at g =', [(int(g), round(float(row[g]), 1)) for g in top])
        for g in sorted({0, q // 4, q // 2, q, (3 * q) // 2, 2 * q, gmax}):
            if g <= gmax:
                print(f'   g={g:5d}: log2 row0 mass {row[g]:9.1f}')
        tot = np.logaddexp2.reduce(row)
        print(f'  log2 total row0 mass (<= gmax) {tot:.2f}')
        # mean g under tilted measure
        w = 2.0 ** (row - row.max())
        print(f'  tilted mean g = {np.sum(np.arange(gmax + 1) * w) / w.sum():.1f}')


if __name__ == '__main__':
    main()
