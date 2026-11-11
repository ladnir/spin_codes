#!/usr/bin/env python3
"""Untrusted q=1 margin screen for the s22 Turbo/RS8 code at any K = 2^e.

Uses the exact integer q1 pipeline of emit_turbo_k20_sparse_cert22 (rigorous
upper bounds, same arithmetic Lean re-checks) with geometry R = K / 2048 * 8,
L = 16 R groups, d = floor(256 L / 10), over a geometric tilt grid, choosing
the best tilt per shell. Prints q1 bits (an upper bound on the whole-margin
contribution of q = 1). No Lean claim.
Usage: screen_turbo_q1_geometry22.py E [E ...]   (E = log2 K)
"""
import json
import sys
from decimal import Decimal, ROUND_HALF_EVEN, localcontext
from fractions import Fraction as Q
from math import comb, log2
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / 'scripts'))
import emit_turbo_k20_sparse_cert22 as P  # noqa: E402
from emit_turbo_profile_family22 import produce_family  # noqa: E402

CANDIDATE = json.loads((ROOT / '.lake/turbo/coarse_packet_profiles.json').read_text())
CACHE = ROOT / '.lake/turbo/q1screen'


def family(lam):
    CACHE.mkdir(parents=True, exist_ok=True)
    path = CACHE / f'fam_{lam.numerator}_{lam.denominator}.json'
    if path.exists():
        return json.loads(path.read_text())
    zden = 1 << 80
    with localcontext() as ctx:
        ctx.prec = 100
        znum = int(((-Decimal(lam.numerator) / Decimal(lam.denominator)).exp() * Decimal(zden))
                   .to_integral_value(rounding=ROUND_HALF_EVEN))
    data = produce_family(CANDIDATE, znum, zden, 'Spin.Turbo.Screen')
    path.write_text(json.dumps(data))
    return data


def screen(e):
    K = 1 << e
    L = K // 128          # groups (128 message bits each)
    R = L // 16
    P.R, P.L = R, L
    P.D = (L * 256) // 10
    cnt, den = P.rs8_counts()
    best = [None] * 65
    # tilt grid: lambda * d roughly 30..120 nats for the dominant shells
    grid = sorted({Q(round(x * (1 << 20) / P.D), 1 << 20) for x in
                   [2 ** (k / 4) for k in range(12, 35)] if round(x * (1 << 20) / P.D) > 0})
    for lam in grid:
        fam = family(lam)
        _, _, result = P.regional(fam)
        rows = P.q1_rows(result[0], result[1])
        _, _, zlast = P.power_trace(fam['z_num'])
        for v in range(8, 65):
            t = Q(cnt[v], den) / comb(64, v) / Q(zlast) * sum(rows[64][v])
            if best[v] is None or t < best[v][0]:
                best[v] = (t, lam)
    q1 = L * (sum(Q(cnt[v], den) for v in range(8)) + sum(b[0] for b in best[8:]))
    bits = -log2(q1.numerator) + log2(q1.denominator)
    used = sorted({str(b[1]) for b in best[8:]}, key=lambda s: Q(s))
    print(f'K=2^{e}: L={L} d={P.D} q1 bits {bits:.4f}; tilts used {used}', flush=True)
    return bits


if __name__ == '__main__':
    for e in map(int, sys.argv[1:]):
        screen(e)
