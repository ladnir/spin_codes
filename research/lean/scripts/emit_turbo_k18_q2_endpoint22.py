"""Emit the K18 q2 common-marker endpoint (untrusted producer; Lean checks all).

Marker t = 24 on the tilt-3 (lambda = 1/625) regional kernels, a positive
Collatz vector for the 11x11 marker matrix, exact outer weight sum.
"""
import json
import sys
from fractions import Fraction as Q
from math import comb, factorial, log2
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / 'scripts'))
import emit_turbo_k18_sparse_cert22 as P  # noqa: E402

OUT = ROOT / 'SpinCodes/Turbo/Generated/K18'
TEMPLATE = ROOT / 'scripts/templates/k18_q2_endpoint.lean.tmpl'
K = 256
TM = 24
M8 = ((1 << 32) - 1) ** 8
RHO_BITS = 64
VEC_BITS = 60


def main():
    fam = json.loads((ROOT / '.lake/turbo/profile_family_k18_tilt3.json').read_text())
    _, _, result = P.regional(fam)
    S = P.S
    mint = [[result[0][i][j] + 2 * TM * result[1][i][j] + TM * TM * result[2][i][j]
             for j in range(11)] for i in range(11)]
    mf = np.array([[mint[i][j] / S for j in range(11)] for i in range(11)])
    eig = max(np.linalg.eigvals(mf).real)
    best = None
    for gap in (1e-9, 1e-6, 1e-4, 1e-3):
        f = np.linalg.solve((eig * (1 + gap)) * np.eye(11) - mf, np.ones(11))
        if not np.all(f > 0):
            continue
        f = f / f.max()
        vec = [max(1, int(x * (1 << VEC_BITS))) for x in f]
        rhoD = 1 << RHO_BITS
        rhoN = max(-((-rhoD * sum(mint[i][j] * vec[j] for j in range(11))) // (S * vec[i]))
                   for i in range(11))
        bound = Q(rhoN, rhoD) ** 64 * Q(vec[0], min(vec))
        if best is None or bound < best[0]:
            best = (bound, vec, rhoN, rhoD)
    bound, vec, rhoN, rhoD = best
    cnt, den = P.rs8_counts()
    Wd = M8 * factorial(64) * TM ** 64
    nterms = [cnt[v] * factorial(v) * factorial(64 - v) * TM ** (64 - v) for v in range(65)]
    Wn = sum(nterms)
    assert Q(Wn, Wd) == sum(Q(cnt[v], M8) / comb(64, v) / Q(TM) ** v for v in range(65))
    man = json.loads((OUT / 'manifest.json').read_text())
    zlast = man['tilts']['3']['zlast']
    C2 = comb(2048, 2)
    num = C2 * rhoN ** 64 * vec[0] * Wn ** 2 * S * (1 << K)
    dd = rhoD ** 64 * min(vec) * Wd ** 2 * zlast
    b2 = -((-num) // dd)
    q2 = Q(b2, 1 << K)
    print('q2 bits:', -log2(q2.numerator) + log2(q2.denominator))
    lst = lambda xs: '[' + ','.join(map(str, xs)) + ']'
    text = TEMPLATE.read_text(encoding='utf-8')
    subs = {'@MINT@': 'Array.mk ' + lst([x for r in mint for x in r]),
            '@VEC@': 'Array.mk ' + lst(vec), '@RHON@': str(rhoN), '@RHOD@': str(rhoD),
            '@VMIN@': str(min(vec)), '@NTERMS@': lst(nterms), '@WN@': str(Wn),
            '@B2@': str(b2), '@K@': str(K), '@D@': str(P.D), '@R@': str(P.R), '@TM@': str(TM),
            '@T3@': 'Spin.Turbo.ProfileFamilyK18Tilt322'}
    for key, val in subs.items():
        text = text.replace(key, val)
    assert '@' not in text.replace('@[', '')
    (OUT / 'Q2Endpoint.lean').write_text(text, encoding='utf-8', newline='\n')
    man['q2_endpoint'] = dict(K=K, b2=b2, bits=-log2(q2.numerator) + log2(q2.denominator))
    (OUT / 'manifest.json').write_text(json.dumps(man, indent=1, default=str) + '\n')


if __name__ == '__main__':
    main()
