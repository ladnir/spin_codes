#!/usr/bin/env python3
"""Geometry-generic tail planner/producer for the s22 K-ladder (untrusted; Lean checks all).

    python scripts/ladder_tail22.py E plan|families|pairs|exact|marker|emit ...

plan: exact pieces on geometric q-ranges below the marker crossover q*, each
with (tilt, alpha) chosen by an exact-integer screen at its endpoints; above
q*, marker *range* pieces (two endpoint certificates each, justified by
TailMarkerRange22.term_le_endpoints) chosen greedily by a float screen.
The other stages reuse scripts/emit_turbo_k20_tail22.py with this geometry.
"""
from __future__ import annotations

import json
import sys
import types
from decimal import Decimal, ROUND_HALF_EVEN, localcontext
from fractions import Fraction as Q
from math import comb, lgamma, log, log2
from pathlib import Path

import numpy as np
from flint import fmpz_poly

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / 'scripts'))
sys.set_int_max_str_digits(0)
import emit_turbo_k18_tail_exact22 as X  # noqa: E402
from emit_turbo_profile_family22 import produce_family  # noqa: E402

LN2 = log(2)
lch = lambda n, k: (lgamma(n + 1) - lgamma(k + 1) - lgamma(n - k + 1)) / LN2
NL = chr(10)


def geom(e):
    L = (1 << e) // 128
    return L, L // 32, (L * 256) // 10


def load_K(e):
    L, NP, D = geom(e)
    src = (ROOT / 'scripts/emit_turbo_k20_tail22.py').read_text(encoding='utf-8')
    subs = [("OUT = ROOT / 'SpinCodes/Turbo/Generated/K20Tail'", f"OUT = ROOT / 'SpinCodes/Turbo/Generated/G{e}Tail'"),
            ("CACHE = ROOT / '.lake/turbo/k20tail'", f"CACHE = ROOT / '.lake/turbo/ladder/g{e}/tail'"),
            ("GROUPS, NPAIRS, DCUT = 8192, 256, 209715", f"GROUPS, NPAIRS, DCUT = {L}, {NP}, {D}"),
            ("NSQ_EXACT = 8", f"NSQ_EXACT = {NP.bit_length() - 1}"),
            ("NSQ_MARK = 14", f"NSQ_MARK = {(64 * NP).bit_length() - 1}"),
            ("        for q in range(p['qlo'], p['qhi'] + 1):\n            Cq = comb(GROUPS, q)\n            vtr, G, yn = prefactor",
             "        for q in ((p['qlo'], p['qhi']) if p.get('range') else range(p['qlo'], p['qhi'] + 1)):\n            Cq = comb(GROUPS, q)\n            vtr, G, yn = prefactor"),
            ('K20Tail', f'G{e}Tail')]
    for a, b in subs:
        assert a in src, a[:80]
        src = src.replace(a, b)
    m = types.ModuleType(f'ladder_tail_{e}')
    m.__file__ = str(ROOT / 'scripts/emit_turbo_k20_tail22.py')
    exec(compile(src, m.__file__, 'exec'), m.__dict__)
    plan = ROOT / f'.lake/turbo/ladder/g{e}/tail_plan.json'
    m.pieces = lambda: json.loads(plan.read_text())
    sys.modules['emit_turbo_k20_tail22'] = m     # the Lean emitter imports it by this name
    m.CACHE.mkdir(parents=True, exist_ok=True)
    return m


# ---------------------------------------------------------------- screens
CAND = None
FCACHE = ROOT / '.lake/turbo/q1screen/hi'


def fam(lam):
    global CAND
    FCACHE.mkdir(parents=True, exist_ok=True)
    p = FCACHE / f'{lam.numerator}_{lam.denominator}.json'
    if p.exists():
        return json.loads(p.read_text())
    if CAND is None:
        CAND = json.loads((ROOT / '.lake/turbo/coarse_packet_profiles.json').read_text())
    zden = 1 << 80
    with localcontext() as ctx:
        ctx.prec = 100
        znum = int(((-Decimal(lam.numerator) / Decimal(lam.denominator)).exp() * Decimal(zden))
                   .to_integral_value(rounding=ROUND_HALF_EVEN))
    d = produce_family(CAND, znum, zden, 'S', endpoint_bits=384)
    p.write_text(json.dumps(d))
    return d


def paired_poly(f, al):
    pi = X.phys_int(f)
    den = float(16 ** 16) * 2.0 ** 384
    B = [np.array([[X.row_sel(pi, a, i, k) for k in range(11)] for i in range(11)], dtype=float) / den
         for a in range(17)]
    P = np.zeros((33, 11, 11))
    for a in range(17):
        for b in range(17):
            P[a + b] += comb(16, a) * comb(16, b) * np.power(B[a] @ B[b], al)
    return P


SP = 1100


def exact_term(e, q, f, al):
    """log2 estimate of the exact-coefficient term (integer fixed point, upper rounding)."""
    L, NP, D = geom(e)
    P = paired_poly(f, al)
    M = [[fmpz_poly([int(np.ceil(P[g, i, j] * 2.0 ** 60)) << (SP - 60) for g in range(min(q + 1, 33))])
          for j in range(11)] for i in range(11)]
    for _ in range(NP.bit_length() - 1):
        out = [[None] * 11 for _ in range(11)]
        for i in range(11):
            for j in range(11):
                acc = fmpz_poly([])
                for k in range(11):
                    acc += fmpz_poly((M[i][k] * M[k][j]).coeffs()[:q + 1])
                out[i][j] = fmpz_poly([-((-int(x)) >> SP) for x in acc.coeffs()])
        M = out
    I = [[int(M[i][j].coeffs()[q]) if len(M[i][j].coeffs()) > q else 0 for j in range(11)] for i in range(11)]
    mb = max(x.bit_length() for r in I for x in r)
    if mb == 0:
        return -1e18
    sh = max(0, mb - 60)
    C = np.array([[float(x >> sh) for x in r] for r in I])
    la = log2(C.max()) + sh - SP - lch(32 * NP, q)
    A = C / C.max()
    # log2 of sum_k (A^64)_{0k} by renormalized repeated squaring (no underflow)
    lg = 0.0
    for _ in range(6):
        A = A @ A
        m = A.max()
        if m <= 0:
            return -1e18
        A = A / m
        lg = 2 * lg + log2(m)
    s0 = A[0].sum()
    if s0 <= 0:
        return -1e18
    s0 = 2.0 ** 0 * s0
    lgs0 = lg + log2(s0)
    z = f['z_num'] / 2 ** 80
    lgX = q * (256 - 8 * log2(65535)) - D * log2(z)
    return lch(L, q) + al * lgX + 64 * la + lgs0


def marker_rho(f, al, ts):
    P = paired_poly(f, al)
    out = {}
    for t in ts:
        top = 32 if t > 1 else 0          # factor out t^top to avoid overflow
        M = sum((t ** (g - top)) * P[g] for g in range(33))
        r = max(abs(np.linalg.eigvals(M)))
        out[t] = top * log2(t) + (log2(r) if r > 0 else -1e9)
    return out


def marker_term(e, q, f, al, t, lr):
    L, NP, D = geom(e)
    z = f['z_num'] / 2 ** 80
    lgX = q * (256 - 8 * log2(65535)) - D * log2(z)
    # Lean (TailMarkerRange22): term = C * X^al * mass / (C^64 * t^(64q)), i.e. C^-63 overall.
    return lch(L, q) + al * lgX - 64 * lch(L, q) - 64 * q * log2(t) + 64 * NP * lr


def lam_of(c, D):
    return Q(round(Q(c) / D * (1 << 40)), 1 << 40)


def fmt(lam):
    return f'{lam.numerator}/{lam.denominator}'


# ---------------------------------------------------------------- plan
def stage_plan(e, target):
    L, NP, D = geom(e)
    LD = [round(10 * 2 ** (k / 2)) for k in range(0, 38)]   # lambda * d (nats), sqrt(2) steps
    ALS = [Q(1), Q(4, 5), Q(7, 10), Q(3, 5)]
    ts = [2 ** (x / 8) for x in range(-96, 48)] + [2.0 ** k for k in range(7, 53)]
    # marker sets
    msets = []
    for c in LD:
        lam = lam_of(c, D)
        f = fam(lam)
        for al in ALS:
            msets.append((lam, al, f, marker_rho(f, float(al), ts)))
    def mbest(q):
        best = (1e18, None)
        for lam, al, f, rho in msets:
            for t in ts:
                v = marker_term(e, q, f, float(al), t, rho[t])
                if v < best[0]:
                    best = (v, (lam, al, t))
        return best
    # crossover: smallest q with marker <= target for q and beyond (checked on a grid)
    qstar = None
    for q in sorted({int(round(2 ** (x / 4))) for x in range(6, 4 * (L.bit_length()))} | {L}):
        if q < 3 or q > L:
            continue
        if mbest(q)[0] <= target:
            qstar = q
            break
    assert qstar is not None
    print(f'G{e}: marker crossover q* ~ {qstar}', flush=True)
    # refine q* downward to the first failing integer
    lo = max(3, qstar // 2)
    while lo < qstar and mbest(qstar - 1)[0] <= target:
        qstar -= 1
    pieces = []
    # exact pieces below q*
    bounds = [3, 6, 11, 20, 35, 64, 124, 248, 496, 992, 1984, 3968, 7936]
    bounds = [b for b in bounds if b < qstar] + [qstar]
    work = [(lo, hi1 - 1) for lo, hi1 in zip(bounds, bounds[1:])]
    work.reverse()
    while work:
        lo, hi = work.pop()
        best = (1e18, None)
        for c in LD:
            lam = lam_of(c, D)
            f = fam(lam)
            for al in ALS:
                v = max(exact_term(e, lo, f, float(al)), exact_term(e, hi, f, float(al)))
                if v < best[0]:
                    best = (v, (lam, al))
        if best[0] > target and hi > lo:
            mid = int(round((lo * (hi + 1)) ** 0.5))
            mid = min(max(mid, lo + 1), hi)
            print(f'  split exact {lo}..{hi} (worst {best[0]:.1f}) at {mid}', flush=True)
            work.append((mid, hi))
            work.append((lo, mid - 1))
            continue
        lam, al = best[1]
        pieces.append(dict(kind='exact', tilt=fmt(lam), alpha=f'{al.numerator}/{al.denominator}', p=None,
                           qlo=lo, qhi=hi))
        print(f'  exact {lo}..{hi}: tilt {fmt(lam)} alpha {al} worst {best[0]:.1f}', flush=True)
        assert best[0] <= target, (lo, hi, best)
    # marker ranges from q* to L: greedy, endpoints only (log-convex in q)
    q0 = qstar
    while q0 <= L:
        lo_q1, hi_q1 = q0, L
        chosen = None
        while lo_q1 <= hi_q1:
            mid = (lo_q1 + hi_q1) // 2
            ok = None
            for lam, al, f, rho in msets:
                for t in ts:
                    v = max(marker_term(e, q0, f, float(al), t, rho[t]),
                            marker_term(e, mid, f, float(al), t, rho[t]))
                    if v <= target and (ok is None or v < ok[0]):
                        ok = (v, lam, al, t)
            if ok:
                chosen = (mid, ok)
                lo_q1 = mid + 1
            else:
                hi_q1 = mid - 1
        assert chosen is not None, q0
        q1, (v, lam, al, t) = chosen
        tq = Q(t).limit_denominator(256)
        pieces.append(dict(kind='marker', tilt=fmt(lam), alpha=f'{al.numerator}/{al.denominator}', p=None,
                           t=f'{tq.numerator}/{tq.denominator}', qlo=q0, qhi=q1, range=True))
        print(f'  marker {q0}..{q1}: tilt {fmt(lam)} alpha {al} t {tq} worst {v:.1f}', flush=True)
        q0 = q1 + 1
    path = ROOT / f'.lake/turbo/ladder/g{e}/tail_plan.json'
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(pieces, indent=1))


def stage_endpoint(e):
    """Emit Generated/G{e}/Endpoint.lean with the largest provable 5-decimal margin."""
    import re
    from decimal import getcontext
    L, NP, D = geom(e)
    gen = ROOT / f'SpinCodes/Turbo/Generated/G{e}'
    aSum = int(re.search(r'def aSum : Int := (\d+)', (gen / 'Q1Endpoint.lean').read_text(encoding='utf-8')).group(1))
    b2 = int(re.search(r'def b2 : Int := (\d+)', (gen / 'Q2Endpoint.lean').read_text(encoding='utf-8')).group(1))
    pieces = json.loads((ROOT / f'SpinCodes/Turbo/Generated/G{e}Tail/pieces.json').read_text())
    cache = ROOT / f'.lake/turbo/ladder/g{e}/tail'
    tot = 0
    for pc in pieces:
        d = json.loads((cache / f"{pc['kind']}{pc['piece']}.json").read_text())
        recs = {int(q): r['a'] for q, r in d['recs'].items()}
        if pc.get('range'):
            tot += max(recs.values()) * (pc['qhi'] - pc['qlo'] + 1)
        else:
            tot += sum(recs.values())
    T = L * aSum + b2 + tot
    getcontext().prec = 60
    lam = Decimal(256) - Decimal(T).ln() / Decimal(2).ln()
    target = (lam * 100000).to_integral_value(rounding='ROUND_FLOOR') / Decimal(100000)
    mn = None
    for n in range(1, 20001):
        m = int((lam * n).to_integral_value(rounding='ROUND_FLOOR'))
        if Decimal(m) / n >= target:
            mn = (m, n)
            break
    m, n = mn
    assert T ** n * 2 ** m <= 2 ** (256 * n)
    tstr = str(target)
    src = (ROOT / 'scripts/emit_turbo_k20_endpoint22.py').read_text(encoding='utf-8')
    subs = [("'SpinCodes/Turbo/Generated/K20/Endpoint.lean'", f"'SpinCodes/Turbo/Generated/G{e}/Endpoint.lean'"),
            ("'SpinCodes/Turbo/Generated/K20Tail/pieces.json'", f"'SpinCodes/Turbo/Generated/G{e}Tail/pieces.json'"),
            ("NPAIRS, DCUT, GROUPS, R = 256, 209715, 8192, 512", f"NPAIRS, DCUT, GROUPS, R = {NP}, {D}, {L}, {L // 16}"),
            ("cache = ROOT / '.lake/turbo/k20tail'", f"cache = ROOT / '.lake/turbo/ladder/g{e}/tail'"),
            ("""        for q, r in d['recs'].items():
            A[int(q)] = r['a']""", """        rr = {int(q): r['a'] for q, r in d['recs'].items()}
        if p.get('range'):
            for q in range(p['qlo'], p['qhi'] + 1):
                A[q] = max(rr.values())
        else:
            A.update(rr)"""),
            ('.toNat ^ 1622 * 2 ^ 91249 ≤ 2 ^ (256 * 1622)', f'.toNat ^ {n} * 2 ^ {m} ≤ 2 ^ (256 * {n})'),
            ('nat_div_two_pow_le_two_rpow_neg {T} 256 91249 1622', f'nat_div_two_pow_le_two_rpow_neg {{T}} 256 {m} {n}'),
            ('(2 : ℝ) ^ (-(((91249 : ℕ) : ℝ) / ((1622 : ℕ) : ℝ))) ≤ (2 : ℝ) ^ (-(56.25709 : ℝ))',
             f'(2 : ℝ) ^ (-((({m} : ℕ) : ℝ) / (({n} : ℕ) : ℝ))) ≤ (2 : ℝ) ^ (-({tstr} : ℝ))'),
            ('(2 : ℝ) ^ (-(56.25709 : ℝ)) := by', f'(2 : ℝ) ^ (-({tstr} : ℝ)) := by'),
            ('K20Tail', f'G{e}Tail'), ('K20Q1', f'G{e}Q1'), ('K20Q2', f'G{e}Q2'),
            ('Generated.K20.', f'Generated.G{e}.'), ('K20Final', f'G{e}Final'), ('k20_endpoint', f'g{e}_endpoint'),
            ('N = 2097152, cutoff floor(N/10) = 209715', f'N = {1 << (e + 1)}, cutoff floor(N/10) = {D}')]
    for a, b in subs:
        assert a in src, a[:70]
        src = src.replace(a, b)
    exec(compile(src, f'ladder_endpoint_{e}', 'exec'), {'__name__': '__main__', '__file__': str(ROOT / 'scripts/emit_turbo_k20_endpoint22.py')})
    print(f'G{e}: endpoint total bits {float(lam):.6f}; statement 2^-{tstr} via {m}/{n}', flush=True)


if __name__ == '__main__':
    e = int(sys.argv[1])
    cfg = json.loads((ROOT / f'.lake/turbo/ladder/g{e}/config.json').read_text())
    target = -(cfg['q1_bits'] + 20)
    for st in sys.argv[2:]:
        if st == 'plan':
            stage_plan(e, target)
            continue
        if st == 'endpoint':
            stage_endpoint(e)
            continue
        K = load_K(e)
        if st == 'families':
            K.stage_families(K.pieces())
        elif st == 'pairs':
            K.stage_pairs(K.pieces())
        elif st == 'exact':
            K.stage_exact(K.pieces())
        elif st == 'marker':
            K.stage_marker(K.pieces())
        elif st == 'emit':
            import emit_turbo_k20_tail_lean22 as LE  # noqa
            LE.emit_all(K.pieces())
