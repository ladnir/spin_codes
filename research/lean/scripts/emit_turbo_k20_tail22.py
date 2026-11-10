#!/usr/bin/env python3
"""Untrusted producer of the K20 tail certificates (q = 3..8192).

Geometry: R = 512 (n = 256 two-step pairs per region, 64 regions), 8192 groups,
cutoff d = 209715. Pieces follow `.lake/turbo/k20_coarse_full_diagnostic.json`:
* exact pieces (q <= 615): `TailTermFrac22.tail_term_fracAB_le` (alpha = a/b,
  packed squaring chain, 8 squarings to 256 pairs);
* dense pieces (q >= 616): `TailMarkerFracAB22.tail_marker_fracAB_le`
  (paired marker, 14 squarings to 64*256).
Everything emitted is re-checked by Lean; nothing here is trusted.

Stages (cached under .lake/turbo/k20tail/):
  families  -> physical families at the 13 tilts (Hi Fix chain modules)
  pairs     -> R, Y_ab per (family, alpha)
  exact     -> base polynomials, squaring stages, per-q certificates
  marker    -> marker bases, chains, per-q certificates
  emit      -> Lean modules under SpinCodes/Turbo/Generated/K20Tail
"""
from __future__ import annotations

import json
import re
import sys
from decimal import Decimal, ROUND_HALF_EVEN, localcontext
from fractions import Fraction as Q
from math import comb, log2
from pathlib import Path

sys.set_int_max_str_digits(0)
from flint import fmpz, fmpz_mat

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / 'scripts'))
import emit_turbo_k18_tail_exact22 as X  # noqa: E402  (phys_int, row_sel, cdiv, lean_mat/poly)
from emit_turbo_profile_family22 import produce_family, emit_lean  # noqa: E402

OUT = ROOT / 'SpinCodes/Turbo/Generated/K20Tail'
CACHE = ROOT / '.lake/turbo/k20tail'
GROUPS, NPAIRS, DCUT = 8192, 256, 209715
S, KT, DEN_BITS, Z_BITS = 704, 256, 384, 80
NSQ_EXACT = 8          # 2^8 = 256 pairs
NSQ_MARK = 14          # 2^14 = 64 * 256
ENDPOINT_BITS = 384
NL = chr(10)


def cdiv(a, b):
    return -((-a) // b)


def iroot_ceil(n, b):
    if n <= 0:
        return 0
    r = int(fmpz(n).root(b))
    while r ** b < n:
        r += 1
    while r > 0 and (r - 1) ** b >= n:
        r -= 1
    return r


def pieces():
    d = json.loads((ROOT / '.lake/turbo/k20_coarse_full_diagnostic.json').read_text())
    out = []
    for w in d['witnesses']:
        kind = 'exact' if w['kind'] == 'exact_regional' else 'marker'
        out.append(dict(kind=kind, tilt=w['tilt'], alpha=w.get('alpha', '1'),
                        p=w.get('marker_probability'), qlo=int(w['minimum']), qhi=int(w['maximum'])))
    plan = CACHE / 'plan.json'
    if plan.exists():
        # Replace the q=248..615 exact witness (precision floor at degree 615) by
        # exact 248..EXACT_HI plus marker chunks chosen in `stage_plan`.
        extra = json.loads(plan.read_text())
        lo, hi = extra[0]['qlo'], extra[-1]['qhi']
        out = [x for x in out if x['qhi'] < lo] + extra + [x for x in out if x['qlo'] > hi]
    assert out[0]['qlo'] == 3 and out[-1]['qhi'] == GROUPS
    for a, b in zip(out, out[1:]):
        assert a['qhi'] + 1 == b['qlo'], (a, b)
    return out


EXACT_HI = 315


def stage_plan():
    """Choose marker chunks for q = EXACT_HI+1..615 (float screening only;
    the Lean certificates are recomputed exactly and re-checked)."""
    import numpy as np
    from math import lgamma, log
    ps = pieces.__wrapped__() if hasattr(pieces, '__wrapped__') else None
    base = json.loads((ROOT / '.lake/turbo/k20_coarse_full_diagnostic.json').read_text())['witnesses']
    w8 = [w for w in base if int(w['minimum']) == 248][0]
    cands = []
    for tilt, alpha in (('17/200', '7/10'), ('3/20', '3/5')):
        fidx = None
        key = None
        for f in range(1, 14):
            fam = json.loads((CACHE / f'family{f}.json').read_text())
            if fam['selection_only_lambda'] == tilt:
                fidx = f
        key = f'{fidx}_{alpha.replace("/", "o")}'
        pr = json.loads((CACHE / f'pairs_{key}.json').read_text())
        z = json.loads((CACHE / f'family{fidx}.json').read_text())['z_num'] / 2 ** 80
        cands.append((tilt, alpha, pr, z))
    s = 2.0 ** S
    LN2 = log(2)
    lch = lambda n, k: (lgamma(n + 1) - lgamma(k + 1) - lgamma(n - k + 1)) / LN2
    grid = list(range(8, 160))          # t = k / 256
    rho = {}
    for ci, (tilt, alpha, pr, z) in enumerate(cands):
        Ys = [[np.array(pr['Y'][a][b], dtype=float) / s for b in range(17)] for a in range(17)]
        for k in grid:
            t = k / 256
            M = sum(comb(16, a) * comb(16, b) * t ** (a + b) * Ys[a][b] for a in range(17) for b in range(17))
            rho[ci, k] = log2(max(abs(np.linalg.eigvals(M))))
    def bound(ci, k, q):
        tilt, alpha, pr, z = cands[ci]
        al = Q(alpha)
        lgC = lch(GROUPS, q)
        lgX = q * (256 - 8 * log2(65535)) - DCUT * log2(z)
        return lgC + float(al) * lgX - 63 * lgC - 64 * q * log2(k / 256) + 16384 * rho[ci, k]
    out = [dict(kind='exact', tilt=w8['tilt'], alpha=w8['alpha'], p=None, qlo=248, qhi=EXACT_HI)]
    q0 = EXACT_HI + 1
    while q0 <= 615:
        best = None
        for q1 in range(min(615, q0 + 39), q0 - 1, -1):
            cand = min(((max(bound(ci, k, q) for q in range(q0, q1 + 1)), ci, k)
                        for ci in range(len(cands)) for k in grid))
            if cand[0] < -80:
                best = cand
                break
        assert best is not None, q0
        val, ci, k = best
        tilt, alpha, _, _ = cands[ci]
        out.append(dict(kind='marker', tilt=tilt, alpha=alpha, p=None, t=f'{k}/256', qlo=q0, qhi=q1))
        print(f'plan {q0}..{q1}: tilt {tilt} alpha {alpha} t {k}/256 worst {val:.1f}', flush=True)
        q0 = q1 + 1
    (CACHE / 'plan.json').write_text(json.dumps(out, indent=1))


def tilts(ps):
    seen = []
    for p in ps:
        if p['tilt'] not in seen:
            seen.append(p['tilt'])
    return seen


# ---------------------------------------------------------------- families
SMOKE = {'Check': 'PhysicalProfileFamilySmoke22Check', 'BudgetCheck': 'PhysicalSmallBirth22BudgetCheck',
         'BudgetSmoke': 'PhysicalSmallBirth22BudgetSmoke', 'BirthFamily': 'PhysicalProfileBirthFamilySmoke22',
         'Family': 'PhysicalProfileFamilySmoke22'}
CHAIN = ['PhysicalProfileFamilyDefs22', 'PhysicalSmallBirth22BudgetDefs', 'PhysicalSmallBirth22Budget',
         'PhysicalProfileBirthFamily22', 'PhysicalProfileCompleteRows22']


def fam_ns(f):
    return f'Spin.Turbo.Hi.ProfileFamilyK20Tail{f}22'


def fam_mod(f):
    return f'SpinCodes.Turbo.Generated.K20Tail.Tail{f}'


def port_smoke(kind, f):
    t = (ROOT / f'SpinCodes/Turbo/{SMOKE[kind]}.lean').read_text(encoding='utf-8')
    mod = fam_mod(f)
    names = {'PhysicalProfileFamilySmoke22Data': f'{mod}Data',
             'PhysicalProfileFamilySmoke22Check': f'{mod}Check',
             'PhysicalSmallBirth22BudgetCheck': f'{mod}BudgetCheck',
             'PhysicalSmallBirth22BudgetSmoke': f'{mod}BudgetSmoke',
             'PhysicalProfileBirthFamilySmoke22': f'{mod}BirthFamily',
             'PhysicalProfileFamilySmoke22': f'{mod}Family'}
    for old in sorted(names, key=len, reverse=True):
        t = re.sub(rf'^import SpinCodes\.Turbo\.{old}$', f'import {names[old]}', t, flags=re.M)
    for c in CHAIN:
        t = re.sub(rf'^import SpinCodes\.Turbo\.{c}$', f'import SpinCodes.Turbo.Hi.{c}', t, flags=re.M)
    t = re.sub(r'\bSpin\.Numeric\b(?!Hi)', 'Spin.NumericHi', t)
    t = t.replace('namespace Spin.Turbo.ProfileFamilySmoke22', f'namespace {fam_ns(f)}')
    t = t.replace('end Spin.Turbo.ProfileFamilySmoke22', f'end {fam_ns(f)}')
    assert 'ProfileFamilySmoke22' not in t
    (OUT / f'Tail{f}{kind}.lean').write_text(t, encoding='utf-8', newline=NL)


def stage_families(ps):
    OUT.mkdir(parents=True, exist_ok=True)
    (OUT.parent / '.gitignore').write_text('*' + NL + '!.gitignore' + NL, encoding='utf-8')
    CACHE.mkdir(parents=True, exist_ok=True)
    candidate = json.loads((ROOT / '.lake/turbo/coarse_packet_profiles.json').read_text())
    for f, ts in enumerate(tilts(ps), start=1):
        path = CACHE / f'family{f}.json'
        if path.exists():
            continue
        lam = Q(ts)
        zden = 1 << Z_BITS
        with localcontext() as ctx:
            ctx.prec = 100
            znum = int(((-Decimal(lam.numerator) / Decimal(lam.denominator)).exp() * Decimal(zden))
                       .to_integral_value(rounding=ROUND_HALF_EVEN))
        data = produce_family(candidate, znum, zden, fam_ns(f), endpoint_bits=ENDPOINT_BITS)
        data['selection_only_lambda'] = ts
        emit_lean(data, OUT / f'Tail{f}Data.lean', hi=True)
        for kind in SMOKE:
            port_smoke(kind, f)
        path.write_text(json.dumps(data) + NL)
        print(f'family {f}: tilt {ts}', flush=True)


# ---------------------------------------------------------------- pairs
def alpha_ab(a):
    fr = Q(a)
    return fr.numerator, fr.denominator


def pair_key(p, fidx):
    return f'{fidx}_{p["alpha"].replace("/", "o")}'


def stage_pairs(ps):
    fidx = {t: i for i, t in enumerate(tilts(ps), start=1)}
    s = 1 << S
    den = 16 ** 16 * (1 << DEN_BITS)
    done = set()
    for p in ps:
        f = fidx[p['tilt']]
        key = pair_key(p, f)
        path = CACHE / f'pairs_{key}.json'
        if key in done or path.exists():
            done.add(key)
            continue
        ea, eb = alpha_ab(p['alpha'])
        fam = json.loads((CACHE / f'family{f}.json').read_text())
        pi = X.phys_int(fam)
        R = [[[X.row_sel(pi, a, i, k) for k in range(11)] for i in range(11)] for a in range(17)]
        Rm = [fmpz_mat(M) for M in R]
        Y = [[None] * 17 for _ in range(17)]
        for a in range(17):
            for b in range(17):
                P = Rm[a] * Rm[b]
                # (P/den^2)^(ea/eb) <= Y/s  <=  P^ea s^eb <= Y^eb den^(2 ea)
                Y[a][b] = [[iroot_ceil(cdiv(int(P[i, k]) ** ea * s ** eb, den ** (2 * ea)), eb)
                            for k in range(11)] for i in range(11)]
        path.write_text(json.dumps(dict(R=R, Y=Y, ea=ea, eb=eb, f=f)))
        done.add(key)
        print(f'pairs {key}', flush=True)


# ---------------------------------------------------------------- shared numerics
def ztrace(znum):
    m, tr = 1, []
    for b in bin(DCUT)[2:]:
        P = m * m * (znum if b == '1' else 1)
        sh = max(0, P.bit_length() - 512)
        m = P >> sh
        tr.append((b == '1', sh, m))
    E = 0
    for bit, sh, _ in tr:
        E = 2 * E + (Z_BITS if bit else 0) - sh
    return tr, m, E


def scaled_trace(base, bits, d):
    m, tr = 1, []
    for b in bin(d)[2:]:
        P = m * m * (base if b == '1' else 1)
        sh = max(0, P.bit_length() - 512)
        m = P >> sh
        tr.append((b == '1', sh, m))
    E = 0
    for bit, sh, _ in tr:
        E = 2 * E + (bits if bit else 0) - sh
    return tr, m, E


def square_chain(U, steps, E0):
    s = 1 << S
    stages, E = [], E0
    for _ in range(steps):
        P = X.mm(U, U)
        mx = max(max(r) for r in P)
        cur = max(1, mx.bit_length()) - S
        shift = cur - (S + 8)
        up, down = (shift, 0) if shift >= 0 else (0, -shift)
        U = [[cdiv(x * (1 << down), s * (1 << up)) for x in r] for r in P]
        stages.append((down, up, U))
        E = 2 * E + up - down
    return stages, E


def prefactor(q, ea, eb, m, Ez):
    vtr, mv, Ev = scaled_trace(65535, 0, 8 * q)
    T = 256 * q + Ez + Ev
    G = (ea * T) // eb - 1200
    need = cdiv(1 << (ea * T - eb * G), (m * mv) ** ea)
    yn = iroot_ceil(need, eb)
    return vtr, G, yn


# ---------------------------------------------------------------- exact pieces
def poly_sq_trunc(P, D):
    s = 1 << S
    F = [fmpz_mat(M) for M in P]
    out = []
    for d in range(D + 1):
        acc = fmpz_mat(11, 11)
        for r in range(d + 1):
            acc += F[r] * F[d - r]
        out.append([[cdiv(int(acc[i, j]), s) for j in range(11)] for i in range(11)])
    return out


def stage_exact(ps, only=None):
    fidx = {t: i for i, t in enumerate(tilts(ps), start=1)}
    for pidx, p in enumerate(ps, start=1):
        if p['kind'] != 'exact' or (only and pidx not in only):
            continue
        path = CACHE / f'exact{pidx}.json'
        prior = json.loads(path.read_text()) if path.exists() else None
        if prior is not None and prior.get('shifted'):
            continue
        f = fidx[p['tilt']]
        pr = json.loads((CACHE / f'pairs_{pair_key(p, f)}.json').read_text())
        ea, eb, Y = pr['ea'], pr['eb'], pr['Y']
        D = p['qhi']
        base = []
        for d in range(D + 1):
            acc = [[0] * 11 for _ in range(11)]
            for a in range(17):
                bb = d - a
                if 0 <= bb <= 16:
                    c = comb(16, a) * comb(16, bb)
                    acc = [[acc[i][k] + c * Y[a][bb][i][k] for k in range(11)] for i in range(11)]
            base.append(acc)
        scache = CACHE / f'exact{pidx}_stages.json'
        if prior is not None:
            stages = prior['stages']
        elif scache.exists():
            stages = json.loads(scache.read_text())
        else:
            stages, P = [], base
            for t in range(NSQ_EXACT):
                P = poly_sq_trunc(P, D)
                stages.append(P)
                print(f'  piece {pidx}: squaring {t + 1}/{NSQ_EXACT}', flush=True)
            scache.write_text(json.dumps(stages))
        last = stages[-1]
        fam = json.loads((CACHE / f'family{f}.json').read_text())
        tr, m, Ez = ztrace(fam['z_num'])
        recs = {}
        worst = 0.0
        for q in range(p['qlo'], p['qhi'] + 1):
            Cq = comb(GROUPS, q)
            mx = max(max(r) for r in last[q])
            sh = (max(1, mx.bit_length()) - Cq.bit_length()) - (S + 8)   # U0 ~ 2^(S+8)
            up0, dn0 = (sh, 0) if sh >= 0 else (0, -sh)
            U0 = [[cdiv(x * (1 << dn0), Cq * (1 << up0)) for x in r] for r in last[q]]
            rst, E = X.region_chain(U0, 1 << S, up0 - dn0)
            sumU = sum(rst[-1][2][0])
            vtr, G, yn = prefactor(q, ea, eb, m, Ez)
            ex = G + E + KT - S
            num = Cq * yn * sumU * (Q(2) ** ex)
            a = cdiv(num.numerator, num.denominator)
            recs[q] = dict(U0=U0, st=rst, vtr=vtr, G=G, yn=yn, a=a, dn0=dn0, up0=up0)
            if q % 25 == 0 or q == p['qlo']:
                print(f'    q {q}: log2 a-256 = {log2(max(a, 1)) - KT:.2f}, log2 U0max = {log2(max(max(r) for r in U0)) - S:.2f}', flush=True)
            worst = max(worst, log2(max(a, 1)) - KT) if q != p['qlo'] else log2(max(a, 1)) - KT
        path.write_text(json.dumps(dict(base=base, stages=stages, ztr=tr, recs=recs, D=D, ea=ea, eb=eb, f=f,
                                        shifted=True)))
        print(f'piece {pidx} exact q {p["qlo"]}..{p["qhi"]}: worst term 2^{worst:.2f}', flush=True)


# ---------------------------------------------------------------- marker pieces
def marker_t(p):
    if p.get('t'):
        t = Q(p['t'])
        return t.numerator, t.denominator, 0, 0
    pr = Q(p['p'])
    if pr == 1:
        return 1, 1, 24, 0
    t = pr / (1 - pr)
    return t.numerator, t.denominator, 0, 0


def stage_marker(ps, only=None):
    fidx = {t: i for i, t in enumerate(tilts(ps), start=1)}
    s = 1 << S
    for pidx, p in enumerate(ps, start=1):
        if p['kind'] != 'marker' or (only and pidx not in only):
            continue
        path = CACHE / f'marker{pidx}.json'
        if path.exists():
            continue
        f = fidx[p['tilt']]
        pr = json.loads((CACHE / f'pairs_{pair_key(p, f)}.json').read_text())
        ea, eb, Y = pr['ea'], pr['eb'], pr['Y']
        tn0, td0, kp, kn = marker_t(p)
        tn, td = tn0 << kp, td0 << kn
        tot = [[0] * 11 for _ in range(11)]
        for a in range(17):
            for b in range(17):
                c = comb(16, a) * comb(16, b) * tn ** (a + b) * td ** (32 - a - b)
                for i in range(11):
                    for k in range(11):
                        tot[i][k] += c * Y[a][b][i][k]
        # M_t <= tot / (sy td^32) with sy = s; U0 = ceil(2^down s tot / (s td^32 2^up))
        mx = max(max(r) for r in tot)
        cur = max(1, mx.bit_length()) - (32 * td.bit_length())   # ~ log2(s * M_t)
        shift = cur - (S + 8)                                      # U0 ~ 2^(S+8), as in the chain
        up, down = (shift, 0) if shift >= 0 else (0, -shift)
        U0 = [[cdiv((1 << down) * x, td ** 32 * (1 << up)) for x in r] for r in tot]
        rst, E = square_chain(U0, NSQ_MARK, up - down)
        sumU = sum(rst[-1][2][0])
        fam = json.loads((CACHE / f'family{f}.json').read_text())
        tr, m, Ez = ztrace(fam['z_num'])
        recs = {}
        worst = -1e9
        for q in range(p['qlo'], p['qhi'] + 1):
            Cq = comb(GROUPS, q)
            vtr, G, yn = prefactor(q, ea, eb, m, Ez)
            T = G + E + KT + 64 * q * kn - S - 64 * q * kp
            # a >= yn sumU td0^(64q) 2^T / (C^63 tn0^(64q))
            lg = (log2(yn) + log2(sumU) + 64 * q * log2(td0) + T - 63 * log2(Cq) - 64 * q * log2(tn0))
            if lg < -16:
                a = 1
            else:
                num = yn * sumU * td0 ** (64 * q)
                dd = Cq ** 63 * tn0 ** (64 * q)
                a = cdiv(num * (1 << T), dd) if T >= 0 else cdiv(num, dd * (1 << (-T)))
            recs[q] = dict(vtr=vtr, G=G, yn=yn, a=a)
            worst = max(worst, lg)
        path.write_text(json.dumps(dict(U0=U0, up=up, down=down, st=rst, tn0=tn0, td0=td0, kp=kp, kn=kn,
                                        ztr=tr, recs=recs, ea=ea, eb=eb, f=f)))
        print(f'piece {pidx} marker q {p["qlo"]}..{p["qhi"]}: worst log2 term/2^-256 {worst:.1f}', flush=True)


def main():
    ps = pieces()
    stage = sys.argv[1] if len(sys.argv) > 1 else 'all'
    only = set(map(int, sys.argv[2:]))
    if stage == 'plan':
        stage_plan()
        return
    if stage in ('families', 'all'):
        stage_families(ps)
    if stage in ('pairs', 'all'):
        stage_pairs(ps)
    if stage in ('exact', 'all'):
        stage_exact(ps, only)
    if stage in ('marker', 'all'):
        stage_marker(ps, only)
    if stage == 'emit':
        import emit_turbo_k20_tail_lean22 as L
        L.emit_all(ps)


if __name__ == '__main__':
    main()
