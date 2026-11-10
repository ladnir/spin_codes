"""Untrusted producer/emitter: K18 tail family 7 (tilt 17/200, α = 4/5, q = 62..143).

Certificates for `TailTermFrac22.tail_term_frac_le`: row matrices R (scale den),
fifth-root bounds Y_ab ≥ (R_a R_b)^(4/5) at scale s, paired base polynomial,
6 squarings (64 pairs) truncated at D = 143, per q domination, region chain,
65535^(8q) trace, prefactor split (G, yn), final numerator. Lean checks all.
"""
import json
import sys
from fractions import Fraction as Q
from math import comb, log2
from pathlib import Path

from flint import fmpz, fmpz_mat

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / 'scripts'))
import emit_turbo_k18_tail_exact22 as X  # noqa: E402

OUT = ROOT / 'SpinCodes/Turbo/Generated/K18Tail'
K, QLO, QHI = 7, 62, 143
S, KT, DEN_BITS = 704, 256, 384
D = QHI


def cdiv(a, b):
    return -((-a) // b)


def iroot5_ceil(n):
    if n <= 0:
        return 0
    r = int(fmpz(n).root(5))
    while r ** 5 < n:
        r += 1
    return r


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


def main():
    s = 1 << S
    den = 16 ** 16 * (1 << DEN_BITS)
    fam = json.loads((ROOT / f'.lake/turbo/profile_family_k18_tail{K}.json').read_text())
    pi = X.phys_int(fam)
    R = [[[X.row_sel(pi, a, i, k) for k in range(11)] for i in range(11)] for a in range(17)]
    Rm = [fmpz_mat(M) for M in R]
    Y = [[None] * 17 for _ in range(17)]
    for a in range(17):
        for b in range(17):
            P = Rm[a] * Rm[b]
            Y[a][b] = [[iroot5_ceil(cdiv(int(P[i, k]) ** 4 * s ** 5, den ** 8)) for k in range(11)]
                       for i in range(11)]
    base = []
    for d in range(D + 1):
        acc = [[0] * 11 for _ in range(11)]
        for a in range(17):
            bb = d - a
            if 0 <= bb <= 16:
                c = comb(16, a) * comb(16, bb)
                acc = [[acc[i][k] + c * Y[a][bb][i][k] for k in range(11)] for i in range(11)]
        base.append(acc)
    cache = ROOT / f'.lake/turbo/tail_stages_{K}.json'
    if cache.exists():
        stages = json.loads(cache.read_text())
    else:
        stages, P = [], base
        for _ in range(6):
            P = X.poly_sq(P, D, s)
            stages.append(P)
        cache.write_text(json.dumps(stages))
    last = stages[-1]
    tr, m, Ez = X.ztrace(fam['z_num'])
    diag = json.loads((ROOT / '.lake/turbo/k18_coarse_full_diagnostic.json').read_text())
    sys.path.insert(0, str(ROOT.parent / 'workstreams/k18_margin60'))
    from replay_fractional import fraction, bits as rbits
    recs = {}
    for q in range(QLO, QHI + 1):
        Cq = comb(2048, q)
        U0 = [[cdiv(x, Cq) for x in r] for r in last[q]]
        rst, E = X.region_chain(U0, s)
        sumU = sum(rst[-1][2][0])
        vtr, mv, Ev = scaled_trace(65535, 0, 8 * q)
        T = 256 * q + Ez + Ev
        G = (4 * T) // 5 - 1200
        need = cdiv(1 << (4 * T - 5 * G), (m * mv) ** 4)
        yn = iroot5_ceil(need)
        ex = G + E + KT - S
        num = Cq * yn * sumU * (Q(2) ** ex)
        a = cdiv(num.numerator, num.denominator)
        recs[q] = dict(U0=U0, st=rst, vtr=vtr, G=G, yn=yn, a=a)
        if q in (QLO, (QLO + QHI) // 2, QHI) or q % 20 == 0:
            ours = KT - log2(a) if a > 0 else 999
            print(f'q {q}: ours {ours:.3f} bits, replay {rbits(fraction(diag["occupancy_uppers"][str(q)])):.3f}',
                  flush=True)
    emit(fam, R, Y, base, stages, tr, recs)


def emit(fam, R, Y, base, stages, tr, recs):
    ns = f'Spin.Turbo.Hi.ProfileFamilyK18Tail{K}22'
    pre = f'SpinCodes.Turbo.Generated.K18Tail.Tail{K}'
    w = lambda name, text: (OUT / f'Tail{K}{name}.lean').write_text(text, encoding='utf-8', newline='\n')
    lm, lp = X.lean_mat, X.lean_poly
    w('FrR', f'import SpinCodes.Turbo.SparseListDefs22\n\nnamespace {ns}.Fr\n'
      f'def R : List (List (List Int)) := {lp(R)}\nend {ns}.Fr\n')
    for a in range(17):
        w(f'FrY{a}', f'import SpinCodes.Turbo.SparseListDefs22\n\nnamespace {ns}.Fr\n'
          f'def Y{a} : List (List (List Int)) := {lp(Y[a])}\nend {ns}.Fr\n')
    w('FrBase', f'import SpinCodes.Turbo.SparseListDefs22\n\nnamespace {ns}.Fr\n'
      f'def base : List (List (List Int)) := {lp(base)}\nend {ns}.Fr\n')
    for t, P in enumerate(stages):
        w(f'FrSq{t}', f'import SpinCodes.Turbo.SparseListDefs22\n\nnamespace {ns}.Fr\n'
          f'def sq{t} : List (List (List Int)) := {lp(P)}\nend {ns}.Fr\n')
    imps = '\n'.join([f'import {pre}FrR', f'import {pre}FrBase'] + [f'import {pre}FrY{a}' for a in range(17)] +
                     [f'import {pre}FrSq{t}' for t in range(6)] + ['import SpinCodes.Turbo.TailTermFrac22'])
    ys = '[' + ', '.join(f'Y{a}' for a in range(17)) + ']'
    trace = '[' + ', '.join(f'(false, sq{t}, [])' for t in range(6)) + ']'
    w('FrChain', f'''{imps}

namespace {ns}.Fr
open Spin.Turbo

def Ys : List (List (List (List Int))) := {ys}
def trace : List (Bool × List (List (List Int)) × List (List (List Int))) := {trace}

set_option maxHeartbeats 0 in
theorem rpowChecked : pairedRpowCheckL22 (16 ^ 16 * 2 ^ {DEN_BITS}) (2 ^ {S}) R Ys = true := by decide +kernel
set_option maxHeartbeats 0 in
theorem baseChecked : pairedBaseCheckDL22 {D} Ys base = true := by decide +kernel
set_option maxHeartbeats 0 in
theorem chainChecked : sparseRegionalPowerCheckDL22 {D} (2 ^ {S}) base base trace = true := by decide +kernel

end {ns}.Fr
''')
    ztr = '[' + ', '.join(f'({str(b).lower()}, {sh}, {mm})' for b, sh, mm in tr) + ']'
    qd = []
    for q, r in recs.items():
        st = '[' + ', '.join(f'({dn}, {up}, {lm(U)})' for dn, up, U in r['st']) + ']'
        vt = '[' + ', '.join(f'({str(b).lower()}, {sh}, {mm})' for b, sh, mm in r['vtr']) + ']'
        qd.append(f'def U0_{q} : List (List Int) := {lm(r["U0"])}\ndef st_{q} : List (Nat × Nat × List (List Int)) := {st}\n'
                  f'def vt_{q} : List (Bool × Nat × Int) := {vt}\ndef G_{q} : Int := {r["G"]}\n'
                  f'def yn_{q} : Int := {r["yn"]}\ndef a_{q} : Int := {r["a"]}\n')
    w('FrQData', f'import SpinCodes.Turbo.SparseListDefs22\n\nnamespace {ns}.Fr\n'
      f'def ztrace : List (Bool × Nat × Int) := {ztr}\n\n' + '\n'.join(qd) + f'end {ns}.Fr\n')
    report = {q: r['a'] for q, r in recs.items()}
    (OUT / 'frac_manifest.json').write_text(json.dumps(dict(q=report), indent=1) + '\n')


if __name__ == '__main__':
    main()
