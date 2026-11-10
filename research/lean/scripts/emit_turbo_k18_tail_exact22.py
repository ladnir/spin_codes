"""Untrusted producer: exact-coefficient K18 tail certificates (alpha = 1).

Per family: base (all 17 selected rows, scale 2^S), 7 squarings of A(x)
truncated at degree D = q_max (R = 128), and per q: U0 = ceil(last_q / C(2048,q)),
six scaled region squarings, final numerator a_q over 2^KT. Uses the same
z as the family (z_num / 2^80) and a scaled z^d trace. Lean checks everything.
"""
import json
import sys
from fractions import Fraction as Q
from math import comb, log2
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / 'SpinCodes/Turbo/Generated/K18Tail'
S = 704                 # s = 16^16 * 2^384 * 2^256
DEN_BITS = 384
MULT_BITS = S - 64 - DEN_BITS
KT = 256
D_CUT = 52428
FAMILIES = {1: (3, 5), 2: (6, 10), 3: (11, 19), 4: (20, 34), 5: (35, 51), 6: (52, 61),
            8: (144, 157), 9: (158, 185), 10: (186, 224), 11: (225, 279)}


def cdiv(a, b):
    return -((-a) // b)


def phys_int(fam):
    den = 1 << DEN_BITS
    births = {0: [den] + [0] * 10, 1: fam['birth_one_caps'], 2: fam['birth_two_caps']}
    g = lambda l, i: l[i] if i < len(l) else 0
    def f(j, i, k):
        if i == 0:
            if j < 3:
                return g(births[j], k)
            if k == 0:
                return g(fam['zero_hold_caps'], j)
            if k == 1:
                return g(fam['first_caps'][j], 0)
            return 0
        if k == 0:
            return g(fam['return_caps'][j], i)
        if k == 2:
            return g(fam['first_caps'][j], i)
        return 0
    return f


def row_sel(pi, d, i, k):
    return sum(comb(d, x) * 15 ** x * 16 ** (16 - d) * pi(x, i, k) for x in range(17))


def mm(A, B):
    return [[sum(A[i][k] * B[k][j] for k in range(11)) for j in range(11)] for i in range(11)]


def madd(A, B):
    return [[A[i][j] + B[i][j] for j in range(11)] for i in range(11)]


ZERO = [[0] * 11 for _ in range(11)]


def poly_sq(P, D, s):
    from flint import fmpz_mat
    F = [fmpz_mat(M) for M in P]
    out = []
    for d in range(D + 1):
        acc = fmpz_mat(11, 11)
        for r in range(d + 1):
            acc += F[r] * F[d - r]
        out.append([[cdiv(int(acc[i, j]), s) for j in range(11)] for i in range(11)])
    return out


def ztrace(znum):
    m, tr = 1, []
    for b in bin(D_CUT)[2:]:
        P = m * m * (znum if b == '1' else 1)
        sh = max(0, P.bit_length() - 512)
        m = P >> sh
        tr.append((b == '1', sh, m))
    E = 0
    for bit, sh, _ in tr:
        E = 2 * E + (80 if bit else 0) - sh
    return tr, m, E


def region_chain(U, s, E0=0):
    stages, E = [], E0
    for _ in range(6):
        P = mm(U, U)
        mx = max(max(r) for r in P)
        cur = max(1, mx.bit_length()) - S
        shift = cur - (S + 8)
        up, down = (shift, 0) if shift >= 0 else (0, -shift)
        U = [[cdiv(x * (1 << down), s * (1 << up)) for x in r] for r in P]
        stages.append((down, up, U))
        E = 2 * E + up - down
    return stages, E


def lean_mat(M):
    return '[' + ','.join('[' + ','.join(map(str, r)) + ']' for r in M) + ']'


def lean_poly(P):
    return '[' + ',\n '.join(lean_mat(M) for M in P) + ']'


def main(only=None):
    s = 1 << S
    diag = json.loads((ROOT / '.lake/turbo/k18_coarse_full_diagnostic.json').read_text())
    sys.path.insert(0, str(ROOT.parent / 'workstreams/k18_margin60'))
    man = {}
    for k, (qlo, qhi) in FAMILIES.items():
        if only and k not in only:
            continue
        fam = json.loads((ROOT / f'.lake/turbo/profile_family_k18_tail{k}.json').read_text())
        pi = phys_int(fam)
        D = qhi
        binom = [comb(16, d) for d in range(17)]
        base = [[[binom[d] * row_sel(pi, d, i, kk) * (1 << MULT_BITS) for kk in range(11)]
                 for i in range(11)] if d <= 16 else ZERO for d in range(D + 1)]
        cache = ROOT / f'.lake/turbo/tail_stages_{k}.json'
        if cache.exists():
            stages = json.loads(cache.read_text())
        else:
            stages, P = [], base
            for _ in range(7):
                P = poly_sq(P, D, s)
                stages.append(P)
            cache.write_text(json.dumps(stages))
        last = stages[-1]
        tr, m, Ez = ztrace(fam['z_num'])
        rec = {}
        for q in range(qlo, qhi + 1):
            Cq = comb(2048, q)
            mx = max(max(r) for r in last[q])
            shift = (S + 8) - (max(1, mx.bit_length()) - Cq.bit_length())
            down0, up0 = (shift, 0) if shift >= 0 else (0, -shift)
            U0 = [[cdiv(x * (1 << down0), Cq * (1 << up0)) for x in r] for r in last[q]]
            rst, E = region_chain(U0, s, up0 - down0)
            sumU = sum(rst[-1][2][0])
            T = 256 * q + Ez + E + KT - S
            num = Cq * sumU * (Q(2) ** T)
            a = cdiv(num.numerator, num.denominator * 65535 ** (8 * q) * m) if False else \
                -((-(num.numerator)) // (num.denominator * 65535 ** (8 * q) * m))
            rec[q] = dict(U0=U0, stages=rst, E=E, a=a, down0=down0, up0=up0)
        man[k] = dict(qlo=qlo, qhi=qhi, D=D, m=m, Ez=Ez,
                      bits={q: -log2(r['a']) + KT if r['a'] > 0 else None for q, r in rec.items()})
        print(f'family {k} tilt {fam["selection_only_lambda"]}: q {qlo}..{qhi}', flush=True)
        for q in (qlo, (qlo + qhi) // 2, qhi):
            from replay_fractional import fraction, bits as rbits
            print(f'  q {q}: ours {man[k]["bits"][q]:.3f} bits, replay {rbits(fraction(diag["occupancy_uppers"][str(q)])):.3f}',
                  flush=True)
        emit(k, fam, base, stages, tr, rec, D, qlo, qhi)
    (OUT / 'exact_manifest.json').write_text(json.dumps(man, indent=1) + '\n')


def emit(k, fam, base, stages, tr, rec, D, qlo, qhi):
    ns = f'Spin.Turbo.Hi.ProfileFamilyK18Tail{k}22'
    pre = f'SpinCodes.Turbo.Generated.K18Tail.Tail{k}'
    w = lambda name, text: (OUT / f'Tail{k}{name}.lean').write_text(text, encoding='utf-8', newline='\n')
    w('ExBase', f'import SpinCodes.Turbo.SparseListDefs22\n\nnamespace {ns}.Ex\n'
      f'def base : List (List (List Int)) := {lean_poly(base)}\nend {ns}.Ex\n')
    for t, P in enumerate(stages):
        w(f'ExSq{t}', f'import SpinCodes.Turbo.SparseListDefs22\n\nnamespace {ns}.Ex\n'
          f'def sq{t} : List (List (List Int)) := {lean_poly(P)}\nend {ns}.Ex\n')
    imports = '\n'.join([f'import {pre}ExBase'] + [f'import {pre}ExSq{t}' for t in range(7)] +
                        [f'import {pre}BudgetCheck', 'import SpinCodes.Turbo.Hi.RowSelected22',
                         'import SpinCodes.Turbo.SparseListPolyD22'])
    trace = '[' + ', '.join(f'(false, sq{t}, [])' for t in range(7)) + ']'
    w('ExChain', f'''{imports}

namespace {ns}.Ex
open Spin.Turbo.Hi

def births (n : Nat) : List Int :=
  match n with
  | 0 => birthZeroCaps22
  | 1 => birthOneCaps
  | _ => birthTwoCaps

def physInt (j i k : Nat) : Int := profileRowIntAll22 births firstCaps returnCaps zeroHoldCaps j i k

def row (d i k : Nat) : Int := rowSel physInt d i k

def trace : List (Bool × List (List (List Int)) × List (List (List Int))) := {trace}

set_option maxHeartbeats 0 in
theorem baseChecked : sparseRegionalBaseCheckDL22 {D} row (2 ^ {MULT_BITS}) base = true := by
  decide +kernel

set_option maxHeartbeats 0 in
theorem chainChecked : sparseRegionalPowerCheckDL22 {D} (2 ^ {S}) base base trace = true := by
  decide +kernel

#print axioms baseChecked
#print axioms chainChecked
end {ns}.Ex
''')
    ztr = '[' + ', '.join(f'({str(b).lower()}, {sh}, {mm_})' for b, sh, mm_ in tr) + ']'
    qdefs = []
    for q in range(qlo, qhi + 1):
        r = rec[q]
        st = '[' + ', '.join(f'({dn}, {up}, {lean_mat(U)})' for dn, up, U in r['stages']) + ']'
        qdefs.append(f'def U0_{q} : List (List Int) := {lean_mat(r["U0"])}\n'
                     f'def dn0_{q} : Nat := {r["down0"]}\ndef up0_{q} : Nat := {r["up0"]}\n'
                     f'def st_{q} : List (Nat × Nat × List (List Int)) := {st}\n'
                     f'def a_{q} : Int := {r["a"]}\n')
    w('ExQData', f'import SpinCodes.Turbo.SparseListDefs22\n\nnamespace {ns}.Ex\n'
      f'def ztrace : List (Bool × Nat × Int) := {ztr}\n\n' + '\n'.join(qdefs) + f'end {ns}.Ex\n')


if __name__ == '__main__':
    main(set(map(int, sys.argv[1:])) or None)
