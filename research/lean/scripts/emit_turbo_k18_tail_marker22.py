"""Untrusted producer/emitter: K18 marker-range tail terms (q = 280..2048, α = 1).

Per (family, t) one shared certificate: base U0 with A_t ≤ 2^(up-down)·U0/s,
13 scaled squarings (A_t^8192). Per q: best grid t, final numerator a_q / 2^KT
(checked by `tailMarkerFinalCheck22`). Lean checks everything.
"""
import json
import sys
from fractions import Fraction as Q
from math import comb, log2
from pathlib import Path

from flint import fmpz_mat

sys.set_int_max_str_digits(0)

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / 'scripts'))
import emit_turbo_k18_tail_exact22 as X  # noqa: E402

OUT = ROOT / 'SpinCodes/Turbo/Generated/K18Tail'
S, KT, DEN_BITS = 704, 256, 384
SEGMENTS = [(11, 280, 281), (12, 282, 344), (13, 345, 427), (14, 428, 546), (15, 547, 769),
            (16, 770, 1499), (17, 1500, 2048)]
GRID = 8             # t = 2^(j/8)
LB2 = 256 - 8 * log2(65535)


def cdiv(a, b):
    return -((-a) // b)


def t_parts(j, grid=GRID):
    """t = 2^(j/GRID) ≈ tn0·2^kp / (td0·2^kn), tn0, td0 ≤ 2^12."""
    ip, fr = divmod(j, grid)
    f = Q(2 ** (fr / grid)).limit_denominator(1 << 12)
    tn0, td0 = f.numerator, f.denominator
    kp, kn = (ip, 0) if ip >= 0 else (0, -ip)
    return tn0, td0, kp, kn


def certificate(pi, j, s, grid=GRID):
    den = 16 ** 16 * (1 << DEN_BITS)
    tn0, td0, kp, kn = t_parts(j, grid)
    tn, td = tn0 << kp, td0 << kn
    rows = [[[X.row_sel(pi, a, i, k) for k in range(11)] for i in range(11)] for a in range(17)]
    num = [[sum(comb(16, a) * tn ** a * td ** (16 - a) * rows[a][i][k] for a in range(17))
            for k in range(11)] for i in range(11)]
    dd = den * td ** 16
    mx = max(max(r) for r in num)
    shift = (S + 8) - (max(1, mx.bit_length()) - dd.bit_length())    # value·2^shift ≈ 2^(S+8)
    down, up = (shift, 0) if shift >= 0 else (0, -shift)
    # 2^down·s·num ≤ dd·2^up·U0 and A_t ≈ num/dd, U0 ≈ value·2^(down-up)·s
    U0 = [[cdiv((1 << down) * s * x, dd * (1 << up)) for x in r] for r in num]
    stages, E, U = [], up - down, U0
    for _ in range(13):
        P = fmpz_mat(U) * fmpz_mat(U)
        Pl = [[int(P[i, k]) for k in range(11)] for i in range(11)]
        mxp = max(max(r) for r in Pl)
        cur = max(1, mxp.bit_length()) - S
        sh = cur - (S + 8)
        u2, d2 = (sh, 0) if sh >= 0 else (0, -sh)
        U = [[cdiv(x * (1 << d2), s * (1 << u2)) for x in r] for r in Pl]
        stages.append((d2, u2, U))
        E = 2 * E + u2 - d2
    sumU = sum(U[0])
    return dict(tn0=tn0, td0=td0, kp=kp, kn=kn, down=down, up=up, U0=U0, stages=stages, E=E, sumU=sumU)


def final(q, c, Ez, m):
    C = comb(2048, q)
    T = 256 * q + KT + 64 * q * c['kn'] + Ez + c['E'] - (S + 64 * q * c['kp'])
    num = c['sumU'] * c['td0'] ** (64 * q) * (Q(2) ** T)
    den = 65535 ** (8 * q) * m * C ** 63 * c['tn0'] ** (64 * q)
    return -((-num.numerator) // (num.denominator * den))


def est_bits(q, c, Ez, m):
    """log2-domain estimate of the final bound (for t selection only)."""
    C = comb(2048, q)
    l = (log2(c['sumU']) + c['E'] - S + 64 * q * (log2(c['td0']) + c['kn'] - log2(c['tn0']) - c['kp'])
         + 256 * q - 8 * q * log2(65535) - log2(m) + Ez - 63 * log2(C))
    return -l


def main():
    s = 1 << S
    fams = {k: json.loads((ROOT / f'.lake/turbo/profile_family_k18_tail{k}.json').read_text())
            for k, _, _ in SEGMENTS}
    cover = json.loads((ROOT / '.lake/turbo/marker_cover.json').read_text())
    man = {}
    for k, qlo, qhi in SEGMENTS:
        fam = fams[k]
        pi = X.phys_int(fam)
        tr, m, Ez = X.ztrace(fam['z_num'])
        grid = 64 if k <= 15 else GRID
        global_grid = grid
        cand = {}
        if k == 17:
            certs = {j: certificate(pi, j, s, grid) for j in range(0, grid * 40)}
            for q in range(qlo, qhi + 1):
                cand[q] = sorted(certs, key=lambda jj: -est_bits(q, certs[jj], Ez, m))[:3]
        else:
            for q in range(qlo, qhi + 1):
                j0 = round(grid * log2(cover[str(q)][2]))
                cand[q] = [j0 + dj for dj in range(-3, 4)]
            js = sorted({j for v in cand.values() for j in v})
            certs = {j: certificate(pi, j, s, grid) for j in js}
        recs, choice = {}, {}
        for q in range(qlo, qhi + 1):
            best = None
            for j in cand[q]:
                a = final(q, certs[j], Ez, m)
                if best is None or a < best[0]:
                    best = (a, j)
            recs[q], choice[q] = best
        used = set(choice.values())
        bits = {q: KT - log2(a) if a > 0 else 99999 for q, a in recs.items()}
        worst = min(bits.values())
        print(f'family {k} ({fam["selection_only_lambda"]}) q {qlo}..{qhi}: {len(used)} chains, '
              f'worst {worst:.1f} bits at q {min(bits, key=bits.get)}', flush=True)
        man[k] = dict(qlo=qlo, qhi=qhi, chains=sorted(used), bits=bits)
        emit(k, fam, tr, {j: certs[j] for j in used}, choice, recs, qlo, qhi)
    (OUT / 'marker_manifest.json').write_text(json.dumps(man, indent=1) + '\n')


def lm(M):
    return X.lean_mat(M)


def emit(k, fam, tr, certs, choice, recs, qlo, qhi):
    ns = f'Spin.Turbo.Hi.ProfileFamilyK18Tail{k}22'
    pre = f'SpinCodes.Turbo.Generated.K18Tail.Tail{k}'
    w = lambda name, text: (OUT / f'Tail{k}{name}.lean').write_text(text, encoding='utf-8', newline='\n')
    den = '(16 ^ 16 * endpointDen)'
    for j, c in certs.items():
        tag = f'm{j}' if j >= 0 else f'n{-j}'
        st = '[' + ', '.join(f'({d}, {u}, {lm(U)})' for d, u, U in c['stages']) + ']'
        w(f'Mk{tag}', f'''import {pre}MkRows
import SpinCodes.Turbo.TailMarker22

namespace {ns}.Mk
open Spin.Turbo Spin.Turbo.Hi

def U0_{tag} : List (List Int) := {lm(c['U0'])}
def st_{tag} : List (Nat × Nat × List (List Int)) := {st}

set_option maxHeartbeats 0 in
theorem base_{tag} : tailMarkerBaseCheck22 row ({c['tn0']} * 2 ^ {c['kp']}) ({c['td0']} * 2 ^ {c['kn']})
    {den} (2 ^ {S}) {c['down']} {c['up']} U0_{tag} = true := by decide +kernel

set_option maxHeartbeats 0 in
theorem chain_{tag} : sparseSquareChainCheckL22 (2 ^ {S}) U0_{tag} st_{tag} = true := by decide +kernel

end {ns}.Mk
''')
    ztr = '[' + ', '.join(f'({str(b).lower()}, {sh}, {mm})' for b, sh, mm in tr) + ']'
    w('MkRows', f'''import {pre}Family
import SpinCodes.Turbo.Hi.RowSelected22
import SpinCodes.Turbo.TailMarker22

noncomputable section
namespace {ns}.Mk
open Spin.Turbo Spin.Turbo.Hi

def births (n : Nat) : List Int :=
  match n with
  | 0 => birthZeroCaps22
  | 1 => birthOneCaps
  | _ => birthTwoCaps

def physInt (j i k : Nat) : Int := profileRowIntAll22 births firstCaps returnCaps zeroHoldCaps j i k

def row (d i k : Nat) : Int := rowSel physInt d i k

def ztrace : List (Bool × Nat × Int) := {ztr}

theorem hsmall (j : Fin 3) (i : Fin 11) :
    smallBirthSmokeCaps22 j i = ((births j.val).getD i.val 0 : ℝ) / endpointDen := by
  fin_cases j <;> simp [smallBirthSmokeCaps22, births]

theorem hphys (x : Fin 17) (i k : Fin 11) :
    physicalRows x i k = (physInt x.val i.val k.val : ℝ) / endpointDen :=
  rows_eq smallBirthSmokeCaps22 firstCaps returnCaps zeroHoldCaps endpointDen births hsmall x i k

theorem hrow (d : ℕ) (hd : d ≤ 16) (i k : Fin 11) :
    selectedRows ⟨d, by omega⟩ i k ≤ (row d i.val k.val : ℝ) / ((16 ^ 16 * endpointDen : ℤ) : ℝ) := by
  show turboThinnedRow22 ⟨d, by omega⟩ birthSmokeZ22 physicalRows i k ≤ _
  rw [selected_eq birthSmokeZ22 physicalRows physInt endpointDen (by decide +kernel) hphys d
    (by omega) i k]
  unfold row
  norm_num

theorem hz : birthSmokeZ22 = (zNum : ℝ) / 2 ^ 80 := by
  unfold birthSmokeZ22
  norm_num [zDen]

set_option maxHeartbeats 0 in
theorem zChecked : sparseScaledPowerCheck22 zNum 80 52428 ztrace = true := by decide +kernel

def binomQ (q : ℕ) : ℤ := (Nat.factorial 2048 / (Nat.factorial q * Nat.factorial (2048 - q)) : ℕ)

theorem binomQ_eq (q : ℕ) (hq : q ≤ 2048) : (binomQ q : ℝ) = (Nat.choose 2048 q : ℝ) := by
  unfold binomQ
  rw [← Nat.choose_eq_factorial_div_factorial hq]
  simp

theorem binomQ_pos (q : ℕ) (hq : q ≤ 2048) : 0 < binomQ q := by
  have h : (0 : ℝ) < binomQ q := by rw [binomQ_eq q hq]; exact_mod_cast Nat.choose_pos hq
  exact_mod_cast h

end {ns}.Mk
''')
    tags = {j: (f'm{j}' if j >= 0 else f'n{-j}') for j in certs}
    imports = '\n'.join(f'import {pre}Mk{t}' for t in tags.values())
    terms = []
    for q in range(qlo, qhi + 1):
        j = choice[q]
        c, tag = certs[j], tags[j]
        terms.append(f'''set_option maxHeartbeats 0 in
theorem mterm_{q} :
    (Nat.choose 2048 {q} : ℝ) * ((turboOldUniformFactor22 ^ {q} / birthSmokeZ22 ^ 52428) ^ (1 : ℝ) *
      ∑ k : Fin 11, ((pairedRegionDP22 selectedRows 1 64 {q}) ^ 64) 0 k) ≤ (({recs[q]} : ℤ) : ℝ) / 2 ^ {KT} :=
  TailMarker22.tail_marker_le selectedRows selectedRows_nonneg {q} (binomQ {q}) (binomQ_eq {q} (by norm_num))
    (binomQ_pos {q} (by norm_num)) row {den} (by decide +kernel) hrow {c['tn0']} {c['td0']} (by norm_num)
    (by norm_num) {c['kp']} {c['kn']} (2 ^ {S}) {S} rfl {c['down']} {c['up']} U0_{tag} base_{tag} st_{tag}
    rfl chain_{tag} birthSmokeZ22 zNum hz ztrace zChecked {recs[q]} {KT} (by decide +kernel)
''')
    amatch = '\n'.join(f'  | {q} => {recs[q]}' for q in range(qlo, qhi + 1))
    # A literal match with hundreds of arms makes `exact mterm_q` unfold deep
    # Nat.decEq chains (recursion limit); marker bounds are uniformly 1.
    adef = ('def aMk (_q : ℕ) : Int := 1' if all(recs[q] == 1 for q in range(qlo, qhi + 1))
            else 'def aMk (q : ℕ) : Int :=\n  match q with\n' + amatch + '\n  | _ => 0')
    cases = '\n'.join(f'  · exact mterm_{q}' for q in range(qlo, qhi + 1))
    w('MkTerm', f'''{imports}

noncomputable section
namespace {ns}.Mk
open Spin.Turbo Spin.Turbo.Hi

{chr(10).join(terms)}
{adef}

theorem mkTerm (q : ℕ) (hq : q ∈ Finset.Icc {qlo} {qhi}) :
    (Nat.choose 2048 q : ℝ) * ((turboOldUniformFactor22 ^ q / birthSmokeZ22 ^ 52428) ^ (1 : ℝ) *
      ∑ k : Fin 11, ((pairedRegionDP22 selectedRows 1 64 q) ^ 64) 0 k) ≤ (aMk q : ℝ) / 2 ^ {KT} := by
  simp only [Finset.mem_Icc] at hq
  obtain ⟨h1, h2⟩ := hq
  interval_cases q
{cases}

#print axioms mkTerm
end {ns}.Mk
''')


if __name__ == '__main__':
    main()
