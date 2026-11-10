"""Emit the K20 q1 endpoint module from the certificate manifest (untrusted).

Text generation only; every statement is checked by Lean.
"""
import json
import sys
from fractions import Fraction as Q
from math import comb, log2
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / 'scripts'))
import emit_turbo_k20_sparse_cert22 as P  # noqa: E402

OUT = ROOT / 'SpinCodes/Turbo/Generated/K20'
TEMPLATE = ROOT / 'scripts/templates/k20_q1_endpoint.lean.tmpl'
K = 256
M8 = ((1 << 32) - 1) ** 8
T = {k: f'Spin.Turbo.ProfileFamilyK20Tilt{k}22' for k in (1, 2, 3, 4)}
CONDS = {1: '8 ≤ v ∧ v ≤ 13', 2: '14 ≤ v ∧ v ≤ 21', 3: '22 ≤ v ∧ v ≤ 57', 4: '58 ≤ v'}


def tilt_of(v):
    return 1 if v <= 13 else (2 if v <= 21 else (3 if v <= 57 else 4))


def check_block(k):
    return f'''def lastT{k} : Int := sparsePowerTraceLast22 {T[k]}.Cert.s {T[k]}.Cert.ztrace

def sumRowT{k} (v : ℕ) : Int := sparseSum11Int22 (fun j => sparseLMat22 {T[k]}.Cert.rows64 v j)

def checkT{k} : Bool := (List.range 65).all fun v =>
  if {CONDS[k]} then decide (cnt v * sumRowT{k} v * 2 ^ {K} ≤ aNum v * (M8 * binom v * lastT{k}))
  else true

theorem checkT{k}_ok : checkT{k} = true := by decide +kernel
'''


def tilt_block(k):
    return f'''theorem tiltTerm{k} (points : Fin 8 → Field16) (hp : Function.Injective points)
    (v : ℕ) (hv : v < 65) (hc : {CONDS[k]}) :
    (rs8ExpectedPacketCount points v / (Nat.choose 64 v : ℝ)) / {T[k]}.birthSmokeZ22 ^ {P.D} *
      ∑ k : Fin 11, orderedMatrixCoefficient sparseBoolDegree22
        (sparseOneRegionChoice22 {T[k]}.physicalRows {P.R}) 64 v 0 k ≤
      (aNum v : ℝ) / 2 ^ {K} := by
  have hchk := List.all_eq_true.mp checkT{k}_ok v (List.mem_range.mpr hv)
  rw [if_pos hc] at hchk
  have hz := {T[k]}.Cert.zpow_lower
  have hs : (0 : ℝ) < ({T[k]}.Cert.s : ℝ) := by
    exact_mod_cast (show (0 : ℤ) < {T[k]}.Cert.s by decide +kernel)
  have hcoef := {T[k]}.Cert.shell_le v hv
  have hsum : (∑ j : Fin 11, (sparseLMat22 {T[k]}.Cert.rows64 v j : ℝ)) = (sumRowT{k} v : ℝ) := by
    unfold sumRowT{k}
    rw [SparseFixedSound22.sum11_cast]
  rw [hsum] at hcoef
  exact shell_term _ _ _ _ (cnt v) (sumRowT{k} v) (binom v) lastT{k} {T[k]}.Cert.s (aNum v)
    (count_cast points hp v) (binom_eq v (by omega)) (binom_pos v (by omega)) (cnt_nonneg v)
    hs hz.1 hz.2
    (coef_nonneg _ (fun b i j => by
      cases b <;> exact sparseRegionalDP22_nonneg _ {T[k]}.physicalRows_nonneg _ _ i j) v)
    hcoef (of_decide_eq_true hchk)
'''


def zpos_block():
    out = []
    for k in (1, 2, 3, 4):
        out.append(f'''theorem zpos{k} : 0 < {T[k]}.birthSmokeZ22 := by
  unfold {T[k]}.birthSmokeZ22
  exact div_pos (by exact_mod_cast (show (0 : ℤ) < {T[k]}.zNum by decide +kernel))
    (by exact_mod_cast (show (0 : ℤ) < {T[k]}.zDen by decide +kernel))

theorem zle{k} : {T[k]}.birthSmokeZ22 ≤ 1 := by
  unfold {T[k]}.birthSmokeZ22
  rw [div_le_one (by exact_mod_cast (show (0 : ℤ) < {T[k]}.zDen by decide +kernel))]
  exact_mod_cast (show {T[k]}.zNum ≤ {T[k]}.zDen by decide +kernel)
''')
    return '\n'.join(out)


def main():
    man = json.loads((OUT / 'manifest.json').read_text())
    cnt, den = P.rs8_counts()
    assert den == M8
    a = []
    for v in range(65):
        if v < 8:
            a.append(-((-cnt[v] * (1 << K)) // M8))
        else:
            t = man['tilts'][str(tilt_of(v))]
            num = cnt[v] * t['terminal'][v] * (1 << K)
            dd = M8 * comb(64, v) * t['zlast']
            a.append(-((-num) // dd))
    q = Q(8192 * sum(a), 1 << K)
    print('q1 bits (rounded numerators):', -log2(q.numerator) + log2(q.denominator))
    text = TEMPLATE.read_text(encoding='utf-8')
    subs = {'@ALIT@': '[' + ','.join(map(str, a)) + ']', '@ASUM@': str(sum(a)),
            '@CHECKS@': '\n'.join(check_block(k) for k in (1, 2, 3, 4)),
            '@TILTS@': '\n'.join(tilt_block(k) for k in (1, 2, 3, 4)),
            '@ZPOS@': zpos_block(), '@K@': str(K), '@D@': str(P.D), '@R@': str(P.R),
            '@T1@': T[1], '@T2@': T[2], '@T3@': T[3], '@T4@': T[4]}
    for key, val in subs.items():
        text = text.replace(key, val)
    assert '@' not in text.replace('@[', '')
    (OUT / 'Q1Endpoint.lean').write_text(text, encoding='utf-8', newline='\n')
    man['q1_endpoint'] = dict(K=K, aSum=sum(a))
    (OUT / 'manifest.json').write_text(json.dumps(man, indent=1, default=str) + '\n')


if __name__ == '__main__':
    main()
