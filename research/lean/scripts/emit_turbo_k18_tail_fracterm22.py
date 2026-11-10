"""Emit Tail7FrRows / Tail7FrTerm (semantic layer for family 7, alpha = 4/5). Text only."""
import json
from pathlib import Path
ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / 'SpinCodes/Turbo/Generated/K18Tail'
K, QLO, QHI, D, S, KT = 7, 62, 143, 143, 704, 256
ns = f'Spin.Turbo.Hi.ProfileFamilyK18Tail{K}22'
pre = f'SpinCodes.Turbo.Generated.K18Tail.Tail{K}'
man = json.loads((OUT / 'frac_manifest.json').read_text())['q']
rows = f'''import {pre}Family
import {pre}FrR
import {pre}FrQData
import SpinCodes.Turbo.Hi.RowSelected22
import SpinCodes.Turbo.TailTermFrac22

noncomputable section
namespace {ns}.Fr
open Spin.Turbo Spin.Turbo.Hi

def births (n : Nat) : List Int :=
  match n with
  | 0 => birthZeroCaps22
  | 1 => birthOneCaps
  | _ => birthTwoCaps

def physInt (j i k : Nat) : Int := profileRowIntAll22 births firstCaps returnCaps zeroHoldCaps j i k

def row (d i k : Nat) : Int := rowSel physInt d i k

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

set_option maxHeartbeats 0 in
theorem hRlit : R = pairedRowMatsL22 row := by decide +kernel

theorem hR (a : Fin 17) (i k : Fin 11) :
    selectedRows a i k ≤ (sparseLMat22 (R.getD a.val []) i.val k.val : ℝ) / ((16 ^ 16 * endpointDen : ℤ) : ℝ) := by
  rw [hRlit]
  exact PairedFracChain22.rowMats_hyp selectedRows row (16 ^ 16 * endpointDen) hrow a i k

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

end {ns}.Fr
'''
(OUT / f'Tail{K}FrRows.lean').write_text(rows, encoding='utf-8', newline='\n')
terms = []
for q in range(QLO, QHI + 1):
    terms.append(f'''set_option maxHeartbeats 0 in
theorem fterm_{q} :
    (Nat.choose 2048 {q} : ℝ) * ((turboOldUniformFactor22 ^ {q} / birthSmokeZ22 ^ 52428) ^ ((4 : ℝ) / 5) *
      ∑ k : Fin 11, ((pairedRegionDP22 selectedRows ((4 : ℝ) / 5) 64 {q}) ^ 64) 0 k) ≤ (a_{q} : ℝ) / 2 ^ {KT} :=
  TailTermFrac22.tail_term_frac_le selectedRows selectedRows_nonneg {D} {q} (by norm_num) (2 ^ {S}) {S} rfl
    (16 ^ 16 * endpointDen) R hR Ys rpowChecked base baseChecked trace (by decide +kernel) chainChecked
    (binomQ {q}) (binomQ_eq {q} (by norm_num)) (binomQ_pos {q} (by norm_num)) U0_{q} (by decide +kernel)
    (by decide +kernel) st_{q} rfl (by decide +kernel) birthSmokeZ22 zNum hz ztrace zChecked vt_{q}
    (by decide +kernel) yn_{q} G_{q} (by decide +kernel) a_{q} {KT} (by decide +kernel)
''')
amatch = '\n'.join(f'  | {q} => a_{q}' for q in range(QLO, QHI + 1))
cases = '\n'.join(f'  · exact fterm_{q}' for q in range(QLO, QHI + 1))
nl = '\n'
term = f'''import {pre}FrRows
import {pre}FrChain

noncomputable section
namespace {ns}.Fr
open Spin.Turbo Spin.Turbo.Hi

{nl.join(terms)}
def aFr (q : ℕ) : Int :=
  match q with
{amatch}
  | _ => 0

theorem frTerm (q : ℕ) (hq : q ∈ Finset.Icc {QLO} {QHI}) :
    (Nat.choose 2048 q : ℝ) * ((turboOldUniformFactor22 ^ q / birthSmokeZ22 ^ 52428) ^ ((4 : ℝ) / 5) *
      ∑ k : Fin 11, ((pairedRegionDP22 selectedRows ((4 : ℝ) / 5) 64 q) ^ 64) 0 k) ≤ (aFr q : ℝ) / 2 ^ {KT} := by
  simp only [Finset.mem_Icc] at hq
  obtain ⟨h1, h2⟩ := hq
  interval_cases q
{cases}

#print axioms frTerm
end {ns}.Fr
'''
(OUT / f'Tail{K}FrTerm.lean').write_text(term, encoding='utf-8', newline='\n')
print('ok')
