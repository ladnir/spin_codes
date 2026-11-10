"""Emit per-family semantic tail modules Tail{k}ExTerm.lean (text only; Lean checks)."""
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / 'SpinCodes/Turbo/Generated/K18Tail'
S, MULT_BITS, KT = 704, 256, 256
FAMILIES = {1: (3, 5), 2: (6, 10), 3: (11, 19), 4: (20, 34), 5: (35, 51), 6: (52, 61),
            8: (144, 157), 9: (158, 185), 10: (186, 224), 11: (225, 279)}


def emit(k):
    qlo, qhi = FAMILIES[k]
    D = qhi
    ns = f'Spin.Turbo.Hi.ProfileFamilyK18Tail{k}22'
    pre = f'SpinCodes.Turbo.Generated.K18Tail.Tail{k}'
    terms = []
    for q in range(qlo, qhi + 1):
        terms.append(f'''set_option maxHeartbeats 0 in
theorem term_{q} :
    (Nat.choose 2048 {q} : ℝ) * ((turboOldUniformFactor22 ^ {q} / birthSmokeZ22 ^ 52428) ^ (1 : ℝ) *
      ∑ k : Fin 11, ((pairedRegionDP22 selectedRows 1 64 {q}) ^ 64) 0 k) ≤ (a_{q} : ℝ) / 2 ^ {KT} :=
  TailTerm22.tail_term_le selectedRows selectedRows_nonneg {D} {q} (by norm_num) (2 ^ {S}) {S} rfl
    (sparseRegionalLastL22 base trace) hlast (binomQ {q}) (binomQ_eq {q} (by norm_num))
    (binomQ_pos {q} (by norm_num)) U0_{q} (by decide +kernel) dn0_{q} up0_{q} (by decide +kernel) st_{q} rfl
    (by decide +kernel) birthSmokeZ22 zNum hz ztrace zChecked a_{q} {KT} (by decide +kernel)
''')
    cases = '\n'.join(f'  · exact term_{q}' for q in range(qlo, qhi + 1))
    amatch = '\n'.join(f'  | {q} => a_{q}' for q in range(qlo, qhi + 1))
    text = f'''import {pre}ExChain
import {pre}ExQData
import {pre}Family
import SpinCodes.Turbo.TailTerm22

/-! Generated semantic tail terms for family {k} (q = {qlo}..{qhi}, α = 1). -/
noncomputable section
namespace {ns}.Ex
open Spin.Turbo Spin.Turbo.Hi

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

theorem hscale : (16 ^ 16 * endpointDen * 2 ^ {MULT_BITS} : ℤ) = 2 ^ {S} := by decide +kernel

theorem hlast : SparseListPolyD22.PolyBoundD selectedRows (2 ^ {S}) 128
    (sparseRegionalLastL22 base trace) {D} := by
  have hrowT : ∀ d (hd : d ≤ 16), ∀ i k : Fin 11, selectedRows ⟨d, by omega⟩ i k ≤
      (rowT d i.val k.val : ℝ) / ((16 ^ 16 * endpointDen : ℤ) : ℝ) := by
    intro d hd i k
    rw [rowT_eq d hd i.val i.isLt k.val k.isLt]
    exact hrow d hd i k
  have hb := SparseListPolyD22.baseD_sound {D} selectedRows selectedRows_nonneg rowT
    (16 ^ 16 * endpointDen) (2 ^ {MULT_BITS}) (by decide +kernel) (by decide +kernel) hrowT base
    baseChecked
  rw [hscale] at hb
  have hc := SparseListPolyD22.chainD_sound {D} selectedRows selectedRows_nonneg (2 ^ {S})
    (by decide +kernel) base hb trace base 1 hb chainChecked
  have he : trace.foldl sparseRegionalStepL22 1 = 128 := by decide +kernel
  rw [he] at hc
  exact hc

def binomQ (q : ℕ) : ℤ := (Nat.factorial 2048 / (Nat.factorial q * Nat.factorial (2048 - q)) : ℕ)

theorem binomQ_eq (q : ℕ) (hq : q ≤ 2048) : (binomQ q : ℝ) = (Nat.choose 2048 q : ℝ) := by
  unfold binomQ
  rw [← Nat.choose_eq_factorial_div_factorial hq]
  simp

theorem binomQ_pos (q : ℕ) (hq : q ≤ 2048) : 0 < binomQ q := by
  have h : (0 : ℝ) < binomQ q := by rw [binomQ_eq q hq]; exact_mod_cast Nat.choose_pos hq
  exact_mod_cast h

theorem hz : birthSmokeZ22 = (zNum : ℝ) / 2 ^ 80 := by
  unfold birthSmokeZ22
  norm_num [zDen]

set_option maxHeartbeats 0 in
theorem zChecked : sparseScaledPowerCheck22 zNum 80 52428 ztrace = true := by decide +kernel

{chr(10).join(terms)}
def aFam (q : ℕ) : Int :=
  match q with
{amatch}
  | _ => 0

theorem famTerm (q : ℕ) (hq : q ∈ Finset.Icc {qlo} {qhi}) :
    (Nat.choose 2048 q : ℝ) * ((turboOldUniformFactor22 ^ q / birthSmokeZ22 ^ 52428) ^ (1 : ℝ) *
      ∑ k : Fin 11, ((pairedRegionDP22 selectedRows 1 64 q) ^ 64) 0 k) ≤ (aFam q : ℝ) / 2 ^ {KT} := by
  simp only [Finset.mem_Icc] at hq
  obtain ⟨h1, h2⟩ := hq
  interval_cases q
{cases}

#print axioms famTerm
end {ns}.Ex
'''
    (OUT / f'Tail{k}ExTerm.lean').write_text(text, encoding='utf-8', newline='\n')


if __name__ == '__main__':
    for k in map(int, sys.argv[1:] or FAMILIES):
        emit(k)
    print('ok')
