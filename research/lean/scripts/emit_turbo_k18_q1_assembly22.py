"""Emit the per-tilt semantic assembly for the generated K18 q1 certificates.

Untrusted text generation only: every emitted statement is checked by Lean.
"""
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / 'SpinCodes/Turbo/Generated/K18'
D = 52428
R = 128


def emit(k: int) -> None:
    ns = f'Spin.Turbo.ProfileFamilyK18Tilt{k}22'
    imports = '\n'.join(f'import SpinCodes.Turbo.Generated.K18.Tilt{k}Steps{c}' for c in range(8))
    rows_cases = '\n'.join(f'  | {n} => rows{n}' for n in range(65))
    step_cases = '\n'.join(f'  · exact step{n}' for n in range(64))
    text = f'''{imports}
import SpinCodes.Turbo.PhysicalProfileFamilyK18Tilt{k}22
import SpinCodes.Turbo.SparseListQ1Regional22

/-! Generated semantic assembly (tilt {k}): integer certificates ⇒ bounds on
the actual normalized regional kernels, every q1 shell coefficient, and z^d. -/
noncomputable section
namespace {ns}.Cert
open Spin.Turbo Spin.Turbo.SparseListSound22 Spin.Turbo.SparseListQ1Regional22

def rows : ℕ → List (List Int)
{rows_cases}
  | _ => []

theorem hstep : ∀ n, n < 64 →
    sparseQ1StepCheckL22 s res0 res1 (rows n) (rows (n + 1)) = true := by
  intro n hn
  interval_cases n
{step_cases}

def birthsF (j : Fin 3) : List Int := births j.val

theorem hsmall (j : Fin 3) (i : Fin 11) :
    smallBirthSmokeCaps22 j i = ((birthsF j).getD i.val 0 : ℝ) / endpointDen := by
  fin_cases j <;> simp [smallBirthSmokeCaps22, birthsF, births]

theorem rowsInt (d : ℕ) (hd : d ≤ 2) (i k : Fin 11) :
    physicalRows ⟨d, by omega⟩ i k = (row d i.val k.val : ℝ) / endpointDen :=
  SparseRowBase22.profileFamilyRows22_small_eq smallBirthSmokeCaps22 firstCaps returnCaps
    zeroHoldCaps endpointDen birthsF hsmall d (by omega) i k

theorem regional :
    (∀ d, d < 3 → MatShape (result.getD d [])) ∧
      ∀ d, d ≤ 2 → ∀ i j : Fin 11,
        sparseRegionalDP22 physicalRows {R} d i j ≤ (sparseLPoly22 result d i j : ℝ) / s :=
  certificate_sound physicalRows physicalRows_nonneg {R} s endpointDen mult (by decide +kernel)
    (by decide +kernel) row rowsInt base trace result regionalChecked

theorem regionalChoice_nonneg (b : Bool) (i j : Fin 11) :
    0 ≤ sparseOneRegionChoice22 physicalRows {R} b i j := by
  cases b <;> exact sparseRegionalDP22_nonneg physicalRows physicalRows_nonneg {R} _ i j

theorem q1 : ∀ v, v < 65 → ∀ j : Fin 11,
    orderedMatrixDP sparseBoolDegree22 (sparseOneRegionChoice22 physicalRows {R}) 64 v 0 j ≤
      (sparseLMat22 rows64 v j : ℝ) / s :=
  (q1Rows_sound s (by decide +kernel) res0 res1 (regional.1 0 (by norm_num)) (regional.1 1 (by norm_num))
    (sparseOneRegionChoice22 physicalRows {R}) regionalChoice_nonneg
    (fun i j => regional.2 0 (by norm_num) i j) (fun i j => regional.2 1 (by norm_num) i j)
    rows 64 rfl hstep 64 le_rfl).2

/-- Every q1 shell coefficient mass is bounded by the materialized terminal row. -/
theorem shell_le (v : ℕ) (hv : v < 65) :
    ∑ j : Fin 11, orderedMatrixCoefficient sparseBoolDegree22
        (sparseOneRegionChoice22 physicalRows {R}) 64 v 0 j ≤
      (∑ j : Fin 11, (sparseLMat22 rows64 v j : ℝ)) / s := by
  rw [orderedMatrixCoefficient_eq_dp, Finset.sum_div]
  exact Finset.sum_le_sum fun j _ => q1 v hv j

theorem z_eq : (zbaseInt : ℝ) / s = birthSmokeZ22 := by
  unfold birthSmokeZ22
  rw [div_eq_div_iff (by norm_num [s]) (by norm_num [zDen])]
  norm_num [zbaseInt, s, zNum, zDen]

/-- Strictly positive lower bound for z^d with the family's exact rational z. -/
theorem zpow_lower :
    0 < (sparsePowerTraceLast22 s ztrace : ℝ) / s ∧
      (sparsePowerTraceLast22 s ztrace : ℝ) / s ≤ birthSmokeZ22 ^ {D} := by
  have h := SparseFixedSound22.powerCertificate_sound s zbaseInt {D} ztrace powerChecked
  rw [z_eq] at h
  exact h

#print axioms shell_le
#print axioms zpow_lower
end {ns}.Cert
'''
    (OUT / f'Tilt{k}Q1.lean').write_text(text, encoding='utf-8', newline='\n')


if __name__ == '__main__':
    for k in (1, 2, 3):
        emit(k)
    print('emitted')
