import SpinCodes.Turbo.Hi.PhysicalProfileAssembly22
import SpinCodes.Turbo.Hi.PhysicalProfileFamilyExtract22
import SpinCodes.Turbo.Hi.PhysicalProfileFamilyFourier22

/-!
# Complete checked profile families for physical and selected steps

The family Boolean authenticates all 17 occupations. The small-birth census
and coarse birth bound provide the remaining birth cones. The selected-step
statement then uses the already proved binomial thinning of the actual input
law. No numerical moment or feedback-probability hypothesis remains.
-/
noncomputable section
namespace Spin.Turbo.Hi
open Spin.NumericHi Spin.NumericHi.Fix Spin.UniformMixer

def profileFamilyRows22 (birth : Fin 17 → Fin 11 → ℝ)
    (first returns : List (List Int)) (den : Int) (j : Fin 17) :
    Matrix (Fin 11) (Fin 11) ℝ :=
  profileOutwardRows22 (birth j)
    (profileEndpoint22 (first.getD j.val []) den)
    (profileEndpoint22 (returns.getD j.val []) den)

theorem profileFamilyRows22_sound {box : Fix} {z : ℝ}
    (hzbox : Mem box z) (hz : 0 ≤ z)
    (first second returns : List (List Int)) (active zeroHold : List Int) (den : Int)
    (h : profileFamilyCheck22 profileNumericEntries22 box first second returns active zeroHold den = true)
    (birth : Fin 17 → Fin 11 → ℝ)
    (hbirth : ∀ j, physicalCone22 (birth j) (turboExactBirthMass22 j z)) (j : Fin 17) :
    physicalRowSound22 (Matrix.of (turboPhysicalKernel22 j z))
      (profileFamilyRows22 birth first returns den j) := by
  have hj := profileFamilyCheck22_row _ _ _ _ _ _ _ _ h j
  obtain ⟨hs, hfirst, hsecond, hfeedback, hzero, hreturn⟩ :=
    profileFamilyRowCheck22_parts _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ hj
  have hone : Mem (Fix.ofInt 1) (1 : ℝ) := by
    simpa only [Int.cast_one] using ofInt_mem 1
  have hlow := physicalZeroHold22_cached_lower hone j (den-active.getD j.val 0) den hs.1 hfeedback
  have hd : (den : ℝ) ≠ 0 := ne_of_gt (by exact_mod_cast hs.1)
  have hfrac : ((den-active.getD j.val 0 : Int) : ℝ)/(den : ℝ) =
      1 - (active.getD j.val 0 : ℝ)/(den : ℝ) := by
    push_cast
    rw [sub_div, div_self hd]
  rw [hfrac] at hlow
  have hf : physicalFeedbackProbability22 (packetStepLaw22 j) ≤
      (active.getD j.val 0 : ℝ)/(den : ℝ) := by
    rw [physicalFeedbackProbability22_zero_hold]
    linarith
  exact profileCheckedRows22_sound j hzbox hz _ _ _ _ _ hs.1
    hfirst hsecond hreturn hf (birth j) (hbirth j)

theorem profileFamilySelectedRows22_sound {box : Fix} {z : ℝ}
    (hzbox : Mem box z) (hz : 0 ≤ z)
    (first second returns : List (List Int)) (active zeroHold : List Int) (den : Int)
    (h : profileFamilyCheck22 profileNumericEntries22 box first second returns active zeroHold den = true)
    (birth : Fin 17 → Fin 11 → ℝ)
    (hbirth : ∀ j, physicalCone22 (birth j) (turboExactBirthMass22 j z)) (j : Fin 17) :
    physicalRowSound22 (Matrix.of (turboSelectedKernel22 j z))
      (turboThinnedRow22 j z (profileFamilyRows22 birth first returns den)) :=
  turboThinnedRow22_sound j z _
    (profileFamilyRows22_sound hzbox hz first second returns active zeroHold den h birth hbirth)

private theorem endpoint_nonneg22 (caps : List Int) (den : Int) (hd : 0 < den)
    (hlen : caps.length = 11) (hn : caps.all (fun n => decide (0 ≤ n)) = true)
    (i : Fin 11) : 0 ≤ profileEndpoint22 caps den i := by
  have hi : i.val < caps.length := by rw [hlen]; exact i.isLt
  have hm : caps.getD i.val 0 ∈ caps := by
    rw [List.getD_eq_getElem _ _ hi]
    exact List.getElem_mem _
  have hh : 0 ≤ caps.getD i.val 0 := of_decide_eq_true (List.all_eq_true.mp hn _ hm)
  exact div_nonneg (by exact_mod_cast hh) (by exact_mod_cast hd.le)

theorem profileFamilyRows22_nonneg {box : Fix} {z : ℝ}
    (first second returns : List (List Int)) (active zeroHold : List Int) (den : Int)
    (h : profileFamilyCheck22 profileNumericEntries22 box first second returns active zeroHold den = true)
    (birth : Fin 17 → Fin 11 → ℝ)
    (hbirth : ∀ j, physicalCone22 (birth j) (turboExactBirthMass22 j z))
    (j : Fin 17) (i k : Fin 11) :
    0 ≤ profileFamilyRows22 birth first returns den j i k := by
  have hj := profileFamilyCheck22_row _ _ _ _ _ _ _ _ h j
  obtain ⟨hs, hf, _, _, _, hr⟩ :=
    profileFamilyRowCheck22_parts _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ hj
  simp only [profileCachedCapsCheck22, Bool.and_eq_true, decide_eq_true_eq] at hf
  simp only [profileReturnCapsCheck22, Bool.and_eq_true, decide_eq_true_eq] at hr
  have hfirst := endpoint_nonneg22 _ den hs.1 hf.1.1 hf.1.2 i
  have hreturn := endpoint_nonneg22 _ den hs.1 hr.1.1 hr.1.2 i
  have hb : 0 ≤ birth j k := by
    obtain ⟨parts, _, hp⟩ := hbirth j
    exact (Finset.sum_nonneg (fun q _ => (hp k).1 q)).trans (hp k).2.2.1
  change 0 ≤ if i = 0 then birth j k else
    if k = 0 then profileEndpoint22 (returns.getD j.val []) den i else
    if k = 2 then profileEndpoint22 (first.getD j.val []) den i else 0
  split_ifs
  · exact hb
  · exact hreturn
  · exact hfirst
  · exact le_rfl

theorem profileFamilySelectedRows22_nonneg {box : Fix} {z : ℝ}
    (first second returns : List (List Int)) (active zeroHold : List Int) (den : Int)
    (h : profileFamilyCheck22 profileNumericEntries22 box first second returns active zeroHold den = true)
    (birth : Fin 17 → Fin 11 → ℝ)
    (hbirth : ∀ j, physicalCone22 (birth j) (turboExactBirthMass22 j z))
    (j : Fin 17) (i k : Fin 11) :
    0 ≤ turboThinnedRow22 j z (profileFamilyRows22 birth first returns den) i k := by
  change 0 ≤ ∑ ell : Fin 17,
    thinningWeight22 j ell * profileFamilyRows22 birth first returns den ell i k
  exact Finset.sum_nonneg (fun ell _ => mul_nonneg (thinningWeight22_nonneg j ell)
    (profileFamilyRows22_nonneg first second returns active zeroHold den h birth hbirth ell i k))

#print axioms profileFamilyRows22_sound
#print axioms profileFamilySelectedRows22_sound
#print axioms profileFamilyRows22_nonneg
#print axioms profileFamilySelectedRows22_nonneg
end Spin.Turbo.Hi
