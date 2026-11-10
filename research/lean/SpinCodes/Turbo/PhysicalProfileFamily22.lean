import SpinCodes.Turbo.PhysicalProfileFamilyDefs22
import SpinCodes.Turbo.PhysicalProfileNumericCache22
import SpinCodes.Turbo.PhysicalProfileNumericTable22

noncomputable section
namespace Spin.Turbo
open Spin.Numeric Spin.Numeric.Fix Spin.UniformMixer

theorem profileCachedCapsCheck22_entry (entries : List (Nat × List Nat)) (cache : ProfileFixCache22)
    (bank : List (Fix × Fix)) (j : Nat) (caps : List Int) (den : Int)
    (h : profileCachedCapsCheck22 entries cache bank j caps den = true) (i : Fin 11) :
    (if i.val = 2 then profileFixUpperCheck22
      (profileFixClassMean22 entries i.val (fun p => profileFixCachedMoment22 cache bank p j))
      (caps.getD i.val 0) den
    else profileFixClassUpperCheck22 entries i.val
      (fun p => profileFixCachedMoment22 cache bank p j) (caps.getD i.val 0) den) = true := by
  simp only [profileCachedCapsCheck22, Bool.and_eq_true] at h
  exact List.all_eq_true.mp h.2 i.val (List.mem_range.mpr i.isLt)

theorem profileCachedCapsCheck22_actual (caps : List Int) {den : Int} (hd : 0 < den)
    {box : Fix} {z : ℝ} (hz : Mem box z) (j : Fin 17)
    (h : profileCachedCapsCheck22 profileNumericEntries22
      (profileFixCache22 profileNumericEntries22 (profileFixMomentBank22 box))
      (profileFixMomentBank22 box) j.val caps den = true)
    (i : Fin 11) (hi : i ≠ 2) (q : State22) (hq : physicalClassAllows22 i q) :
    profileNumericMoment22 (List.ofFn (physicalProfile22 q)) z j ≤
      (caps.getD i.val 0 : ℝ)/(den : ℝ) := by
  have hh := profileCachedCapsCheck22_entry _ _ _ _ _ _ h i
  have hi' : i.val ≠ 2 := by intro he; exact hi (Fin.ext he)
  simp only [hi', ite_false] at hh
  exact profileNumericClassUpperCheck22_actual i _ _
    (fun p => profileFixCachedMoment22_sound profileNumericEntries22 p hz j) hd hh q hq

theorem profileCachedCapsCheck22_mean (caps : List Int) {den : Int} (hd : 0 < den)
    {box : Fix} {z : ℝ} (hz : Mem box z) (j : Fin 17)
    (h : profileCachedCapsCheck22 profileNumericEntries22
      (profileFixCache22 profileNumericEntries22 (profileFixMomentBank22 box))
      (profileFixMomentBank22 box) j.val caps den = true) :
    physicalProfileClassMean22 2 (fun p => profileNumericMoment22 (List.ofFn p) z j) ≤
      (caps.getD 2 0 : ℝ)/(den : ℝ) := by
  have hh := profileCachedCapsCheck22_entry _ _ _ _ _ _ h (2 : Fin 11)
  change profileFixUpperCheck22 (profileFixClassMean22 profileNumericEntries22 2 (fun p => profileFixCachedMoment22 (profileFixCache22 profileNumericEntries22 (profileFixMomentBank22 box)) (profileFixMomentBank22 box) p j.val)) (caps.getD 2 0) den = true at hh
  exact profileFixUpperCheck22_sound
    (profileNumericClassMean22_actual_mem 2 _ _
      (fun p => profileFixCachedMoment22_sound profileNumericEntries22 p hz j)) hd hh

theorem profileReturnCap_common_den_sound (mn sn fn rn den live : Int)
    (hd : 0 < den) (hl : 0 < live) (hr : 0 ≤ rn)
    (h : mn ≤ rn*live ∨ sn*fn ≤ (rn*live)^2) :
    min ((mn : ℝ)/(den : ℝ)) (Real.sqrt (((sn : ℝ)/(den : ℝ))*((fn : ℝ)/(den : ℝ)))) /
      (live : ℝ) ≤ (rn : ℝ)/(den : ℝ) := by
  have hdR : (0 : ℝ) < (den : ℝ) := by exact_mod_cast hd
  have hlR : (0 : ℝ) < (live : ℝ) := by exact_mod_cast hl
  apply (div_le_iff₀ hlR).mpr
  rcases h with hm | hs
  · apply (min_le_left _ _).trans
    have hh : (mn : ℝ) ≤ (rn : ℝ)*(live : ℝ) := by exact_mod_cast hm
    convert div_le_div_of_nonneg_right hh hdR.le using 1 <;> ring
  · apply (min_le_right _ _).trans
    apply Real.sqrt_le_iff.mpr
    refine ⟨mul_nonneg (div_nonneg (by exact_mod_cast hr) hdR.le) hlR.le, ?_⟩
    have hh : (sn : ℝ)*(fn : ℝ) ≤ ((rn : ℝ)*(live : ℝ))^2 := by exact_mod_cast hs
    have hd2 : (0 : ℝ) ≤ (den : ℝ)^2 := sq_nonneg _
    convert div_le_div_of_nonneg_right hh hd2 using 1 <;>
      field_simp [hdR.ne'] <;> ring

theorem profileReturnCapsCheck22_sound (first second returns : List Int) (active den : Int)
    (hd : 0 < den) (h : profileReturnCapsCheck22 first second returns active = true)
    (i : Fin 11) (hi : i ≠ 0) :
    min ((first.getD i.val 0 : ℝ)/(den : ℝ))
      (Real.sqrt (((second.getD i.val 0 : ℝ)/(den : ℝ))*((active : ℝ)/(den : ℝ)))) /
        (4194303 : ℝ) ≤ (returns.getD i.val 0 : ℝ)/(den : ℝ) := by
  simp only [profileReturnCapsCheck22, Bool.and_eq_true, decide_eq_true_eq] at h
  have hidx : i.val < returns.length := by rw [h.1.1]; exact i.isLt
  have hm : returns.getD i.val 0 ∈ returns := by
    rw [List.getD_eq_getElem _ _ hidx]
    exact List.getElem_mem _
  have hr : 0 ≤ returns.getD i.val 0 := of_decide_eq_true (List.all_eq_true.mp h.1.2 _ hm)
  have hival : i.val ≠ 0 := by intro he; exact hi (Fin.ext he)
  have hk : i.val-1 < 10 := by omega
  have he : i.val-1+1 = i.val := by omega
  have hh := List.all_eq_true.mp h.2 (i.val-1) (List.mem_range.mpr hk)
  simp only [he, decide_eq_true_eq] at hh
  apply profileReturnCap_common_den_sound _ _ _ _ den 4194303 hd (by decide) hr
  simpa only [pow_two] using hh

#print axioms profileCachedCapsCheck22_actual
#print axioms profileCachedCapsCheck22_mean
#print axioms profileReturnCapsCheck22_sound
end Spin.Turbo

