import SpinCodes.Turbo.PhysicalProfileRows22
import SpinCodes.Turbo.PhysicalProfileFamily22

/-!
# Checked profile caps imply sound physical transfer rows

The historical `turboPhysicalKernel22_outward_cauchy_rows_sound` supplies the
physical kernel argument. Here its moment hypotheses are discharged by the
integer profile checks, with a mean for the uniformly dominated class and
maxima for the other classes. The return cap is already divided by the number
of nonzero states. Birth and feedback bounds are supplied by their separate
authenticated calculations; this theorem alone is not a finite certificate.
-/
noncomputable section
namespace Spin.Turbo
open Spin.Numeric Spin.Numeric.Fix Spin.UniformMixer

def profileEndpoint22 (caps : List Int) (den : Int) (i : Fin 11) : ℝ :=
  (caps.getD i.val 0 : ℝ) / (den : ℝ)

def profileOutwardRows22 (birth first returns : Fin 11 → ℝ) :
    Matrix (Fin 11) (Fin 11) ℝ := Matrix.of fun i k =>
  if i = 0 then birth k else if k = 0 then returns i else if k = 2 then first i else 0

private theorem cachedCaps_nonneg22 (entries : List (Nat × List Nat))
    (cache : ProfileFixCache22) (bank : List (Fix × Fix)) (j : Nat)
    (caps : List Int) (den : Int) (hd : 0 < den)
    (h : profileCachedCapsCheck22 entries cache bank j caps den = true) (i : Fin 11) :
    0 ≤ profileEndpoint22 caps den i := by
  simp only [profileCachedCapsCheck22, Bool.and_eq_true, decide_eq_true_eq, and_assoc] at h
  have hi : i.val < caps.length := by rw [h.1]; exact i.isLt
  have hm : caps.getD i.val 0 ∈ caps := by
    rw [List.getD_eq_getElem _ _ hi]
    exact List.getElem_mem _
  have hn : 0 ≤ caps.getD i.val 0 :=
    of_decide_eq_true (List.all_eq_true.mp h.2.1 _ hm)
  exact div_nonneg (by exact_mod_cast hn) (by exact_mod_cast hd.le)

theorem profileCheckedRows22_sound (j : Fin 17) {box : Fix} {z : ℝ}
    (hzbox : Mem box z) (hz : 0 ≤ z)
    (first second returns : List Int) (active den : Int) (hd : 0 < den)
    (hfirst : profileCachedCapsCheck22 profileNumericEntries22
      (profileFixCache22 profileNumericEntries22 (profileFixMomentBank22 box))
      (profileFixMomentBank22 box) j.val first den = true)
    (hsecond : profileCachedCapsCheck22 profileNumericEntries22
      (profileFixCache22 profileNumericEntries22 (profileFixMomentBank22 (Fix.mul box box)))
      (profileFixMomentBank22 (Fix.mul box box)) j.val second den = true)
    (hreturn : profileReturnCapsCheck22 first second returns active = true)
    (hfeedback : physicalFeedbackProbability22 (packetStepLaw22 j) ≤
      (active : ℝ)/(den : ℝ))
    (birth : Fin 11 → ℝ) (hbirth : physicalCone22 birth (turboExactBirthMass22 j z)) :
    physicalRowSound22 (Matrix.of (turboPhysicalKernel22 j z))
      (profileOutwardRows22 birth (profileEndpoint22 first den) (profileEndpoint22 returns den)) := by
  have hfirst0 := cachedCaps_nonneg22 _ _ _ _ _ _ hd hfirst
  have hsecond0 := cachedCaps_nonneg22 _ _ _ _ _ _ hd hsecond
  have hsq : Mem (Fix.mul box box) (z*z) := Fix.mul_mem hzbox hzbox
  apply turboPhysicalKernel22_outward_cauchy_rows_sound j z hz birth
    (profileEndpoint22 first den) (profileEndpoint22 second den)
    ((active : ℝ)/(den : ℝ)) hbirth
  · intro i _; exact hsecond0 i
  · intro i hi b μ hμ
    apply physicalProfile22_first_weighted j z hz i hi _ (hfirst0 i) _ _ b μ hμ
    · intro hi2 q hq
      exact profileCachedCapsCheck22_actual first hd hzbox j hfirst i hi2 q hq
    · intro hi2
      subst i
      exact profileCachedCapsCheck22_mean first hd hzbox j hfirst
  · intro i hi b μ hμ
    apply physicalProfile22_second_weighted j z i hi _ (hsecond0 i) _ _ b μ hμ
    · intro hi2 q hq
      exact profileCachedCapsCheck22_actual second hd hsq j hsecond i hi2 q hq
    · intro hi2
      subst i
      exact profileCachedCapsCheck22_mean second hd hsq j hsecond
  · exact hfeedback
  · intro i k
    by_cases hi : i = 0
    · simp only [physicalExactCauchyRows22, profileOutwardRows22, Matrix.of_apply,
        hi, ite_true, le_refl]
    · simp only [physicalExactCauchyRows22, profileOutwardRows22, Matrix.of_apply,
        hi, ite_false, physicalLiveCapRow22]
      by_cases hk : k = 0
      · simp only [hk, ite_true]
        have hl : (liveCount State22 : ℝ) = 4194303 := by
          rw [state22_liveCount]
          norm_num
        rw [hl]
        exact profileReturnCapsCheck22_sound first second returns active den hd hreturn i hi
      · simp only [hk, ite_false, le_refl]

#print axioms profileCheckedRows22_sound
end Spin.Turbo
