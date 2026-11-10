import SpinCodes.Turbo.PhysicalProfileNumericCacheDefs22
import SpinCodes.Turbo.PhysicalProfileNumericAggregate22

noncomputable section
namespace Spin.Turbo
open Spin.Numeric Spin.Numeric.Fix

theorem profileFixCache22_lookup (entries : List (Nat × List Nat)) (bank : List (Fix × Fix))
    (profile : List Nat) (cs : List Fix)
    (h : (profileFixCache22 entries bank).lookup profile = some cs) :
    cs = profileFixCoefficients22 (profileFixFactors22 profile bank) := by
  induction entries generalizing profile cs with
  | nil => simp [profileFixCache22] at h
  | cons e entries ih =>
    by_cases hp : profile = e.2
    · simp only [profileFixCache22, List.map_cons, List.lookup_cons, hp, beq_self_eq_true, ite_true,
        Option.some.injEq] at h
      simpa only [hp] using h.symm
    · have hb : (profile == e.2) = false := beq_eq_false_iff_ne.mpr hp
      simp only [profileFixCache22, List.map_cons, List.lookup_cons, hb, Bool.false_eq_true, ite_false] at h
      exact ih profile cs h

theorem profileFixCachedRaw22_eq (entries : List (Nat × List Nat)) (bank : List (Fix × Fix))
    (profile : List Nat) (j : Nat) :
    profileFixCachedRaw22 (profileFixCache22 entries bank) bank profile j =
      (profileFixCoefficients22 (profileFixFactors22 profile bank)).getD j Fix.zero := by
  unfold profileFixCachedRaw22
  cases hl : (profileFixCache22 entries bank).lookup profile with
  | none => rfl
  | some cs => rw [profileFixCache22_lookup entries bank profile cs hl]

theorem profileFixCachedMoment22_sound (entries : List (Nat × List Nat)) (profile : List Nat)
    {box : Fix} {z : ℝ} (hz : Mem box z) (j : Fin 17) :
    Mem (profileFixCachedMoment22 (profileFixCache22 entries (profileFixMomentBank22 box))
      (profileFixMomentBank22 box) profile j.val) (profileNumericMoment22 profile z j) := by
  unfold profileFixCachedMoment22
  rw [profileFixCachedRaw22_eq]
  exact profileFixMoment22_sound profile hz j

theorem profileFixCachedAllFourier22_sound (entries : List (Nat × List Nat))
    {box : Fix} {z : ℝ} (hz : Mem box z) (j : Fin 17) :
    Mem (profileFixCachedAllFourier22 entries (profileFixCache22 entries (profileFixFourierBank22 box))
      (profileFixFourierBank22 box) j.val)
      (((entries.map (fun e => (e.1 : ℝ) *
        (profileBankPolynomial22 e.2 (profileRealFourierBank22 z)).coeff j.val)).sum) /
          ((4194304 : ℝ) * ((Nat.choose 16 j.val : ℝ)*15^j.val))) := by
  unfold profileFixCachedAllFourier22
  simp_rw [profileFixCachedRaw22_eq]
  exact profileFixAllFourier22_mem entries hz j

#print axioms profileFixCachedMoment22_sound
#print axioms profileFixCachedAllFourier22_sound
end Spin.Turbo
