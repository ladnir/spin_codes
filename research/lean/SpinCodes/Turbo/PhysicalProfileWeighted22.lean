import SpinCodes.Turbo.PhysicalProfileNumericTable22
import SpinCodes.Turbo.PhysicalCauchy22

/-!
# Profile bounds for weighted physical pieces

The historical live-row proof uses a maximum for arbitrary supported mass
and a mean for the uniformly capped piece. This adapter preserves that
distinction when the observable is supplied by the authenticated profile
table. In particular, it does not replace the uniform mean by a maximum.
-/
noncomputable section
namespace Spin.Turbo
open Finset Spin.UniformMixer

theorem physicalProfile22_uniform_mean (f : (Fin 5 → Nat) → ℝ) :
    (∑ q ∈ liveStateSet, f (physicalProfile22 q)) / (liveCount State22 : ℝ) =
      physicalProfileClassMean22 2 f := by
  unfold physicalProfileClassMean22
  rw [← physicalProfile22_class_sum]
  have hc : physicalProfileClassCount22 2 = liveCount State22 := by
    rw [physicalProfileClassCount22_checked, state22_liveCount]
    rfl
  rw [hc]
  apply congrArg (fun x : ℝ => x / (liveCount State22 : ℝ))
  have ha (q : State22) : physicalClassAllows22 2 q ↔ q ≠ 0 := by rfl
  simp only [ha]
  classical
  have hs : (liveStateSet : Finset State22) = univ.filter (fun q => q ≠ 0) := by
    ext q
    simp [liveStateSet]
  rw [hs, Finset.sum_filter]

/-- A class cap suffices except for coordinate two, where the actual live
mean exploits the pointwise uniform cap in `physicalPiece22`. -/
theorem physicalProfile22_weighted_bound (i : Fin 11) (hi : i ≠ 0)
    (b : ℝ) (μ : State22 → ℝ) (hμ : physicalPiece22 i b μ)
    (f : (Fin 5 → Nat) → ℝ) (hf : ∀ q, 0 ≤ f (physicalProfile22 q))
    (M : ℝ) (hM : 0 ≤ M)
    (hcap : i ≠ 2 → ∀ q, physicalClassAllows22 i q → f (physicalProfile22 q) ≤ M)
    (hmean : i = 2 → physicalProfileClassMean22 2 f ≤ M) :
    (∑ q : State22, μ q * f (physicalProfile22 q)) ≤ b*M := by
  by_cases hi2 : i = 2
  · subst i
    apply physicalUniformPiece22_mean_bound hμ _ hf M
    rw [physicalProfile22_uniform_mean]
    exact hmean rfl
  · exact physicalLivePiece22_max_bound hi hμ _ M hM
      (fun q hq => hcap hi2 q (hμ.2.1 q hq))

#print axioms physicalProfile22_uniform_mean
#print axioms physicalProfile22_weighted_bound
end Spin.Turbo
