import SpinCodes.Structured.ConcretePlacementLimitVanishing
import Mathlib.Analysis.BoxIntegral.UnitPartition
import Mathlib.Analysis.Convex.Measure

/-! The ordered-site simplex is a valid domain for the actual lattice Riemann sum. -/
noncomputable section
namespace Spin.Structured.Placement
open Finset Routing ConcreteEncoder Filter MeasureTheory Bornology
open scoped Topology Pointwise

def orderedSiteDomain (a : Nat) : Set (Fin a → ℝ) :=
  {x | (∀ i, 0 ≤ x i ∧ x i < 1) ∧ StrictMono x}

theorem orderedSiteDomain_convex (a : Nat) : Convex ℝ (orderedSiteDomain a) := by
  intro x hx y hy s t hs ht hst
  constructor
  · intro i
    change 0 ≤ s * x i + t * y i ∧ s * x i + t * y i < 1
    constructor
    · exact add_nonneg (mul_nonneg hs (hx.1 i).1) (mul_nonneg ht (hy.1 i).1)
    · by_cases hsp : 0 < s
      · nlinarith [mul_pos hsp (sub_pos.mpr (hx.1 i).2), mul_nonneg ht (sub_nonneg.mpr (hy.1 i).2.le)]
      · have hs0 : s = 0 := by linarith
        have ht1 : t = 1 := by linarith
        simpa only [hs0, ht1, zero_mul, one_mul, zero_add] using (hy.1 i).2
  · intro i j hij
    change s * x i + t * y i < s * x j + t * y j
    by_cases hsp : 0 < s
    · nlinarith [mul_pos hsp (sub_pos.mpr (hx.2 hij)), mul_nonneg ht (sub_nonneg.mpr (hy.2 hij).le)]
    · have hs0 : s = 0 := by linarith
      have ht1 : t = 1 := by linarith
      simpa only [hs0, ht1, zero_mul, one_mul, zero_add] using hy.2 hij

theorem orderedSiteDomain_measurable (a : Nat) : MeasurableSet (orderedSiteDomain a) := by
  unfold orderedSiteDomain StrictMono
  simp only [Set.ofPred_and, Set.ofPred_forall]
  apply MeasurableSet.inter
  · exact MeasurableSet.iInter fun i =>
      (measurableSet_le measurable_const (measurable_pi_apply i)).inter
        (measurableSet_lt (measurable_pi_apply i) measurable_const)
  · exact MeasurableSet.iInter fun i => MeasurableSet.iInter fun j =>
      MeasurableSet.iInter fun _ => measurableSet_lt (measurable_pi_apply i) (measurable_pi_apply j)

theorem orderedSiteDomain_bounded (a : Nat) : IsBounded (orderedSiteDomain a) := by
  apply (isCompact_Icc (a := (0 : Fin a → ℝ)) (b := 1)).isBounded.subset
  intro x hx
  exact ⟨fun i => (hx.1 i).1, fun i => (hx.1 i).2.le⟩

def integerGrid (a R : Nat) : Set (Fin a → ℝ) :=
  (R : ℝ)⁻¹ • (Submodule.span ℤ (Set.range (Pi.basisFun ℝ (Fin a))) : Set (Fin a → ℝ))

theorem orderedSite_riemann_tendsto {a : Nat} (F : (Fin a → ℝ) → ℝ) (hF : Continuous F) :
    Tendsto (fun R : Nat => (∑' x : ↑(orderedSiteDomain a ∩ integerGrid a R), F x) / (R : ℝ)^a)
      atTop (𝓝 (∫ x in orderedSiteDomain a, F x)) := by
  simpa only [Fintype.card_fin, integerGrid] using
    tendsto_tsum_div_pow_atTop_integral (orderedSiteDomain a) F hF (orderedSiteDomain_bounded a)
      (orderedSiteDomain_measurable a) ((orderedSiteDomain_convex a).addHaar_frontier volume)

end Spin.Structured.Placement



