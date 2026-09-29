import SpinCodes.Structured.ConcretePlacementLimitProduct

/-! Exact normalization and a quantitative finite-gap to simplex-integrand error. -/
noncomputable section
namespace Spin.Structured.Placement
open Finset ConcreteEncoder FiniteKernel

def simplexProduct (γ : ℝ) {a : Nat} (u : Fin (a + 1) → ℝ) : Matrix (Fin 2) (Fin 2) ℝ :=
  timeProduct γ (fun i => u i.castSucc) (u (Fin.last a))

theorem simplexProduct_rowError {a : Nat} {γ : ℝ} (hγ : 0 ≤ γ)
    (u v : Fin (a + 1) → ℝ) (hu : ∀ i, 0 ≤ u i) (hv : ∀ i, 0 ≤ v i) :
    RowError (simplexProduct γ u) (simplexProduct γ v) (γ * ∑ i, |u i - v i|) := by
  have h := timeProduct_rowError hγ (fun i => u i.castSucc) (fun i => v i.castSucc)
    (u (Fin.last a)) (v (Fin.last a)) (fun i => hu i.castSucc) (fun i => hv i.castSucc)
    (hu _) (hv _)
  simpa only [simplexProduct, Fin.sum_univ_castSucc, add_comm] using h

theorem coarseImpulseProduct_normalized {a : Nat} (c : ℝ) {L : ℝ} (hL : L ≠ 0)
    (gaps : Fin a → Nat) (last : Nat) :
    coarseImpulseProduct c gaps last =
      timeProduct (c * epochMean * L) (fun i => (gaps i : ℝ) / L) ((last : ℝ) / L) := by
  have he (g : Nat) : coarseEmpty c g = timeEmpty (c * epochMean * L) ((g : ℝ) / L) := by
    unfold coarseEmpty timeEmpty
    congr 2
    field_simp
  induction a with
  | zero => exact he last
  | succ a ih =>
    simp only [coarseImpulseProduct, timeProduct, he, ih]
    rfl

def placementProduct {R a : Nat} (θ : ℝ) (B : BlockSubset R a) : Matrix (Fin 2) (Fin 2) ℝ :=
  coarseImpulseProduct (θ / (128 * R)) (fun i => emptyGaps B i.castSucc) (emptyGaps B (Fin.last a))

theorem placementProduct_normalized {R a : Nat} (hR : 0 < R) (θ : ℝ) (B : BlockSubset R a) :
    placementProduct θ B = simplexProduct (θ * epochMean / 128)
      (fun i => (emptyGaps B i : ℝ) / R) := by
  have hR0 : (R : ℝ) ≠ 0 := by exact_mod_cast Nat.ne_of_gt hR
  unfold placementProduct simplexProduct
  rw [coarseImpulseProduct_normalized _ hR0]
  congr 1
  field_simp

theorem placementProduct_simplex_error {R a : Nat} (hR : 0 < R) {θ : ℝ} (hθ : 0 ≤ θ)
    (B : BlockSubset R a) :
    RowError (placementProduct θ B) (simplexProduct (θ * epochMean / 128) (simplexGaps B))
      ((θ * epochMean / 128) * ((a : ℝ) / R)) := by
  rw [placementProduct_normalized hR]
  have hγ : 0 ≤ θ * epochMean / 128 := by unfold epochMean; positivity
  have h := simplexProduct_rowError hγ (fun i => (emptyGaps B i : ℝ) / R) (simplexGaps B)
    (fun i => by positivity) (simplexGaps_nonneg B)
  have he : (∑ i, |(emptyGaps B i : ℝ) / R - simplexGaps B i|) = (a : ℝ) / R := by
    simp_rw [abs_sub_comm ((emptyGaps B _ : ℝ) / R)]
    exact simplexGaps_distance B
  rw [he] at h
  exact h

end Spin.Structured.Placement
