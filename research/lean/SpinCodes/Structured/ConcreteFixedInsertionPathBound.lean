import SpinCodes.Structured.ConcreteFixedInsertionLivePaths

/-! Scalar live-duration bounds control the full nonnegative matrix path integral. -/
noncomputable section
namespace Spin.Structured.Placement
open Finset ConcreteEncoder Matrix MeasureTheory
attribute [local instance] Classical.propDecidable

def twoWeight (v : ℝ) : Fin 2 → ℝ := ![1,v]

theorem rowNormLe_action {v c : ℝ} (hv : 0 ≤ v) {K : Matrix (Fin 2) (Fin 2) ℝ}
    (hK : Spin.RowNormLe v c K) (i : Fin 2) : (K *ᵥ twoWeight v) i ≤ c * twoWeight v i := by
  fin_cases i
  · simp [Matrix.mulVec, dotProduct, Fin.sum_univ_two, twoWeight]
    exact (add_le_add (le_abs_self _) (mul_le_mul_of_nonneg_right (le_abs_self _) hv)).trans hK.1
  · simp [Matrix.mulVec, dotProduct, Fin.sum_univ_two, twoWeight]
    exact (add_le_add (le_abs_self _) (mul_le_mul_of_nonneg_right (le_abs_self _) hv)).trans hK.2

theorem mv_product_nonneg {Q : Nat} {σ v r : ℝ} (hσ : 0 ≤ σ) (hv : 0 < v)
    (u : Fin Q → ℝ) (hu : ∀ k, 0 ≤ u k) : 0 ≤ ∏ k, Spin.mv σ v r (u k) := by
  apply prod_nonneg
  intro k hk
  have huk := hu k
  exact (show 0 ≤ 1+u k*σ*v by positivity).trans (le_max_left _ _)

theorem constant_path_action_bound {Q : Nat} {σ v r : ℝ} (hσ0 : 0 ≤ σ) (hσ1 : σ ≤ 1)
    (hv : 0 < v) (hr : 0 ≤ r) (u : Fin Q → ℝ) (hu : ∀ k, 0 ≤ u k) (i : Fin 2) :
    (diagonalPathProduct (fun _ : Fin (Q+1) => constantDiagonal σ) (fun k => 1+u k • Spin.Pplus r) *ᵥ twoWeight v) i ≤
      (∏ k, Spin.mv σ v r (u k)) * twoWeight v i := by
  rw [constant_pathProduct]
  have hP := Spin.rowNormLe_pathProduct hσ0 hv hr (List.ofFn u) (by
    intro t ht; obtain ⟨k,rfl⟩ := List.mem_ofFn.mp ht; exact hu k)
  have hP' : Spin.RowNormLe v (∏ k, Spin.mv σ v r (u k)) (Spin.pathProduct σ r (List.ofFn u)) := by
    simpa only [List.map_ofFn, List.prod_ofFn, Function.comp_def] using hP
  have h := (Spin.rowNormLe_Dsig hσ0 hσ1 hv.le).mul hv.le (mv_product_nonneg hσ0 hv u hu) hP'
  simpa only [one_mul] using rowNormLe_action hv.le h i

theorem pathWeight_nonneg {Q : Nat} {v r : ℝ} (hv : 0 ≤ v) (hr : 0 ≤ r)
    (u : Fin Q → ℝ) (hu : ∀ k, 0 ≤ u k) (i : Fin 2) (ys : Fin Q → Fin 2) :
    0 ≤ pathWeight (fun k => 1+u k • Spin.Pplus r) (twoWeight v) i ys := by
  have hP (i j : Fin 2) : 0 ≤ Spin.Pplus r i j := by fin_cases i <;> fin_cases j <;> simp [Spin.Pplus,hr]
  have hw (i : Fin 2) : 0 ≤ twoWeight v i := by fin_cases i <;> simp [twoWeight,hv]
  apply mul_nonneg _ (hw _)
  apply prod_nonneg
  intro k hk
  exact add_nonneg (FiniteKernel.Substochastic.one.nonneg _ _) (mul_nonneg (hu k) (hP _ _))

theorem continuous_integrable_ordered {Q : Nat} (F : (Fin Q → ℝ) → ℝ) (hF : Continuous F) :
    IntegrableOn F (orderedSiteDomain Q) := by
  apply hF.integrableOn_Icc.mono_set
  intro x hx
  exact ⟨fun i => (hx.1 i).1, fun i => (hx.1 i).2.le⟩

theorem fullSiteProduct_integrated_path_bound {Q : Nat} {γ σ v r C : ℝ}
    (hσ0 : 0 ≤ σ) (hσ1 : σ ≤ 1) (hv : 0 < v) (hr : 0 ≤ r) (hC : 0 ≤ C)
    (hscalar : ∀ S : Finset (Fin (Q+1)), (Q.factorial : ℝ) *
      (∫ x in orderedSiteDomain Q, Real.exp (-γ * liveDuration S x)) ≤ C * σ^S.card)
    (u : Fin Q → ℝ) (hu : ∀ k, 0 ≤ u k) (i : Fin 2) :
    (Q.factorial : ℝ) * (∫ x in orderedSiteDomain Q,
      (fullSiteProduct γ (Spin.Pplus r) x u *ᵥ twoWeight v) i) ≤
        C * (∏ k, Spin.mv σ v r (u k)) * twoWeight v i := by
  simp_rw [fullSiteProduct_path_expansion]
  have hint (ys : Fin Q → Fin 2) : IntegrableOn
      (fun x : Fin Q → ℝ => Real.exp (-γ * liveDuration (pathLive i ys) x) *
        pathWeight (fun k => 1+u k • Spin.Pplus r) (twoWeight v) i ys) (orderedSiteDomain Q) := by
    apply continuous_integrable_ordered
    exact (Real.continuous_exp.comp (continuous_const.mul (continuous_liveDuration (pathLive i ys)))).mul continuous_const
  rw [integral_finsetSum univ (fun ys hys => hint ys), Finset.mul_sum]
  calc
    _ = ∑ ys : Fin Q → Fin 2, ((Q.factorial : ℝ) *
        (∫ x in orderedSiteDomain Q, Real.exp (-γ * liveDuration (pathLive i ys) x))) *
          pathWeight (fun k => 1+u k • Spin.Pplus r) (twoWeight v) i ys := by
      apply sum_congr rfl
      intro ys hys
      rw [integral_mul_const]
      ring
    _ ≤ ∑ ys : Fin Q → Fin 2, (C * σ^(pathLive i ys).card) *
        pathWeight (fun k => 1+u k • Spin.Pplus r) (twoWeight v) i ys := by
      apply sum_le_sum
      intro ys hys
      exact mul_le_mul_of_nonneg_right (hscalar _) (pathWeight_nonneg hv.le hr u hu i ys)
    _ = C * (diagonalPathProduct (fun _ : Fin (Q+1) => constantDiagonal σ)
        (fun k => 1+u k • Spin.Pplus r) *ᵥ twoWeight v) i := by
      rw [constant_path_expansion, Finset.mul_sum]
      apply sum_congr rfl
      intro ys hys
      ring
    _ ≤ C * ((∏ k, Spin.mv σ v r (u k)) * twoWeight v i) :=
      mul_le_mul_of_nonneg_left (constant_path_action_bound hσ0 hσ1 hv hr u hu i) hC
    _ = _ := by ring

end Spin.Structured.Placement


