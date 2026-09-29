import SpinCodes.Structured.ConcreteScalarBits
import SpinCodes.Structured.ConcretePairing

noncomputable section
attribute [local instance] Classical.propDecidable
namespace Spin.Structured.ConcreteFourier
open Finset

theorem poissonBinom_product_moment {n : ℕ} (p : Fin n → ℝ)
    (hp0 : ∀ i, 0 ≤ p i) (hp1 : ∀ i, p i ≤ 1) (z : Fin n → ℝ) :
    (poissonBinom p hp0 hp1).expect (fun X => ∏ i ∈ X, z i) =
      ∏ i, (1-p i+p i*z i) := by
  have h := Finset.prod_add (fun i => p i*z i) (fun i => 1-p i) (univ : Finset (Fin n))
  simp only [powerset_univ] at h
  unfold FinPMF.expect
  calc
    _ = ∑ X : Finset (Fin n), (∏ i ∈ X, p i*z i) * ∏ i ∈ univ \ X, (1-p i) := by
      apply sum_congr rfl
      intro X _
      rw [poissonBinom_apply, prod_mul_distrib]
      ring
    _ = ∏ i, (p i*z i+(1-p i)) := h.symm
    _ = _ := by apply prod_congr rfl; intro i _; ring

theorem chi_prod {n : ℕ} (v X : Finset (Fin n)) :
    (Spin.chi v X : ℝ) = ∏ i ∈ X, if i ∈ v then (-1 : ℝ) else 1 := by
  rw [Spin.chi]
  push_cast
  rw [Finset.prod_ite]
  simp only [prod_const, one_pow, mul_one]
  congr 1
  congr 1
  ext i
  simp only [mem_inter, mem_filter]
  tauto

/-- Fourier transform of the independent coordinate law, including endpoint probabilities. -/
theorem poissonBinom_character {n : ℕ} (p : Fin n → ℝ)
    (hp0 : ∀ i, 0 ≤ p i) (hp1 : ∀ i, p i ≤ 1) (v : Finset (Fin n)) :
    (poissonBinom p hp0 hp1).expect (fun X => (Spin.chi v X : ℝ)) =
      ∏ i ∈ v, (1-2*p i) := by
  simp_rw [chi_prod]
  rw [poissonBinom_product_moment]
  have h : (∏ i : Fin n, (1-p i+p i*(if i ∈ v then (-1 : ℝ) else 1))) =
      ∏ i : Fin n, if i ∈ v then 1-2*p i else 1 := by
    apply prod_congr rfl
    intro i _
    by_cases hi : i ∈ v <;> simp only [hi, ite_true, ite_false] <;> ring
  rw [h, prod_ite]
  have hv : (univ : Finset (Fin n)).filter (fun i => i ∈ v) = v := by ext; simp
  rw [hv]
  simp

end Spin.Structured.ConcreteFourier

