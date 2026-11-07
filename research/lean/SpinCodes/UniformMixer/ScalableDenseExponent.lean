import SpinCodes.UniformMixer.ScalableDenseTilted
import SpinCodes.Structured.DenseOccupationFixedVertices

/-! The dense Scalable exponent has one 64-bit transfer step per 64 output bits. -/

noncomputable section

namespace Spin.UniformMixer

open Spin.Structured Spin.Structured.DenseOccupationFixed Set

def scalableBoxExponent (m c p y ρ z α x : ℝ) : ℝ :=
  α * (m * x + c) + binKL α p + α * binKL x y +
    Real.log ρ / 64 - (11 / 100) * Real.log z

theorem scalableBoxExponent_convex_alpha (m c p y ρ z x : ℝ)
    (hp0 : 0 < p) (hp1 : p < 1) :
    ConvexOn ℝ (Icc (0 : ℝ) 1)
      (fun α => scalableBoxExponent m c p y ρ z α x) := by
  refine ⟨convex_Icc _ _, ?_⟩
  intro a ha b hb u v hu hv huv
  have h := (convexOn_binKL_left hp0 hp1).2 ha hb hu hv huv
  simp only [smul_eq_mul] at h ⊢
  unfold scalableBoxExponent
  nlinarith [h, congrArg
    (fun t => t * (Real.log ρ / 64 - (11 / 100) * Real.log z)) huv]

theorem scalableBoxExponent_convex_x (m c p y ρ z α : ℝ)
    (hα : 0 ≤ α) (hy0 : 0 < y) (hy1 : y < 1) :
    ConvexOn ℝ (Icc (0 : ℝ) 1)
      (fun x => scalableBoxExponent m c p y ρ z α x) := by
  refine ⟨convex_Icc _ _, ?_⟩
  intro a ha b hb u v hu hv huv
  have h := (convexOn_binKL_left hy0 hy1).2 ha hb hu hv huv
  simp only [smul_eq_mul] at h ⊢
  have hh := mul_le_mul_of_nonneg_left h hα
  unfold scalableBoxExponent
  nlinarith [hh, congrArg
    (fun t => t * (α * c + binKL α p +
      Real.log ρ / 64 - (11 / 100) * Real.log z)) huv]

theorem scalableBoxExponent_le_vertices (m c p y ρ z : ℝ)
    (hp0 : 0 < p) (hp1 : p < 1)
    (hy0 : 0 < y) (hy1 : y < 1)
    {a₀ a₁ x₀ x₁ α x : ℝ}
    (ha₀ : a₀ ∈ Icc (0 : ℝ) 1) (ha₁ : a₁ ∈ Icc (0 : ℝ) 1)
    (hx₀ : x₀ ∈ Icc (0 : ℝ) 1) (hx₁ : x₁ ∈ Icc (0 : ℝ) 1)
    (hα : α ∈ Icc a₀ a₁) (hx : x ∈ Icc x₀ x₁) :
    scalableBoxExponent m c p y ρ z α x ≤
      max (max (scalableBoxExponent m c p y ρ z a₀ x₀)
        (scalableBoxExponent m c p y ρ z a₀ x₁))
        (max (scalableBoxExponent m c p y ρ z a₁ x₀)
          (scalableBoxExponent m c p y ρ z a₁ x₁)) := by
  have h := (scalableBoxExponent_convex_alpha m c p y ρ z x hp0 hp1).le_max_of_mem_Icc
    ha₀ ha₁ hα
  have h₀ := (scalableBoxExponent_convex_x m c p y ρ z a₀ ha₀.1 hy0 hy1).le_max_of_mem_Icc
    hx₀ hx₁ hx
  have h₁ := (scalableBoxExponent_convex_x m c p y ρ z a₁ ha₁.1 hy0 hy1).le_max_of_mem_Icc
    hx₀ hx₁ hx
  exact h.trans (max_le_max h₀ h₁)

theorem scalable_dense_finite_exponent_le
    {N α x p y ρ z m c : ℝ} {R d : ℕ}
    (hN : 0 ≤ N) (hα0 : 0 < α) (hα1 : α ≤ 1)
    (hx0 : 0 < x) (hx1 : x ≤ 1)
    (hp0 : 0 < p) (hp1 : p < 1)
    (hy0 : 0 < y) (hy1 : y < 1)
    (hρ : 0 < ρ) (hz0 : 0 < z) (hz1 : z ≤ 1)
    (hround : 64 * (R : ℝ) = N)
    (hd : (d : ℝ) ≤ (11 / 100) * N) :
    Real.exp (N * α * (m * x + c)) *
      Real.exp (N * binKL (α * x) (p * y)) * ρ ^ R / z ^ d ≤
      Real.exp (N * scalableBoxExponent m c p y ρ z α x) := by
  have heR : ρ ^ R = Real.exp ((R : ℝ) * Real.log ρ) := by
    rw [Real.exp_nat_mul, Real.exp_log hρ]
  have hez : z ^ d = Real.exp ((d : ℝ) * Real.log z) := by
    rw [Real.exp_nat_mul, Real.exp_log hz0]
  rw [heR, hez, ← Real.exp_add, ← Real.exp_add, ← Real.exp_sub]
  apply Real.exp_le_exp.mpr
  have hk := mul_le_mul_of_nonneg_left
    (binKL_chain hα0 hα1 hx0 hx1 hp0 hp1 hy0 hy1) hN
  have hl : Real.log z ≤ 0 := Real.log_nonpos hz0.le hz1
  have hdt := mul_le_mul_of_nonneg_right hd (neg_nonneg.mpr hl)
  unfold scalableBoxExponent
  nlinarith [congrArg (fun t => t * Real.log ρ / 64) hround]

#print axioms scalableBoxExponent_le_vertices
#print axioms scalable_dense_finite_exponent_le

end Spin.UniformMixer
