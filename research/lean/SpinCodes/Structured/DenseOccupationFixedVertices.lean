import SpinCodes.Structured.Convexity

noncomputable section
namespace Spin.Structured.DenseOccupationFixed
open Set

/-- The positive-occupation exponent for one affine outer segment and one fixed witness. -/
def boxExponent (m c p y radius z α x : ℝ) : ℝ :=
  α*(m*x+c) + binKL α p + α*binKL x y + Real.log radius/128 - (11/100)*Real.log z

theorem boxExponent_convex_alpha (m c p y radius z x : ℝ) (hp0 : 0<p) (hp1 : p<1) :
    ConvexOn ℝ (Icc (0:ℝ) 1) (fun α => boxExponent m c p y radius z α x) := by
  refine ⟨convex_Icc _ _, ?_⟩
  intro a ha b hb u v hu hv huv
  have h := (convexOn_binKL_left hp0 hp1).2 ha hb hu hv huv
  simp only [smul_eq_mul] at h ⊢
  unfold boxExponent
  nlinarith [h, congrArg (fun t => t*(Real.log radius/128-(11/100)*Real.log z)) huv]

theorem boxExponent_convex_x (m c p y radius z α : ℝ) (hα : 0≤α)
    (hy0 : 0<y) (hy1 : y<1) :
    ConvexOn ℝ (Icc (0:ℝ) 1) (fun x => boxExponent m c p y radius z α x) := by
  refine ⟨convex_Icc _ _, ?_⟩
  intro a ha b hb u v hu hv huv
  have h := (convexOn_binKL_left hy0 hy1).2 ha hb hu hv huv
  simp only [smul_eq_mul] at h ⊢
  have hh := mul_le_mul_of_nonneg_left h hα
  unfold boxExponent
  nlinarith [hh, congrArg
    (fun t => t*(α*c+binKL α p+Real.log radius/128-(11/100)*Real.log z)) huv]

/-- Four original (α,x) box corners suffice; no rectangular change-of-variables assumption. -/
theorem boxExponent_le_vertices (m c p y radius z : ℝ) (hp0 : 0<p) (hp1 : p<1)
    (hy0 : 0<y) (hy1 : y<1) {a₀ a₁ x₀ x₁ α x : ℝ}
    (ha₀ : a₀∈Icc (0:ℝ) 1) (ha₁ : a₁∈Icc (0:ℝ) 1)
    (hx₀ : x₀∈Icc (0:ℝ) 1) (hx₁ : x₁∈Icc (0:ℝ) 1)
    (hα : α∈Icc a₀ a₁) (hx : x∈Icc x₀ x₁) :
    boxExponent m c p y radius z α x ≤
      max (max (boxExponent m c p y radius z a₀ x₀) (boxExponent m c p y radius z a₀ x₁))
          (max (boxExponent m c p y radius z a₁ x₀) (boxExponent m c p y radius z a₁ x₁)) := by
  have h := (boxExponent_convex_alpha m c p y radius z x hp0 hp1).le_max_of_mem_Icc ha₀ ha₁ hα
  have h₀ := (boxExponent_convex_x m c p y radius z a₀ ha₀.1 hy0 hy1).le_max_of_mem_Icc hx₀ hx₁ hx
  have h₁ := (boxExponent_convex_x m c p y radius z a₁ ha₁.1 hy0 hy1).le_max_of_mem_Icc hx₀ hx₁ hx
  exact h.trans (max_le_max h₀ h₁)

/-- A certified vertex bound controls every point of its original rectangular box. -/
theorem boxExponent_le_of_vertices (m c p y radius z B : ℝ) (hp0 : 0<p) (hp1 : p<1)
    (hy0 : 0<y) (hy1 : y<1) {a₀ a₁ x₀ x₁ α x : ℝ}
    (ha₀ : a₀∈Icc (0:ℝ) 1) (ha₁ : a₁∈Icc (0:ℝ) 1)
    (hx₀ : x₀∈Icc (0:ℝ) 1) (hx₁ : x₁∈Icc (0:ℝ) 1)
    (hα : α∈Icc a₀ a₁) (hx : x∈Icc x₀ x₁)
    (h00 : boxExponent m c p y radius z a₀ x₀≤B)
    (h01 : boxExponent m c p y radius z a₀ x₁≤B)
    (h10 : boxExponent m c p y radius z a₁ x₀≤B)
    (h11 : boxExponent m c p y radius z a₁ x₁≤B) :
    boxExponent m c p y radius z α x≤B :=
  (boxExponent_le_vertices m c p y radius z hp0 hp1 hy0 hy1 ha₀ ha₁ hx₀ hx₁ hα hx).trans
    (max_le (max_le h00 h01) (max_le h10 h11))

#print axioms boxExponent_le_of_vertices
end Spin.Structured.DenseOccupationFixed
