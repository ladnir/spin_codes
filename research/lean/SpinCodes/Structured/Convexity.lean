/-
Convexity of the positive-occupation exponent, and the vertex check.

    "The verifier covers [10⁻⁴,1]×𝒲 by 1,023 rational boxes within the 39
     outer segments.  It fixes a witness on each box and checks four vertices.
     The exponent is convex in (α, αx) on each affine segment."

Two things are needed to justify checking only the corners.

*  **The change of variables.**  In `(u, w) = (α, αx)` the exponent's three
   non-constant terms are `u·â_BA(w/u)`, `D(u‖p)` and `u·D(w/u‖y)`.  The first
   is *linear* exactly when `â_BA` is affine — which is what "on each affine
   segment" buys.  The third is the perspective transform of a convex
   function, which is convex; Mathlib does not have it, so it is proved here.

*  **The vertex principle.**  A convex function on a rectangle is at most the
   largest of its four corner values.  This is the 1-D maximum principle
   applied along each axis in turn.
-/
import SpinCodes.Structured.KLChain

namespace Spin

open Set

private lemma pos_comb {u₁ u₂ a b : ℝ} (hu₁ : 0 < u₁) (hu₂ : 0 < u₂)
    (ha : 0 ≤ a) (hb : 0 ≤ b) (hab : a + b = 1) : 0 < a * u₁ + b * u₂ := by
  rcases eq_or_lt_of_le ha with ha0 | ha0
  · have hb1 : b = 1 := by linarith
    rw [← ha0, hb1]
    simpa using hu₂
  · nlinarith

/-! ## Restricting a convex function to an axis-parallel line -/

lemma convexOn_slice₁ {f : ℝ × ℝ → ℝ} {s₁ s₂ : Set ℝ} (hs₁ : Convex ℝ s₁)
    (hf : ConvexOn ℝ (s₁ ×ˢ s₂) f) {c : ℝ} (hc : c ∈ s₂) :
    ConvexOn ℝ s₁ (fun t => f (t, c)) := by
  refine ⟨hs₁, fun x hx y hy a b ha hb hab => ?_⟩
  have h := hf.2 (Set.mk_mem_prod hx hc) (Set.mk_mem_prod hy hc) ha hb hab
  have hcc : a * c + b * c = c := by
    have : (a + b) * c = c := by rw [hab, one_mul]
    linarith [this]
  simpa [Prod.smul_mk, Prod.mk_add_mk, smul_eq_mul, hcc] using h

lemma convexOn_slice₂ {f : ℝ × ℝ → ℝ} {s₁ s₂ : Set ℝ} (hs₂ : Convex ℝ s₂)
    (hf : ConvexOn ℝ (s₁ ×ˢ s₂) f) {c : ℝ} (hc : c ∈ s₁) :
    ConvexOn ℝ s₂ (fun t => f (c, t)) := by
  refine ⟨hs₂, fun x hx y hy a b ha hb hab => ?_⟩
  have h := hf.2 (Set.mk_mem_prod hc hx) (Set.mk_mem_prod hc hy) ha hb hab
  have hcc : a * c + b * c = c := by
    have : (a + b) * c = c := by rw [hab, one_mul]
    linarith [this]
  simpa [Prod.smul_mk, Prod.mk_add_mk, smul_eq_mul, hcc] using h

/-! ## The vertex principle -/

/-- **Checking four vertices suffices.**  A convex function on a rectangle is
at most the largest of its corner values. -/
theorem le_max_vertices {f : ℝ × ℝ → ℝ} {a₁ b₁ a₂ b₂ : ℝ}
    (hf : ConvexOn ℝ (Icc a₁ b₁ ×ˢ Icc a₂ b₂) f)
    {z : ℝ × ℝ} (hz : z ∈ Icc a₁ b₁ ×ˢ Icc a₂ b₂) :
    f z ≤ max (max (f (a₁, a₂)) (f (a₁, b₂))) (max (f (b₁, a₂)) (f (b₁, b₂))) := by
  obtain ⟨hz1, hz2⟩ := Set.mem_prod.mp hz
  have h1 : a₁ ≤ b₁ := le_trans hz1.1 hz1.2
  have h2 : a₂ ≤ b₂ := le_trans hz2.1 hz2.2
  have hslice : ConvexOn ℝ (Icc a₁ b₁) (fun t => f (t, z.2)) :=
    convexOn_slice₁ (convex_Icc _ _) hf hz2
  have step1 : f z ≤ max (f (a₁, z.2)) (f (b₁, z.2)) := by
    have h := hslice.le_max_of_mem_Icc (left_mem_Icc.mpr h1) (right_mem_Icc.mpr h1) hz1
    simpa using h
  have hA : ConvexOn ℝ (Icc a₂ b₂) (fun t => f (a₁, t)) :=
    convexOn_slice₂ (convex_Icc _ _) hf (left_mem_Icc.mpr h1)
  have hB : ConvexOn ℝ (Icc a₂ b₂) (fun t => f (b₁, t)) :=
    convexOn_slice₂ (convex_Icc _ _) hf (right_mem_Icc.mpr h1)
  have stepA : f (a₁, z.2) ≤ max (f (a₁, a₂)) (f (a₁, b₂)) :=
    hA.le_max_of_mem_Icc (left_mem_Icc.mpr h2) (right_mem_Icc.mpr h2) hz2
  have stepB : f (b₁, z.2) ≤ max (f (b₁, a₂)) (f (b₁, b₂)) :=
    hB.le_max_of_mem_Icc (left_mem_Icc.mpr h2) (right_mem_Icc.mpr h2) hz2
  exact step1.trans (max_le_max stepA stepB)

/-! ## The perspective transform

`(u, w) ↦ u · g(w/u)` is convex on the cone over `g`'s domain.  This is what
makes `α · D(x‖y)` convex in `(α, αx)`. -/

/-- The cone over `t`: the pairs `(u,w)` with `u > 0` and `w/u ∈ t`. -/
def perspectiveDomain (t : Set ℝ) : Set (ℝ × ℝ) :=
  {q : ℝ × ℝ | 0 < q.1 ∧ q.2 / q.1 ∈ t}

/-- The ratio at a convex combination is a convex combination of the ratios,
with weights proportional to the `u`-components. -/
private lemma ratio_comb {u₁ w₁ u₂ w₂ a b : ℝ} (hu : 0 < a * u₁ + b * u₂)
    (hu₁ : 0 < u₁) (hu₂ : 0 < u₂) :
    (a * w₁ + b * w₂) / (a * u₁ + b * u₂)
      = (a * u₁ / (a * u₁ + b * u₂)) * (w₁ / u₁)
        + (b * u₂ / (a * u₁ + b * u₂)) * (w₂ / u₂) := by
  field_simp

private lemma weights_sum {u₁ u₂ a b : ℝ} (hu : 0 < a * u₁ + b * u₂) :
    a * u₁ / (a * u₁ + b * u₂) + b * u₂ / (a * u₁ + b * u₂) = 1 := by
  field_simp

lemma convex_perspectiveDomain {t : Set ℝ} (ht : Convex ℝ t) :
    Convex ℝ (perspectiveDomain t) := by
  rintro ⟨u₁, w₁⟩ ⟨hu₁, hq₁⟩ ⟨u₂, w₂⟩ ⟨hu₂, hq₂⟩ a b ha hb hab
  have hu : 0 < a * u₁ + b * u₂ := pos_comb hu₁ hu₂ ha hb hab
  have hw1 : 0 ≤ a * u₁ / (a * u₁ + b * u₂) := by positivity
  have hw2 : 0 ≤ b * u₂ / (a * u₁ + b * u₂) := by positivity
  have hmem := ht hq₁ hq₂ hw1 hw2 (weights_sum hu)
  simp only [smul_eq_mul] at hmem
  refine ⟨by simpa using hu, ?_⟩
  simpa [ratio_comb hu hu₁ hu₂] using hmem

/-- **The perspective transform of a convex function is convex.** -/
theorem convexOn_perspective {g : ℝ → ℝ} {t : Set ℝ} (hg : ConvexOn ℝ t g) :
    ConvexOn ℝ (perspectiveDomain t) (fun q => q.1 * g (q.2 / q.1)) := by
  refine ⟨convex_perspectiveDomain hg.1, ?_⟩
  rintro ⟨u₁, w₁⟩ ⟨hu₁, hq₁⟩ ⟨u₂, w₂⟩ ⟨hu₂, hq₂⟩ a b ha hb hab
  have hu : 0 < a * u₁ + b * u₂ := pos_comb hu₁ hu₂ ha hb hab
  have hw1 : 0 ≤ a * u₁ / (a * u₁ + b * u₂) := by positivity
  have hw2 : 0 ≤ b * u₂ / (a * u₁ + b * u₂) := by positivity
  have hjensen := hg.2 hq₁ hq₂ hw1 hw2 (weights_sum hu)
  simp only [smul_eq_mul] at hjensen
  have hmul := mul_le_mul_of_nonneg_left hjensen (le_of_lt hu)
  have hR : (a * u₁ + b * u₂)
        * ((a * u₁ / (a * u₁ + b * u₂)) * g (w₁ / u₁)
            + (b * u₂ / (a * u₁ + b * u₂)) * g (w₂ / u₂))
      = a * (u₁ * g (w₁ / u₁)) + b * (u₂ * g (w₂ / u₂)) := by
    field_simp
  rw [hR] at hmul
  simpa [Prod.smul_mk, Prod.mk_add_mk, smul_eq_mul, ← ratio_comb hu hu₁ hu₂] using hmul

/-! ## Convexity of `D(·‖p)` -/

private lemma convexOn_affine (s : Set ℝ) (hs : Convex ℝ s) (c d : ℝ) :
    ConvexOn ℝ s (fun u => c * u + d) := by
  refine ⟨hs, fun x _ y _ a b ha hb hab => ?_⟩
  simp only [smul_eq_mul]
  have : c * (a * x + b * y) + d = a * (c * x + d) + b * (c * y + d) := by
    linear_combination (-d) * hab
  exact le_of_eq this

private lemma convexOn_one_sub_mul_log :
    ConvexOn ℝ (Iic (1 : ℝ)) (fun u => (1 - u) * Real.log (1 - u)) := by
  refine ⟨convex_Iic 1, fun x hx y hy a b ha hb hab => ?_⟩
  have h1 : (1 : ℝ) - x ∈ Ici (0 : ℝ) := by simpa using hx
  have h2 : (1 : ℝ) - y ∈ Ici (0 : ℝ) := by simpa using hy
  have h := Real.convexOn_mul_log.2 h1 h2 ha hb hab
  simp only [smul_eq_mul] at h ⊢
  have hcomb : 1 - (a * x + b * y) = a * (1 - x) + b * (1 - y) := by
    linear_combination (-1 : ℝ) * hab
  rw [hcomb]
  exact h

/-- **`D(·‖p)` is convex.** -/
theorem convexOn_binKL_left {p : ℝ} (hp0 : 0 < p) (hp1 : p < 1) :
    ConvexOn ℝ (Icc (0 : ℝ) 1) (fun u => binKL u p) := by
  have hsplit : ∀ u ∈ Icc (0 : ℝ) 1,
      binKL u p = (u * Real.log u) + ((1 - u) * Real.log (1 - u))
        + ((- Real.log p + Real.log (1 - p)) * u + (- Real.log (1 - p))) := by
    intro u hu
    have hu0 : 0 ≤ u := hu.1
    have hu1 : u ≤ 1 := hu.2
    have e1 : u * Real.log (u / p) = u * Real.log u - u * Real.log p := by
      rcases eq_or_lt_of_le hu0 with h | h
      · simp [← h]
      · rw [Real.log_div (ne_of_gt h) (ne_of_gt hp0)]; ring
    have e2 : (1 - u) * Real.log ((1 - u) / (1 - p))
        = (1 - u) * Real.log (1 - u) - (1 - u) * Real.log (1 - p) := by
      rcases eq_or_lt_of_le hu1 with h | h
      · simp [h]
      · have hpos : (0 : ℝ) < 1 - u := by linarith
        rw [Real.log_div (ne_of_gt hpos) (by linarith : (1 : ℝ) - p ≠ 0)]; ring
    unfold binKL
    rw [e1, e2]
    ring
  have hconv : ConvexOn ℝ (Icc (0 : ℝ) 1)
      (fun u => (u * Real.log u) + ((1 - u) * Real.log (1 - u))
        + ((- Real.log p + Real.log (1 - p)) * u + (- Real.log (1 - p)))) := by
    refine ConvexOn.add (ConvexOn.add ?_ ?_) (convexOn_affine _ (convex_Icc _ _) _ _)
    · exact Real.convexOn_mul_log.subset (fun u hu => hu.1) (convex_Icc _ _)
    · exact convexOn_one_sub_mul_log.subset (fun u hu => hu.2) (convex_Icc _ _)
  refine ⟨convex_Icc _ _, fun x hx y hy a b ha hb hab => ?_⟩
  have hmem : a • x + b • y ∈ Icc (0 : ℝ) 1 := (convex_Icc (0 : ℝ) 1) hx hy ha hb hab
  show binKL (a • x + b • y) p ≤ a • binKL x p + b • binKL y p
  rw [hsplit _ hmem, hsplit x hx, hsplit y hy]
  exact hconv.2 hx hy ha hb hab

end Spin
