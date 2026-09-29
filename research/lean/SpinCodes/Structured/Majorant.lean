/-
The concave outer majorant and the Jensen step (`eq:structured-ba-jensen`).

The paper's certificate does not hand us an abstract concave function: it
represents `â_BA` concretely, as "the lower envelope of rational affine
supports" (39 segments in the refined version), and notes that "taking the
lower envelope preserves concavity".

That is the whole content of the Jensen step, and it is unconditional — it
depends on the *representation*, not on the numerical values of the supports.
So nothing here is certificate-dependent: whatever the 39 rational pairs turn
out to be, concavity and the averaging inequality hold. What the certificate
still owes is the majorization `a_BA ≤ â_BA` on `𝒲` (`eq:ba-spectrum-majorant`,
queued as T7), which is a different statement.
-/
import Mathlib

set_option linter.unusedSectionVars false

namespace Spin.Structured

open Finset

/-- An affine function is concave on any convex set. -/
lemma concaveOn_affine (s : Set ℝ) (hs : Convex ℝ s) (a c : ℝ) :
    ConcaveOn ℝ s (fun x => a * x + c) := by
  refine ⟨hs, fun x _ y _ p q hp hq hpq => ?_⟩
  simp only [smul_eq_mul]
  linear_combination c * hpq

/-- The pointwise `inf'` of finitely many concave functions is concave.
This is the "taking the lower envelope preserves concavity" step. -/
lemma concaveOn_finset_inf' {ι : Type*} (s : Set ℝ) (hs : Convex ℝ s)
    (P : Finset ι) (hP : P.Nonempty) (f : ι → ℝ → ℝ)
    (hf : ∀ i ∈ P, ConcaveOn ℝ s (f i)) :
    ConcaveOn ℝ s (fun x => P.inf' hP (fun i => f i x)) := by
  refine ⟨hs, fun x hx y hy p q hp hq hpq => ?_⟩
  simp only [smul_eq_mul]
  refine Finset.le_inf' hP _ fun i hi => ?_
  have hix : P.inf' hP (fun i => f i x) ≤ f i x := Finset.inf'_le _ hi
  have hiy : P.inf' hP (fun i => f i y) ≤ f i y := Finset.inf'_le _ hi
  have hconc := (hf i hi).2 hx hy hp hq hpq
  simp only [smul_eq_mul] at hconc
  nlinarith [mul_le_mul_of_nonneg_left hix hp, mul_le_mul_of_nonneg_left hiy hq]

/-- **Jensen, in the averaging form the paper uses.**  For a concave `f`,

  `(1/Q) Σ_j f(x_j) ≤ f( (1/Q) Σ_j x_j )`. -/
theorem mean_le_of_concaveOn {ι : Type*} (s : Set ℝ) (f : ℝ → ℝ)
    (hf : ConcaveOn ℝ s f) (t : Finset ι) (ht : t.Nonempty) (x : ι → ℝ)
    (hx : ∀ i ∈ t, x i ∈ s) :
    (∑ i ∈ t, f (x i)) / t.card ≤ f ((∑ i ∈ t, x i) / t.card) := by
  have hcard : (0 : ℝ) < t.card := by
    exact_mod_cast Finset.card_pos.mpr ht
  have hw0 : ∀ i ∈ t, (0 : ℝ) ≤ 1 / t.card := fun i _ => by positivity
  have hw1 : ∑ _i ∈ t, (1 / (t.card : ℝ)) = 1 := by
    rw [Finset.sum_const, nsmul_eq_mul]
    field_simp
  have h := hf.le_map_sum hw0 hw1 hx
  simp only [smul_eq_mul] at h
  calc (∑ i ∈ t, f (x i)) / t.card = ∑ i ∈ t, 1 / (t.card : ℝ) * f (x i) := by
        rw [Finset.sum_div]
        exact Finset.sum_congr rfl fun i _ => by ring
    _ ≤ f (∑ i ∈ t, 1 / (t.card : ℝ) * x i) := h
    _ = f ((∑ i ∈ t, x i) / t.card) := by
        congr 1
        rw [Finset.sum_div]
        exact Finset.sum_congr rfl fun i _ => by ring

/-- The outer majorant `â_BA`, represented exactly as the certificate
represents it: the lower envelope of a finite set of rational affine supports
`(slope, intercept)`. -/
structure ConcaveMajorant where
  supports : Finset (ℚ × ℚ)
  nonempty : supports.Nonempty

namespace ConcaveMajorant

variable (M : ConcaveMajorant)

/-- `â_BA(x) = min over supports of (slope · x + intercept)`. -/
noncomputable def toFun (x : ℝ) : ℝ :=
  M.supports.inf' M.nonempty (fun p => (p.1 : ℝ) * x + (p.2 : ℝ))

theorem concaveOn (s : Set ℝ) (hs : Convex ℝ s) : ConcaveOn ℝ s M.toFun :=
  concaveOn_finset_inf' s hs M.supports M.nonempty _
    (fun p _ => concaveOn_affine s hs _ _)

/-- `â_BA` never exceeds any of its affine supports. -/
theorem toFun_le_support {p : ℚ × ℚ} (hp : p ∈ M.supports) (x : ℝ) :
    M.toFun x ≤ (p.1 : ℝ) * x + (p.2 : ℝ) :=
  Finset.inf'_le _ hp

/-- **`eq:structured-ba-jensen`.**  Concavity of `â_BA` reduces the selected
outer spectrum to the mean row weight:

  `(1/Q) Σ_j â_BA(x_j) ≤ â_BA(x)`,  `x = (1/Q) Σ_j x_j`. -/
theorem jensen {ι : Type*} (s : Set ℝ) (hs : Convex ℝ s) (t : Finset ι)
    (ht : t.Nonempty) (x : ι → ℝ) (hx : ∀ i ∈ t, x i ∈ s) :
    (∑ i ∈ t, M.toFun (x i)) / t.card ≤ M.toFun ((∑ i ∈ t, x i) / t.card) :=
  mean_le_of_concaveOn s M.toFun (M.concaveOn s hs) t ht x hx

/-- The certified window `𝒲 = [0.104, 0.896]` is convex, so the Jensen step
applies on it. -/
theorem jensen_on_window {ι : Type*} (t : Finset ι) (ht : t.Nonempty)
    (x : ι → ℝ) (hx : ∀ i ∈ t, x i ∈ Set.Icc (0.104 : ℝ) 0.896) :
    (∑ i ∈ t, M.toFun (x i)) / t.card ≤ M.toFun ((∑ i ∈ t, x i) / t.card) :=
  M.jensen _ (convex_Icc _ _) t ht x hx

end ConcaveMajorant

end Spin.Structured
