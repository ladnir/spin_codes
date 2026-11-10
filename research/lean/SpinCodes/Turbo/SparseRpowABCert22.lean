import SpinCodes.Turbo.SparseRpowABDefs22
import SpinCodes.Turbo.SparseRpowCert22

/-!
# Integer certificates for entrywise `x ^ (a/b)` upper bounds

Generalizes `SparseRpowCert22` (α = 4/5) to any rational exponent
`α = a/b` with `b ≠ 0`. The kernel check compares `X^a · sy^b ≤ Y^b · sx^a`
entrywise by lockstep list traversal. Soundness turns this into
`(M i j) ^ (a/b) ≤ Y_ij / sy` for any nonnegative real matrix `M`
dominated by `X / sx`.
-/
noncomputable section
namespace Spin.Turbo
namespace SparseRpowABCert22
open SparseListSound22

/-- `x^a ≤ y^b` gives `x^(a/b) ≤ y` for nonnegative reals. -/
theorem rpowAB_le_of_pow (a b : ℕ) (hb : b ≠ 0) {x y : ℝ} (hx : 0 ≤ x) (hy : 0 ≤ y)
    (h : x ^ a ≤ y ^ b) : x ^ ((a : ℝ) / b) ≤ y := by
  have hbR : (b : ℝ) ≠ 0 := by exact_mod_cast hb
  have h1 : (x ^ ((a : ℝ) / b)) ^ b = x ^ a := by
    rw [← Real.rpow_natCast (x ^ ((a : ℝ) / b)) b, ← Real.rpow_mul hx,
      div_mul_cancel₀ _ hbR, Real.rpow_natCast]
  have h0 : 0 ≤ x ^ ((a : ℝ) / b) := Real.rpow_nonneg hx _
  rw [← h1] at h
  exact le_of_pow_le_pow_left₀ hb hy h

/-- Monotonicity of `t ↦ t^(a/b)` on `[0, ∞)`. -/
theorem rpowAB_mono (a b : ℕ) {x x' : ℝ} (hx : 0 ≤ x) (hxx : x ≤ x') :
    x ^ ((a : ℝ) / b) ≤ x' ^ ((a : ℝ) / b) :=
  Real.rpow_le_rpow hx hxx (by positivity)

/-- Integer cross-multiplied comparison of `(n/dn)^a` and `(yn/yd)^b`. -/
theorem powAB_div_le (a b : ℕ) {n dn yn yd : ℤ} (hdn : 0 < dn) (hyd : 0 < yd)
    (h : n ^ a * yd ^ b ≤ yn ^ b * dn ^ a) :
    ((n : ℝ) / dn) ^ a ≤ ((yn : ℝ) / yd) ^ b := by
  have hdn' : (0 : ℝ) < (dn : ℝ) ^ a := pow_pos (by exact_mod_cast hdn) a
  have hyd' : (0 : ℝ) < (yd : ℝ) ^ b := pow_pos (by exact_mod_cast hyd) b
  rw [div_pow, div_pow, div_le_div_iff₀ hdn' hyd']
  exact_mod_cast h

/-- Entry extraction from the lockstep check. -/
theorem rpowAB_check_entry (a b : ℕ) {sx sy : Int} {X Y : List (List Int)}
    (h : sparseRpowABCheckL22 a b sx sy X Y = true) :
    0 < sx ∧ 0 < sy ∧ MatShape X ∧ MatShape Y ∧
      ∀ i j, i < 11 → j < 11 →
        sparseRpowABEntryL22 a b sx sy (sparseLMat22 X i j) (sparseLMat22 Y i j) = true := by
  simp only [sparseRpowABCheckL22, Bool.and_eq_true, decide_eq_true_eq] at h
  obtain ⟨⟨⟨⟨hsx, hsy⟩, hX⟩, hY⟩, hall⟩ := h
  have hXs := matShape_of_check X hX
  have hYs := matShape_of_check Y hY
  refine ⟨hsx, hsy, hXs, hYs, fun i j hi hj => ?_⟩
  have hiX : i < X.length := by rw [hXs.1]; exact hi
  have hiY : i < Y.length := by rw [hYs.1]; exact hi
  have hrow := all_zipWith hall i hiX hiY
  have hjX : j < X[i].length := by rw [hXs.2 _ (List.getElem_mem hiX)]; exact hj
  have hjY : j < Y[i].length := by rw [hYs.2 _ (List.getElem_mem hiY)]; exact hj
  have he := all_zipWith hrow j hjX hjY
  unfold sparseLMat22
  rw [gd X [] hiX, gd Y [] hiY, gd X[i] 0 hjX, gd Y[i] 0 hjY]
  exact he

/-- Soundness of the list check for any nonnegative real matrix dominated by `X / sx`. -/
theorem sparseRpowABCheckL22_sound (a b : ℕ) (hb : b ≠ 0) {sx sy : Int}
    {X Y : List (List Int)}
    (h : sparseRpowABCheckL22 a b sx sy X Y = true) (M : Matrix (Fin 11) (Fin 11) ℝ)
    (hM0 : ∀ i j, 0 ≤ M i j)
    (hMX : ∀ i j : Fin 11, M i j ≤ (sparseLMat22 X i.val j.val : ℝ) / sx) :
    ∀ i j : Fin 11, M i j ^ ((a : ℝ) / b) ≤ (sparseLMat22 Y i.val j.val : ℝ) / sy := by
  obtain ⟨hsx, hsy, -, -, hent⟩ := rpowAB_check_entry a b h
  intro i j
  have he := hent i.val j.val i.isLt j.isLt
  simp only [sparseRpowABEntryL22, Bool.and_eq_true, decide_eq_true_eq] at he
  have hy : (0 : ℝ) ≤ (sparseLMat22 Y i.val j.val : ℝ) / sy :=
    div_nonneg (by exact_mod_cast he.1) (by exact_mod_cast hsy.le)
  have ha : M i j ^ a ≤ ((sparseLMat22 X i.val j.val : ℝ) / sx) ^ a :=
    pow_le_pow_left₀ (hM0 i j) (hMX i j) a
  exact rpowAB_le_of_pow a b hb (hM0 i j) hy (ha.trans (powAB_div_le a b hsx hsy he.2))

/-- The same bound phrased via `matrixEntrywiseRpow`. -/
theorem sparseRpowABCheckL22_entrywiseRpow (a b : ℕ) (hb : b ≠ 0) {sx sy : Int}
    {X Y : List (List Int)}
    (h : sparseRpowABCheckL22 a b sx sy X Y = true) (M : Matrix (Fin 11) (Fin 11) ℝ)
    (hM0 : ∀ i j, 0 ≤ M i j)
    (hMX : ∀ i j : Fin 11, M i j ≤ (sparseLMat22 X i.val j.val : ℝ) / sx) :
    ∀ i j : Fin 11,
      matrixEntrywiseRpow M ((a : ℝ) / b) i j ≤ (sparseLMat22 Y i.val j.val : ℝ) / sy :=
  sparseRpowABCheckL22_sound a b hb h M hM0 hMX

end SparseRpowABCert22
end Spin.Turbo

#print axioms Spin.Turbo.SparseRpowABCert22.rpowAB_le_of_pow
#print axioms Spin.Turbo.SparseRpowABCert22.sparseRpowABCheckL22_entrywiseRpow

/-- Kernel smoke test: `a/b = 7/10` on a tiny diagonal entry pattern. -/
example : Spin.Turbo.sparseRpowABCheckL22 7 10 1 1 (List.replicate 11 (List.replicate 11 1))
    (List.replicate 11 (List.replicate 11 1)) = true := by decide +kernel
