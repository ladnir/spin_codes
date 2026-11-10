import SpinCodes.Turbo.SparseListSound22
import SpinCodes.Turbo.ConditionalFractionalMatrix

/-!
# Integer certificates for entrywise `x ^ (4/5)` upper bounds

The kernel check compares `X^4 · sy^5 ≤ Y^5 · sx^4` entrywise by lockstep
list traversal (`zipWith`), never by indexed access. Soundness turns this
into `(M i j) ^ (4/5) ≤ Y_ij / sy` for any nonnegative real matrix `M`
dominated by `X / sx`.
-/
namespace Spin.Turbo

/-- One entry: `0 ≤ y` and `x^4 · sy^5 ≤ y^5 · sx^4`. -/
def sparseRpow45EntryL22 (sx sy : Int) (x y : Int) : Bool :=
  decide (0 ≤ y) && decide (x ^ 4 * sy ^ 5 ≤ y ^ 5 * sx ^ 4)

def sparseRpow45RowL22 (sx sy : Int) (xs ys : List Int) : Bool :=
  (List.zipWith (sparseRpow45EntryL22 sx sy) xs ys).all id

def sparseRpow45CheckL22 (sx sy : Int) (X Y : List (List Int)) : Bool :=
  decide (0 < sx) && decide (0 < sy) && sparseMatShapeL22 X && sparseMatShapeL22 Y &&
    (List.zipWith (sparseRpow45RowL22 sx sy) X Y).all id

end Spin.Turbo

noncomputable section
namespace Spin.Turbo
namespace SparseRpowCert22
open SparseListSound22

/-- `x^4 ≤ y^5` gives `x^(4/5) ≤ y` for nonnegative reals. -/
theorem rpow45_le_of_pow {x y : ℝ} (hx : 0 ≤ x) (hy : 0 ≤ y) (h : x ^ 4 ≤ y ^ 5) :
    x ^ ((4 : ℝ) / 5) ≤ y := by
  have h1 : (x ^ ((4 : ℝ) / 5)) ^ (5 : ℕ) = x ^ (4 : ℕ) := by
    rw [← Real.rpow_natCast (x ^ ((4 : ℝ) / 5)) 5, ← Real.rpow_mul hx,
      show (4 : ℝ) / 5 * ((5 : ℕ) : ℝ) = ((4 : ℕ) : ℝ) by norm_num, Real.rpow_natCast]
  have h0 : 0 ≤ x ^ ((4 : ℝ) / 5) := Real.rpow_nonneg hx _
  rw [← h1] at h
  exact le_of_pow_le_pow_left₀ (by norm_num : (5 : ℕ) ≠ 0) hy h

/-- Monotonicity of `t ↦ t^(4/5)` on `[0, ∞)`. -/
theorem rpow45_mono {x x' : ℝ} (hx : 0 ≤ x) (hxx : x ≤ x') :
    x ^ ((4 : ℝ) / 5) ≤ x' ^ ((4 : ℝ) / 5) :=
  Real.rpow_le_rpow hx hxx (by norm_num)

/-- Integer cross-multiplied comparison of `(n/dn)^4` and `(yn/yd)^5`. -/
theorem pow45_div_le {n dn yn yd : ℤ} (hdn : 0 < dn) (hyd : 0 < yd)
    (h : n ^ 4 * yd ^ 5 ≤ yn ^ 5 * dn ^ 4) :
    ((n : ℝ) / dn) ^ 4 ≤ ((yn : ℝ) / yd) ^ 5 := by
  have hdn' : (0 : ℝ) < (dn : ℝ) ^ 4 := pow_pos (by exact_mod_cast hdn) 4
  have hyd' : (0 : ℝ) < (yd : ℝ) ^ 5 := pow_pos (by exact_mod_cast hyd) 5
  rw [div_pow, div_pow, div_le_div_iff₀ hdn' hyd']
  exact_mod_cast h

/-- Scalar certificate: `(n/dn)^(4/5) ≤ yn/yd`. -/
theorem rpow45_scalar_le {n dn yn yd : ℤ} (hn : 0 ≤ n) (hdn : 0 < dn) (hyn : 0 ≤ yn)
    (hyd : 0 < yd) (h : n ^ 4 * yd ^ 5 ≤ yn ^ 5 * dn ^ 4) :
    ((n : ℝ) / dn) ^ ((4 : ℝ) / 5) ≤ (yn : ℝ) / yd := by
  have hx : (0 : ℝ) ≤ (n : ℝ) / dn :=
    div_nonneg (by exact_mod_cast hn) (by exact_mod_cast hdn.le)
  have hy : (0 : ℝ) ≤ (yn : ℝ) / yd :=
    div_nonneg (by exact_mod_cast hyn) (by exact_mod_cast hyd.le)
  exact rpow45_le_of_pow hx hy (pow45_div_le hdn hyd h)

/-- Entry extraction from the lockstep check. -/
theorem rpow45_check_entry {sx sy : Int} {X Y : List (List Int)}
    (h : sparseRpow45CheckL22 sx sy X Y = true) :
    0 < sx ∧ 0 < sy ∧ MatShape X ∧ MatShape Y ∧
      ∀ i j, i < 11 → j < 11 →
        sparseRpow45EntryL22 sx sy (sparseLMat22 X i j) (sparseLMat22 Y i j) = true := by
  simp only [sparseRpow45CheckL22, Bool.and_eq_true, decide_eq_true_eq] at h
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
theorem sparseRpow45CheckL22_sound {sx sy : Int} {X Y : List (List Int)}
    (h : sparseRpow45CheckL22 sx sy X Y = true) (M : Matrix (Fin 11) (Fin 11) ℝ)
    (hM0 : ∀ i j, 0 ≤ M i j)
    (hMX : ∀ i j : Fin 11, M i j ≤ (sparseLMat22 X i.val j.val : ℝ) / sx) :
    ∀ i j : Fin 11, M i j ^ ((4 : ℝ) / 5) ≤ (sparseLMat22 Y i.val j.val : ℝ) / sy := by
  obtain ⟨hsx, hsy, -, -, hent⟩ := rpow45_check_entry h
  intro i j
  have he := hent i.val j.val i.isLt j.isLt
  simp only [sparseRpow45EntryL22, Bool.and_eq_true, decide_eq_true_eq] at he
  have hy : (0 : ℝ) ≤ (sparseLMat22 Y i.val j.val : ℝ) / sy :=
    div_nonneg (by exact_mod_cast he.1) (by exact_mod_cast hsy.le)
  have h4 : M i j ^ 4 ≤ ((sparseLMat22 X i.val j.val : ℝ) / sx) ^ 4 :=
    pow_le_pow_left₀ (hM0 i j) (hMX i j) 4
  exact rpow45_le_of_pow (hM0 i j) hy (h4.trans (pow45_div_le hsx hsy he.2))

/-- The same bound phrased via `matrixEntrywiseRpow`. -/
theorem sparseRpow45CheckL22_entrywiseRpow {sx sy : Int} {X Y : List (List Int)}
    (h : sparseRpow45CheckL22 sx sy X Y = true) (M : Matrix (Fin 11) (Fin 11) ℝ)
    (hM0 : ∀ i j, 0 ≤ M i j)
    (hMX : ∀ i j : Fin 11, M i j ≤ (sparseLMat22 X i.val j.val : ℝ) / sx) :
    ∀ i j : Fin 11,
      matrixEntrywiseRpow M ((4 : ℝ) / 5) i j ≤ (sparseLMat22 Y i.val j.val : ℝ) / sy :=
  sparseRpow45CheckL22_sound h M hM0 hMX

end SparseRpowCert22
end Spin.Turbo

#print axioms Spin.Turbo.SparseRpowCert22.rpow45_le_of_pow
#print axioms Spin.Turbo.SparseRpowCert22.rpow45_mono
#print axioms Spin.Turbo.SparseRpowCert22.rpow45_scalar_le
#print axioms Spin.Turbo.SparseRpowCert22.sparseRpow45CheckL22_sound
#print axioms Spin.Turbo.SparseRpowCert22.sparseRpow45CheckL22_entrywiseRpow

/-- Kernel smoke test: all-zero 11x11 data. -/
example : Spin.Turbo.sparseRpow45CheckL22 1 1 (List.replicate 11 (List.replicate 11 0))
    (List.replicate 11 (List.replicate 11 0)) = true := by decide
