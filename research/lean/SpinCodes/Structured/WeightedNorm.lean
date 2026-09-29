/-
The weighted row norm of `app:imt-fixed`.

    "Submultiplicativity, ‖D_σ‖_v ≤ 1, and ‖(I+uP_+)D_σ‖_v = m_v(u) prove
     (eq:imt-product-norm).  This reasoning requires nonnegative coefficients,
     not stochasticity of P_+."

With `v = (1, v₁)ᵀ` the weighted row norm of `A` is `max_i (1/v_i) ∑_j |A_ij| v_j`.
Rather than define that maximum, the bound is carried as a predicate,

    RowNormLe v₁ c A  :  ∑_j |A_ij| v_j ≤ c · v_i   for each row `i`,

which is exactly how it is used, and which avoids both a supremum and the
division by `v₁`.  The paper's

    m_v(u) = max{1 + uσv₁,  ur/v₁ + σ(1+u)}

is then the smallest `c` that works for `(I+uP_+)D_σ`, the two arguments of
the maximum being the two rows.

The three ingredients are proved here and composed along the interleaved
product `D_σ(I+u₁P_+)D_σ ⋯ (I+u_QP_+)D_σ`.  Nonnegativity of `u`, `r`, `σ`,
`v₁` is what lets the absolute values be dropped — stochasticity is never used.
-/
import SpinCodes.Structured.Simplex

namespace Spin

open Matrix

/-! ## The norm bound as a predicate -/

/-- `∑_j |A_ij| v_j ≤ c · v_i` for both rows, with `v = (1, v₁)`. -/
def RowNormLe (v1 c : ℝ) (A : Matrix (Fin 2) (Fin 2) ℝ) : Prop :=
  |A 0 0| + |A 0 1| * v1 ≤ c ∧ |A 1 0| + |A 1 1| * v1 ≤ c * v1

/-- One row of the product, against the norm bound on `B`. -/
private lemma row_bound {v1 b : ℝ} (hv1 : 0 ≤ v1) {B : Matrix (Fin 2) (Fin 2) ℝ}
    (hB : RowNormLe v1 b B) (x y : ℝ) :
    |x * B 0 0 + y * B 1 0| + |x * B 0 1 + y * B 1 1| * v1
      ≤ b * (|x| + |y| * v1) := by
  have t1 : |x * B 0 0 + y * B 1 0| ≤ |x| * |B 0 0| + |y| * |B 1 0| := by
    refine (abs_add_le _ _).trans ?_
    rw [abs_mul, abs_mul]
  have t2 : |x * B 0 1 + y * B 1 1| ≤ |x| * |B 0 1| + |y| * |B 1 1| := by
    refine (abs_add_le _ _).trans ?_
    rw [abs_mul, abs_mul]
  have w1 : |x| * (|B 0 0| + |B 0 1| * v1) ≤ |x| * b :=
    mul_le_mul_of_nonneg_left hB.1 (abs_nonneg _)
  have w2 : |y| * (|B 1 0| + |B 1 1| * v1) ≤ |y| * (b * v1) :=
    mul_le_mul_of_nonneg_left hB.2 (abs_nonneg _)
  nlinarith [mul_le_mul_of_nonneg_right t2 hv1, t1, w1, w2]

/-- **Submultiplicativity.** -/
theorem RowNormLe.mul {v1 a b : ℝ} (hv1 : 0 ≤ v1) (hb : 0 ≤ b)
    {A B : Matrix (Fin 2) (Fin 2) ℝ} (hA : RowNormLe v1 a A) (hB : RowNormLe v1 b B) :
    RowNormLe v1 (a * b) (A * B) := by
  have e00 : (A * B) 0 0 = A 0 0 * B 0 0 + A 0 1 * B 1 0 := by
    simp [Matrix.mul_apply, Fin.sum_univ_two]
  have e01 : (A * B) 0 1 = A 0 0 * B 0 1 + A 0 1 * B 1 1 := by
    simp [Matrix.mul_apply, Fin.sum_univ_two]
  have e10 : (A * B) 1 0 = A 1 0 * B 0 0 + A 1 1 * B 1 0 := by
    simp [Matrix.mul_apply, Fin.sum_univ_two]
  have e11 : (A * B) 1 1 = A 1 0 * B 0 1 + A 1 1 * B 1 1 := by
    simp [Matrix.mul_apply, Fin.sum_univ_two]
  constructor
  · rw [e00, e01]
    have h := row_bound hv1 hB (A 0 0) (A 0 1)
    have hfin : b * (|A 0 0| + |A 0 1| * v1) ≤ b * a :=
      mul_le_mul_of_nonneg_left hA.1 hb
    nlinarith
  · rw [e10, e11]
    have h := row_bound hv1 hB (A 1 0) (A 1 1)
    have hfin : b * (|A 1 0| + |A 1 1| * v1) ≤ b * (a * v1) :=
      mul_le_mul_of_nonneg_left hA.2 hb
    nlinarith

/-! ## The two factors -/

/-- `D_σ = diag(1, σ)`. -/
def Dsig (sig : ℝ) : Matrix (Fin 2) (Fin 2) ℝ := !![1, 0; 0, sig]

/-- `I + u P_+` for `P_+ = [[0,1],[r,1]]`. -/
def IuP (r u : ℝ) : Matrix (Fin 2) (Fin 2) ℝ := !![1, u; u * r, 1 + u]

/-- `I + u P_+` really is the identity plus `u` times `P_+`. -/
theorem IuP_eq (r u : ℝ) : IuP r u = 1 + u • Pplus r := by
  ext i j
  fin_cases i <;> fin_cases j <;> simp [IuP, Pplus, Matrix.one_apply] <;> ring

/-- The paper's `m_v(u)`, the larger of the two row requirements. -/
noncomputable def mv (sig v1 r u : ℝ) : ℝ :=
  max (1 + u * sig * v1) (u * r / v1 + sig * (1 + u))

/-- **`‖D_σ‖_v ≤ 1`.** -/
theorem rowNormLe_Dsig {sig v1 : ℝ} (hsig0 : 0 ≤ sig) (hsig1 : sig ≤ 1) (hv1 : 0 ≤ v1) :
    RowNormLe v1 1 (Dsig sig) := by
  constructor
  · show |(1 : ℝ)| + |(0 : ℝ)| * v1 ≤ 1
    simp
  · show |(0 : ℝ)| + |sig| * v1 ≤ 1 * v1
    rw [abs_zero, abs_of_nonneg hsig0]
    nlinarith

/-- **`‖(I+uP_+)D_σ‖_v = m_v(u)`**, in the form the chain uses. -/
theorem rowNormLe_step {sig v1 r u : ℝ} (hsig0 : 0 ≤ sig) (hv1 : 0 < v1)
    (hr : 0 ≤ r) (hu : 0 ≤ u) :
    RowNormLe v1 (mv sig v1 r u) (IuP r u * Dsig sig) := by
  have hprod : IuP r u * Dsig sig = !![1, u * sig; u * r, (1 + u) * sig] := by
    ext i j
    fin_cases i <;> fin_cases j <;>
      simp [IuP, Dsig, Matrix.mul_apply, Fin.sum_univ_two] <;> ring
  rw [hprod]
  constructor
  · show |(1 : ℝ)| + |u * sig| * v1 ≤ mv sig v1 r u
    rw [abs_one, abs_of_nonneg (by positivity : (0 : ℝ) ≤ u * sig)]
    exact le_max_left _ _
  · show |u * r| + |(1 + u) * sig| * v1 ≤ mv sig v1 r u * v1
    rw [abs_of_nonneg (by positivity : (0 : ℝ) ≤ u * r),
      abs_of_nonneg (by nlinarith : (0 : ℝ) ≤ (1 + u) * sig)]
    have hle : u * r / v1 + sig * (1 + u) ≤ mv sig v1 r u := le_max_right _ _
    have hmul := mul_le_mul_of_nonneg_right hle (le_of_lt hv1)
    have heq : (u * r / v1 + sig * (1 + u)) * v1 = u * r + (1 + u) * sig * v1 := by
      field_simp
    linarith [heq ▸ hmul]

/-! ## The interleaved product -/

/-- `(I+u₁P_+)D_σ ⋯ (I+u_QP_+)D_σ`. -/
noncomputable def pathProduct (sig r : ℝ) : List ℝ → Matrix (Fin 2) (Fin 2) ℝ
  | [] => 1
  | u :: us => (IuP r u * Dsig sig) * pathProduct sig r us

theorem rowNormLe_one {v1 : ℝ} :
    RowNormLe v1 1 (1 : Matrix (Fin 2) (Fin 2) ℝ) := by
  constructor
  · norm_num [Matrix.one_apply]
  · norm_num [Matrix.one_apply]

private lemma mv_nonneg {sig v1 r u : ℝ} (hsig0 : 0 ≤ sig) (hv1 : 0 < v1) (hu : 0 ≤ u) :
    0 ≤ mv sig v1 r u := by
  have : (0 : ℝ) ≤ 1 + u * sig * v1 := by positivity
  exact le_trans this (le_max_left _ _)

/-- **`eq:imt-product-norm`, the norm half.**  The interleaved product has
weighted row norm at most `∏ m_v(u_j)`. -/
theorem rowNormLe_pathProduct {sig v1 r : ℝ} (hsig0 : 0 ≤ sig) (hv1 : 0 < v1)
    (hr : 0 ≤ r) : ∀ us : List ℝ, (∀ u ∈ us, 0 ≤ u) →
      RowNormLe v1 ((us.map (mv sig v1 r)).prod) (pathProduct sig r us)
  | [], _ => by simpa [pathProduct] using rowNormLe_one (v1 := v1)
  | u :: us, hall => by
      have hu : 0 ≤ u := hall u (by simp)
      have hrest := rowNormLe_pathProduct hsig0 hv1 hr us
        (fun w hw => hall w (by simp [hw]))
      have hprodnn : 0 ≤ (us.map (mv sig v1 r)).prod := by
        refine List.prod_nonneg ?_
        intro x hx
        obtain ⟨w, hw, rfl⟩ := List.mem_map.mp hx
        exact mv_nonneg hsig0 hv1 (hall w (by simp [hw]))
      have := (rowNormLe_step hsig0 hv1 hr hu).mul hv1.le hprodnn hrest
      simpa [pathProduct, List.map_cons, List.prod_cons] using this

end Spin
