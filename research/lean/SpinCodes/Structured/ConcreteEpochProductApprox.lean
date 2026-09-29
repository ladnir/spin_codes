import SpinCodes.Structured.ConcreteImpulseCoarseProduct

/-! Finite row-error calculus for products of nonnegative weighted kernels. -/
noncomputable section
namespace Spin.Structured.FiniteKernel
variable {α : Type*} [Fintype α] [DecidableEq α]

structure Substochastic (K : Matrix α α ℝ) : Prop where
  nonneg : ∀ q r, 0 ≤ K q r
  row_le : ∀ q, ∑ r, K q r ≤ 1

def rowAbs (K : Matrix α α ℝ) (q : α) : ℝ := ∑ r, |K q r|
def RowError (K L : Matrix α α ℝ) (ε : ℝ) : Prop := ∀ q, rowAbs (K - L) q ≤ ε

theorem Substochastic.one : Substochastic (1 : Matrix α α ℝ) := by
  constructor
  · intro q r; simp only [Matrix.one_apply]; split_ifs <;> norm_num
  · intro q; simp only [Matrix.one_apply, Finset.sum_ite_eq, Finset.mem_univ, ite_true, le_refl]

theorem Substochastic.mul {K L : Matrix α α ℝ} (hK : Substochastic K)
    (hL : Substochastic L) : Substochastic (K * L) := by
  constructor
  · intro q r
    exact Finset.sum_nonneg (fun s _ => mul_nonneg (hK.nonneg q s) (hL.nonneg s r))
  · intro q
    simp only [Matrix.mul_apply]
    rw [Finset.sum_comm]
    simp only [← Finset.mul_sum]
    calc ∑ s, K q s * ∑ r, L s r ≤ ∑ s, K q s * 1 :=
        Finset.sum_le_sum (fun s _ => mul_le_mul_of_nonneg_left (hL.row_le s) (hK.nonneg q s))
      _ ≤ 1 := by simpa only [mul_one] using hK.row_le q

theorem Substochastic.rowAbs_le {K : Matrix α α ℝ} (hK : Substochastic K) (q : α) :
    rowAbs K q ≤ 1 := by
  simpa only [rowAbs, abs_of_nonneg (hK.nonneg q _)] using hK.row_le q

theorem rowAbs_mul_le (K L : Matrix α α ℝ) (q : α) :
    rowAbs (K * L) q ≤ ∑ s, |K q s| * rowAbs L s := by
  unfold rowAbs
  simp only [Matrix.mul_apply]
  calc (∑ r, |∑ s, K q s * L s r|) ≤ ∑ r, ∑ s, |K q s * L s r| :=
      Finset.sum_le_sum (fun r _ => Finset.abs_sum_le_sum_abs _ _)
    _ = _ := by
      rw [Finset.sum_comm]
      simp only [abs_mul, Finset.mul_sum]

theorem RowError.refl (K : Matrix α α ℝ) : RowError K K 0 := by
  intro q
  simp only [rowAbs, Matrix.sub_apply, sub_self, abs_zero, Finset.sum_const_zero, le_refl]

theorem RowError.trans {K L M : Matrix α α ℝ} {ε δ : ℝ}
    (hKL : RowError K L ε) (hLM : RowError L M δ) : RowError K M (ε + δ) := by
  intro q
  calc
    rowAbs (K - M) q ≤ rowAbs (K - L) q + rowAbs (L - M) q := by
      simp only [rowAbs, Matrix.sub_apply, ← Finset.sum_add_distrib]
      exact Finset.sum_le_sum (fun r _ => abs_sub_le _ _ _)
    _ ≤ ε + δ := add_le_add (hKL q) (hLM q)

theorem RowError.mul_right {K L T : Matrix α α ℝ} {ε : ℝ}
    (h : RowError K L ε) (hT : Substochastic T) : RowError (K * T) (L * T) ε := by
  intro q
  rw [← Matrix.sub_mul]
  calc rowAbs ((K - L) * T) q ≤ ∑ s, |(K - L) q s| * rowAbs T s := rowAbs_mul_le _ _ _
    _ ≤ ∑ s, |(K - L) q s| * 1 :=
      Finset.sum_le_sum (fun s _ => mul_le_mul_of_nonneg_left (hT.rowAbs_le s) (abs_nonneg _))
    _ ≤ ε := by simpa only [mul_one, rowAbs] using h q

theorem RowError.mul_left {K L T : Matrix α α ℝ} {ε : ℝ}
    (h : RowError K L ε) (hT : Substochastic T) (hε : 0 ≤ ε) :
    RowError (T * K) (T * L) ε := by
  intro q
  rw [← Matrix.mul_sub]
  calc rowAbs (T * (K - L)) q ≤ ∑ s, |T q s| * rowAbs (K - L) s := rowAbs_mul_le _ _ _
    _ ≤ ∑ s, T q s * ε := by
      apply Finset.sum_le_sum
      intro s _
      rw [abs_of_nonneg (hT.nonneg q s)]
      exact mul_le_mul_of_nonneg_left (h s) (hT.nonneg q s)
    _ = (∑ s, T q s) * ε := (Finset.sum_mul _ _ _).symm
    _ ≤ 1 * ε := mul_le_mul_of_nonneg_right (hT.row_le q) hε
    _ = ε := one_mul _

theorem RowError.mul {K K' L L' : Matrix α α ℝ} {ε δ : ℝ}
    (hK : RowError K K' ε) (hL : RowError L L' δ)
    (hLsub : Substochastic L) (hKsub : Substochastic K') (hδ : 0 ≤ δ) :
    RowError (K * L) (K' * L') (ε + δ) :=
  (hK.mul_right hLsub).trans (hL.mul_left hKsub hδ)

theorem Substochastic.list_prod (Ks : List (Matrix α α ℝ))
    (h : ∀ K ∈ Ks, Substochastic K) : Substochastic Ks.prod := by
  induction Ks with
  | nil => exact Substochastic.one
  | cons K Ks ih =>
    exact (h K (by simp)).mul (ih (fun L hL => h L (by simp [hL])))

theorem product_rowError {n : Nat} (K L : Fin n → Matrix α α ℝ) (ε : Fin n → ℝ)
    (hK : ∀ i, Substochastic (K i)) (hL : ∀ i, Substochastic (L i))
    (hε : ∀ i, 0 ≤ ε i) (h : ∀ i, RowError (K i) (L i) (ε i)) :
    RowError (List.ofFn K).prod (List.ofFn L).prod (∑ i, ε i) := by
  induction n with
  | zero => simpa using RowError.refl (1 : Matrix α α ℝ)
  | succ n ih =>
    rw [List.ofFn_succ, List.ofFn_succ, List.prod_cons, List.prod_cons, Fin.sum_univ_succ]
    apply RowError.mul (h 0)
      (ih (fun i => K i.succ) (fun i => L i.succ) (fun i => ε i.succ)
        (fun i => hK i.succ) (fun i => hL i.succ) (fun i => hε i.succ) (fun i => h i.succ))
    · apply Substochastic.list_prod
      intro A hA
      obtain ⟨i, rfl⟩ := List.mem_ofFn.mp hA
      exact hK i.succ
    · exact hL 0
    · exact Finset.sum_nonneg (fun i _ => hε i.succ)

theorem rowError_of_entry {K L : Matrix α α ℝ} {ε : ℝ}
    (h : ∀ q r, |K q r - L q r| ≤ ε) : RowError K L ((Fintype.card α : ℝ) * ε) := by
  intro q
  calc rowAbs (K - L) q ≤ ∑ _r : α, ε := Finset.sum_le_sum (fun r _ => h q r)
    _ = _ := by simp only [Finset.sum_const, Finset.card_univ, nsmul_eq_mul]

end Spin.Structured.FiniteKernel
