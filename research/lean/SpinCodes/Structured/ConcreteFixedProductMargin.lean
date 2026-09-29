import SpinCodes.Structured.ConcreteFixedRegionComposition

/-! Nonnegative weighted kernel bounds compose without accumulating additive error. -/
noncomputable section
namespace Spin.Structured.FiniteKernel
open Finset
variable {α : Type*} [Fintype α] [DecidableEq α]

structure WeightedBound (w : α → ℝ) (c : ℝ) (K : Matrix α α ℝ) : Prop where
  nonneg : ∀ q r, 0 ≤ K q r
  bound : ∀ q, ∑ r, K q r * w r ≤ c * w q

theorem WeightedBound.one (w : α → ℝ) : WeightedBound w 1 (1 : Matrix α α ℝ) := by
  constructor
  · intro q r; simp only [Matrix.one_apply]; split_ifs <;> norm_num
  · intro q; simp [Matrix.one_apply]

theorem WeightedBound.mul {w : α → ℝ} {a b : ℝ} {K L : Matrix α α ℝ}
    (hK : WeightedBound w a K) (hL : WeightedBound w b L) (hb : 0 ≤ b) :
    WeightedBound w (a * b) (K * L) := by
  constructor
  · intro q r; exact sum_nonneg (fun s _ => mul_nonneg (hK.nonneg q s) (hL.nonneg s r))
  · intro q
    simp only [Matrix.mul_apply, sum_mul]
    rw [sum_comm]
    simp only [mul_assoc, ← mul_sum]
    calc
      _ ≤ ∑ s, K q s * (b * w s) := sum_le_sum (fun s _ =>
        mul_le_mul_of_nonneg_left (hL.bound s) (hK.nonneg q s))
      _ = b * ∑ s, K q s * w s := by rw [mul_sum]; apply sum_congr rfl; intro s _; ring
      _ ≤ b * (a * w q) := mul_le_mul_of_nonneg_left (hK.bound q) hb
      _ = _ := by ring

theorem WeightedBound.ofFn {b : Nat} {w : α → ℝ} {c : ℝ} (hc : 0 ≤ c)
    (K : Fin b → Matrix α α ℝ) (hK : ∀ i, WeightedBound w c (K i)) :
    WeightedBound w (c ^ b) (List.ofFn K).prod := by
  induction b with
  | zero => simpa using WeightedBound.one w
  | succ b ih =>
    rw [List.ofFn_succ, List.prod_cons, pow_succ']
    exact (hK 0).mul (ih (fun i => K i.succ) (fun i => hK i.succ)) (pow_nonneg hc _)

theorem WeightedBound.total {w : α → ℝ} {c m : ℝ} {K : Matrix α α ℝ}
    (hK : WeightedBound w c K) (hm : 0 < m) (hw : ∀ r, m ≤ w r) (q : α) :
    ∑ r, K q r ≤ c * w q / m := by
  apply (le_div_iff₀ hm).mpr
  calc
    _ = ∑ r, K q r * m := sum_mul _ _ _
    _ ≤ ∑ r, K q r * w r := sum_le_sum (fun r _ => mul_le_mul_of_nonneg_left (hw r) (hK.nonneg q r))
    _ ≤ c * w q := hK.bound q

end Spin.Structured.FiniteKernel
namespace Spin.Structured.Placement
open Finset ConcreteEncoder FiniteKernel Filter

theorem actual_fugacity_product_margin {Q : Nat} {θ v c δ : ℝ} (hθ : 0 ≤ θ) (hv : 0 < v)
    (hc : 0 ≤ c) (hδ : 0 < δ) (u : Fin Q → ℝ) (hu : ∀ i, 0 ≤ u i)
    (hK : Spin.RowNormLe v c (fugacityContinuum θ u)) :
    ∀ᶠ R : Nat in atTop, ∀ b : Nat, ∀ marks : Fin b → (Fin Q ↪ Fin (128 * R)),
      WeightedBound (stateWeight v) ((c + δ) ^ b) (List.ofFn (fun j => fugacityKernel θ (marks j) u)).prod := by
  filter_upwards [actual_fugacity_weighted_margin hθ hv hδ u hu hK] with R hR
  intro b marks
  apply WeightedBound.ofFn (by linarith)
  intro j
  exact ⟨fugacityKernel_nonneg θ (marks j) u hu, hR (marks j)⟩

theorem actual_fugacity_product_total {Q : Nat} {θ v c δ : ℝ} (hθ : 0 ≤ θ) (hv : 0 < v)
    (hc : 0 ≤ c) (hδ : 0 < δ) (u : Fin Q → ℝ) (hu : ∀ i, 0 ≤ u i)
    (hK : Spin.RowNormLe v c (fugacityContinuum θ u)) :
    ∀ᶠ R : Nat in atTop, ∀ b : Nat, ∀ marks : Fin b → (Fin Q ↪ Fin (128 * R)),
      ∑ r, (List.ofFn (fun j => fugacityKernel θ (marks j) u)).prod ∅ r ≤
        (c + δ) ^ b / min 1 v := by
  filter_upwards [actual_fugacity_product_margin hθ hv hc hδ u hu hK] with R hR
  intro b marks
  have hw (q : State) : min 1 v ≤ stateWeight v q := by
    unfold stateWeight
    split_ifs
    · exact min_le_left _ _
    · exact min_le_right _ _
  simpa [stateWeight] using
    (hR b marks).total (lt_min (by norm_num) hv) hw ∅

end Spin.Structured.Placement

