import SpinCodes.Structured.Collatz

/-! The fixed-weight occupation matrix of Appendix `imt-finite-transfers`.

The entries below retain the lazy branch when every syndrome vanishes.
The column-action identities expose exactly the expressions used by the
sparse polynomial verifier. Cancellation bounds are explicit inputs; their
validity for the concrete maps is a separate finite-counting obligation.
-/

noncomputable section

namespace Spin.Imt

open Finset
variable {k : ℕ}

@[ext] theorem Coords.ext {c d : Coords k} (hZ : c.Z = d.Z)
    (hD : c.D = d.D) (hS : c.S = d.S) : c = d := by
  cases c; cases d; simp_all

namespace Occupation

/-- Moments and cancellation envelopes for one input weight. -/
structure RowData (k : ℕ) where
  live : ℝ
  zeroMoment : ℝ
  momentD : ℝ
  momentS : Fin k → ℝ
  cancelD : ℝ
  cancelS : Fin k → ℝ

/-- The paper's fixed-weight matrix, with shell cardinalities `count` and
nonzero-state count `M`. -/
def fixed (count : Fin k → ℝ) (M : ℝ) (r : RowData k) : Transfer k where
  rowZ := ⟨(1 - r.live) * r.zeroMoment, r.live * r.zeroMoment, fun _ => 0⟩
  rowD := ⟨r.cancelD / 2 + min r.momentD r.live / (2 * M),
    r.momentD / 2, fun h => count h * r.momentD / (2 * M)⟩
  rowS := fun i =>
    ⟨r.cancelS i / 2 + min (r.momentS i) r.live / (2 * M),
     if r.live = 0 then 0 else r.momentS i / 2,
     fun h => count h * r.momentS i / (2 * M) +
       if r.live = 0 ∧ h = i then r.momentS i / 2 else 0⟩

def liveAverage (count : Fin k → ℝ) (M : ℝ) (w : Coords k) : ℝ :=
  (∑ i, count i * w.S i) / M

theorem fixed_nonneg (count : Fin k → ℝ) (M : ℝ) (r : RowData k)
    (hc : ∀ i, 0 ≤ count i) (hM : 0 ≤ M)
    (hl : 0 ≤ r.live) (hl1 : r.live ≤ 1) (hz : 0 ≤ r.zeroMoment)
    (hmD : 0 ≤ r.momentD) (hmS : ∀ i, 0 ≤ r.momentS i)
    (hcD : 0 ≤ r.cancelD) (hcS : ∀ i, 0 ≤ r.cancelS i) :
    (fixed count M r).Nonneg := by
  have hd : 0 ≤ min r.momentD r.live := le_min hmD hl
  have hs : ∀ i, 0 ≤ min (r.momentS i) r.live := fun i => le_min (hmS i) hl
  dsimp [fixed, Transfer.Nonneg, Coords.Nonneg]
  refine ⟨⟨by positivity, by positivity, fun _ => le_rfl⟩,
    ⟨by positivity, by positivity, fun h => ?_⟩, fun i => ?_⟩
  · have := hc h
    positivity
  have hm := hmS i
  have hcancel := hcS i
  have hmin := hs i
  refine ⟨by positivity, ?_, fun h => ?_⟩
  · split <;> positivity
  · have := hc h
    apply add_nonneg (by positivity)
    split <;> positivity

lemma shell_pairing (count : Fin k → ℝ) (M a : ℝ) (w : Coords k) :
    (∑ i, count i * a / (2 * M) * w.S i) = a * liveAverage count M w / 2 := by
  simp only [liveAverage, div_eq_mul_inv, mul_inv_rev]
  calc
    _ = (∑ i, count i * w.S i) * (a * M⁻¹ * (2 : ℝ)⁻¹) := by
      rw [Finset.sum_mul]
      apply Finset.sum_congr rfl
      intro i _
      ring
    _ = _ := by ring

@[simp] theorem fixed_col_Z (count : Fin k → ℝ) (M : ℝ) (r : RowData k)
    (w : Coords k) :
    ((fixed count M r).applyCol w).Z =
      (1 - r.live) * r.zeroMoment * w.Z + r.live * r.zeroMoment * w.D := by
  simp [fixed, Transfer.applyCol, Coords.dot]

theorem fixed_col_D (count : Fin k → ℝ) (M : ℝ) (r : RowData k)
    (w : Coords k) :
    ((fixed count M r).applyCol w).D =
      (r.cancelD / 2 + min r.momentD r.live / (2 * M)) * w.Z +
        r.momentD * (w.D + liveAverage count M w) / 2 := by
  simp only [fixed, Transfer.applyCol, Coords.dot]
  rw [shell_pairing]
  ring

theorem fixed_col_S (count : Fin k → ℝ) (M : ℝ) (r : RowData k)
    (w : Coords k) (i : Fin k) :
    ((fixed count M r).applyCol w).S i =
      (r.cancelS i / 2 + min (r.momentS i) r.live / (2 * M)) * w.Z +
        r.momentS i * (liveAverage count M w +
          if r.live = 0 then w.S i else w.D) / 2 := by
  classical
  simp only [fixed, Transfer.applyCol, Coords.dot, add_mul, Finset.sum_add_distrib]
  rw [shell_pairing]
  by_cases h : r.live = 0
  · simp [h]
    ring
  · simp [h]
    ring

/-- A polynomial upper bound may replace a maximum or a minimum independently.
Only domination is required; the choice need not stay optimal on the interval. -/
theorem fixed_col_D_le (count : Fin k → ℝ) (M : ℝ) (r : RowData k)
    (w : Coords k) {a c f : ℝ} (hM : 0 ≤ M) (hwZ : 0 ≤ w.Z)
    (hwLive : 0 ≤ w.D + liveAverage count M w)
    (ha : r.momentD ≤ a) (hc : r.cancelD ≤ c)
    (hf : min r.momentD r.live ≤ f) :
    ((fixed count M r).applyCol w).D ≤
      (c / 2 + f / (2 * M)) * w.Z + a * (w.D + liveAverage count M w) / 2 := by
  rw [fixed_col_D]
  apply add_le_add
  · apply mul_le_mul_of_nonneg_right _ hwZ
    exact add_le_add (div_le_div_of_nonneg_right hc (by norm_num))
      (div_le_div_of_nonneg_right hf (by positivity))
  · exact div_le_div_of_nonneg_right (mul_le_mul_of_nonneg_right ha hwLive)
      (by norm_num)

theorem fixed_col_S_le (count : Fin k → ℝ) (M : ℝ) (r : RowData k)
    (w : Coords k) (i : Fin k) {c f : ℝ} (hM : 0 ≤ M) (hwZ : 0 ≤ w.Z)
    (hc : r.cancelS i ≤ c) (hf : min (r.momentS i) r.live ≤ f) :
    ((fixed count M r).applyCol w).S i ≤
      (c / 2 + f / (2 * M)) * w.Z + r.momentS i *
        (liveAverage count M w + if r.live = 0 then w.S i else w.D) / 2 := by
  rw [fixed_col_S]
  apply add_le_add _ le_rfl
  apply mul_le_mul_of_nonneg_right _ hwZ
  exact add_le_add (div_le_div_of_nonneg_right hc (by norm_num))
    (div_le_div_of_nonneg_right hf (by positivity))

/-- Binomial weight of an input with `j` occupied positions among `n`. -/
def probability (n j : ℕ) (β : ℝ) : ℝ :=
  (n.choose j : ℝ) * β ^ j * (1 - β) ^ (n - j)

theorem probability_nonneg (n j : ℕ) {β : ℝ} (h0 : 0 ≤ β) (h1 : β ≤ 1) :
    0 ≤ probability n j β := by
  unfold probability
  positivity

/-- The finite binomial mixture, including both endpoint input weights. -/
def matrix (n : ℕ) (count : Fin k → ℝ) (M β : ℝ)
    (r : Fin (n + 1) → RowData k) : Transfer k where
  rowZ := Coords.sum fun j : Fin (n + 1) => Coords.smul (probability n j β) (fixed count M (r j)).rowZ
  rowD := Coords.sum fun j : Fin (n + 1) => Coords.smul (probability n j β) (fixed count M (r j)).rowD
  rowS := fun i => Coords.sum fun j : Fin (n + 1) =>
    Coords.smul (probability n j β) ((fixed count M (r j)).rowS i)

lemma dot_sum_smul {ι : Type*} [Fintype ι] (p : ι → ℝ)
    (c : ι → Coords k) (w : Coords k) :
    (Coords.sum fun j => Coords.smul (p j) (c j)).dot w =
      ∑ j, p j * (c j).dot w := by
  simp only [Coords.sum, Coords.smul, Coords.dot, Finset.sum_mul,
    Finset.mul_sum, mul_add, Finset.sum_add_distrib, mul_assoc]
  rw [Finset.sum_comm]

theorem matrix_col (n : ℕ) (count : Fin k → ℝ) (M β : ℝ)
    (r : Fin (n + 1) → RowData k) (w : Coords k) :
    (matrix n count M β r).applyCol w = Coords.sum fun j : Fin (n + 1) =>
      Coords.smul (probability n j β) ((fixed count M (r j)).applyCol w) := by
  apply Coords.ext <;> simp only [matrix, Transfer.applyCol, Coords.sum, Coords.smul]
  · exact dot_sum_smul _ _ _
  · exact dot_sum_smul _ _ _
  · funext i
    exact dot_sum_smul _ _ _

/-- Every fixed-weight column bound survives the binomial mixture. -/
theorem matrix_col_le (n : ℕ) (count : Fin k → ℝ) (M β : ℝ)
    (r : Fin (n + 1) → RowData k) (w : Coords k)
    (u : Fin (n + 1) → Coords k) (h0 : 0 ≤ β) (h1 : β ≤ 1)
    (hu : ∀ j, ((fixed count M (r j)).applyCol w).le (u j)) :
    ((matrix n count M β r).applyCol w).le
      (Coords.sum fun j : Fin (n + 1) => Coords.smul (probability n j β) (u j)) := by
  rw [matrix_col]
  refine ⟨?_, ?_, fun i => ?_⟩ <;>
    simp only [Coords.sum, Coords.smul] <;>
    apply Finset.sum_le_sum <;> intro j _
  · exact mul_le_mul_of_nonneg_left (hu j).1 (probability_nonneg _ _ h0 h1)
  · exact mul_le_mul_of_nonneg_left (hu j).2.1 (probability_nonneg _ _ h0 h1)
  · exact mul_le_mul_of_nonneg_left ((hu j).2.2 i) (probability_nonneg _ _ h0 h1)

end Occupation
end Spin.Imt

