import SpinCodes.Structured.Domination

namespace Spin.Structured.ConcreteRoute

open Finset

/-- The fixed wiring from rows to regions. -/
def transpose {L b : ℕ} (rows : Fin L → Finset (Fin b)) : Fin b → Finset (Fin L) :=
  fun j => univ.filter (fun i => j ∈ rows i)

@[simp] theorem mem_transpose {L b : ℕ} (rows : Fin L → Finset (Fin b))
    (i : Fin L) (j : Fin b) : i ∈ transpose rows j ↔ j ∈ rows i := by
  simp [transpose]

@[simp] theorem transpose_transpose {L b : ℕ} (rows : Fin L → Finset (Fin b)) :
    transpose (transpose rows) = rows := by
  funext i
  ext j
  simp

/-- Transposition is a bijection of bit matrices. -/
def transposeEquiv (L b : ℕ) : (Fin L → Finset (Fin b)) ≃ (Fin b → Finset (Fin L)) where
  toFun := transpose
  invFun := transpose
  left_inv := transpose_transpose
  right_inv := transpose_transpose

theorem poissonBinom_eq_prod {n : ℕ} (p : Fin n → ℝ)
    (hp0 : ∀ i, 0 ≤ p i) (hp1 : ∀ i, p i ≤ 1) (S : Finset (Fin n)) :
    (poissonBinom p hp0 hp1).p S = ∏ i, if i ∈ S then p i else 1 - p i := by
  rw [poissonBinom_apply, Finset.prod_ite]
  congr 1 <;> congr 1 <;> ext i <;> simp

/-- Independent Bernoulli cells remain independent after transposition. -/
theorem transpose_density {L b : ℕ} (p : Fin L → ℝ)
    (hp0 : ∀ i, 0 ≤ p i) (hp1 : ∀ i, p i ≤ 1)
    (regions : Fin b → Finset (Fin L)) :
    (piPMF (fun i => poissonBinom (fun _ : Fin b => p i)
      (fun _ => hp0 i) (fun _ => hp1 i))).p (transpose regions) =
    (piPMF (fun _ : Fin b => poissonBinom p hp0 hp1)).p regions := by
  simp only [piPMF_apply, poissonBinom_eq_prod, mem_transpose]
  exact Finset.prod_comm

end Spin.Structured.ConcreteRoute

namespace Spin.Structured.ConcreteRoute

open Finset

/-- Number of input rows carrying at least one bit. -/
def activeRows {L b : ℕ} (rows : Fin L → Finset (Fin b)) : ℕ :=
  (univ.filter (fun i => rows i ≠ ∅)).card

/-- Total Hamming weight before routing. -/
def totalWeight {L b : ℕ} (rows : Fin L → Finset (Fin b)) : ℕ :=
  ∑ i, (rows i).card

/-- Mean density of the complete bit matrix. -/
noncomputable def density {L b : ℕ} (rows : Fin L → Finset (Fin b)) : ℝ :=
  (totalWeight rows : ℝ) / ((L : ℝ) * b)

theorem prod_rowCost {L b : ℕ} (rows : Fin L → Finset (Fin b)) :
    (∏ i, if rows i = ∅ then (1 : ℝ) else (b : ℝ) + 1) =
      ((b : ℝ) + 1) ^ activeRows rows := by
  rw [Finset.prod_ite]
  simp [activeRows]

theorem average_rowDensity {L b : ℕ} (rows : Fin L → Finset (Fin b)) :
    (∑ i, (rows i).card / (b : ℝ)) / (L : ℝ) = density rows := by
  rw [← Finset.sum_div, div_div]
  simp only [density, totalWeight, Nat.cast_sum]
  rw [mul_comm]

end Spin.Structured.ConcreteRoute


namespace Spin.Structured.ConcreteRoute

open Finset

theorem activeRows_le {L b : ℕ} (rows : Fin L → Finset (Fin b)) :
    activeRows rows ≤ L := by
  simpa [activeRows] using Finset.card_le_card
    (Finset.filter_subset (fun i => rows i ≠ ∅) (univ : Finset (Fin L)))

theorem totalWeight_le {L b : ℕ} (rows : Fin L → Finset (Fin b)) :
    totalWeight rows ≤ L * b := by
  calc totalWeight rows ≤ ∑ _ : Fin L, b := by
        apply Finset.sum_le_sum
        intro i _
        simpa using Finset.card_le_univ (rows i)
    _ = L * b := by simp

theorem totalWeight_transpose {L b : ℕ} (rows : Fin L → Finset (Fin b)) :
    totalWeight (transpose rows) = totalWeight rows := by
  have hc {n : ℕ} (S : Finset (Fin n)) :
      S.card = ∑ i : Fin n, if i ∈ S then 1 else 0 := by
    rw [← Finset.sum_filter]
    simp
  simp only [totalWeight, hc, mem_transpose]
  exact Finset.sum_comm

theorem density_nonneg {L b : ℕ} (rows : Fin L → Finset (Fin b)) :
    0 ≤ density rows := by
  unfold density
  positivity

theorem density_le_one {L b : ℕ} (hL : 0 < L) (hb : 0 < b)
    (rows : Fin L → Finset (Fin b)) : density rows ≤ 1 := by
  have h : (totalWeight rows : ℝ) ≤ (L : ℝ) * b := by
    exact_mod_cast totalWeight_le rows
  unfold density
  exact (div_le_one (by positivity)).mpr h

end Spin.Structured.ConcreteRoute
