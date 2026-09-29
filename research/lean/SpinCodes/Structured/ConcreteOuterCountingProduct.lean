import SpinCodes.Structured.ConcreteOuterCounting

noncomputable section
namespace Spin.Structured.ConcreteOuter
open Finset Routing ConcreteMarked

/-- Messages with exactly the specified active outer positions. -/
def activeMessages {L k : ℕ} (S : Finset (Fin L)) : Finset (Fin L → LocalMessage k) :=
  univ.filter (fun x => ∀ i, x i ≠ 0 ↔ i ∈ S)

/-- Counting all these messages and independently row-shuffling their encoded outputs. -/
def tupleCount {L k : ℕ} (seed : Seed k) (S : Finset (Fin L))
    (rows : Fin L → Finset (Fin (k * 24))) : ℝ :=
  ∑ x ∈ activeMessages S, (ConcreteRoute.rowLaw (rowSupports seed x)).p rows

/-- A position contributes the nonzero-message counting law or the deterministic zero row. -/
def positionCount {k : ℕ} (seed : Seed k) (active : Prop) [Decidable active]
    (row : Finset (Fin (k * 24))) : ℝ :=
  if active then shuffledCount seed row else (shuffleLaw ∅).p row

theorem positionCount_eq_sum {k : ℕ} (seed : Seed k) (active : Prop) [Decidable active]
    (row : Finset (Fin (k * 24))) :
    positionCount seed active row = ∑ x : LocalMessage k,
      if (x ≠ 0 ↔ active) then (shuffleLaw (support (encode seed x))).p row else 0 := by
  by_cases h : active
  · simp [positionCount, h, shuffledCount, Finset.sum_filter]
  · have hz : support (encode seed (0 : LocalMessage k)) = ∅ :=
      (support_eq_empty _).mpr (encode_zero seed)
    have hw : wtF (0 : Fin (k * 24) → Bool) = 0 := (wtF_eq_zero _).mpr rfl
    simp [positionCount, h, hw]

/-- For a fixed common seed, summing message tuples factors over their positions. -/
theorem tupleCount_eq_prod {L k : ℕ} (seed : Seed k) (S : Finset (Fin L))
    (rows : Fin L → Finset (Fin (k * 24))) :
    tupleCount seed S rows = ∏ i, positionCount seed (i ∈ S) (rows i) := by
  simp only [positionCount_eq_sum, prod_sum_eq_sum_prod, tupleCount, activeMessages,
    Finset.sum_filter, ConcreteRoute.rowLaw, piPMF_apply, rowSupports]
  apply Finset.sum_congr rfl
  intro x _
  by_cases h : ∀ i, x i ≠ 0 ↔ i ∈ S
  · rw [if_pos h]
    apply Finset.prod_congr rfl
    intro i _
    simp only [if_pos (h i)]
  · rw [if_neg h]
    obtain ⟨i, hi⟩ := not_forall.mp h
    symm
    exact Finset.prod_eq_zero (Finset.mem_univ i) (if_neg hi)

theorem shuffle_empty_eq_iid_zero (b : ℕ) :
    shuffleLaw (∅ : Finset (Fin b)) = iidBits b 0 (by norm_num) (by norm_num) := by
  ext S
  rw [shuffleLaw_apply, iidBits, poissonBinom_const_apply]
  by_cases h : S.card = 0
  · simp [h]
  · simp [h, zero_pow h]

theorem positionCount_nonneg {k : ℕ} (seed : Seed k) (active : Prop) [Decidable active]
    (row : Finset (Fin (k * 24))) : 0 ≤ positionCount seed active row := by
  unfold positionCount
  split_ifs
  · exact shuffledCount_nonneg seed row
  · exact (shuffleLaw ∅).nonneg row

theorem tupleCount_le_fairRows {L k : ℕ} (seed : Seed k) (S : Finset (Fin L)) {B : ℝ}
    (hB : ∀ w ≤ k * 24, (spectrum seed w : ℝ) ≤ B * ((k * 24).choose w : ℝ))
    (rows : Fin L → Finset (Fin (k * 24))) :
    tupleCount seed S rows ≤ ((2 : ℝ)^(k * 24) * B)^S.card * (fairRows S (k * 24)).p rows := by
  rw [tupleCount_eq_prod]
  have hc : (∏ i : Fin L, if i ∈ S then ((2 : ℝ)^(k * 24) * B) else 1) =
      ((2 : ℝ)^(k * 24) * B)^S.card := by
    rw [Finset.prod_ite]
    simp
  rw [← hc, fairRows, piPMF_apply, ← Finset.prod_mul_distrib]
  apply Finset.prod_le_prod₀ (fun i _ => positionCount_nonneg seed (i ∈ S) (rows i))
  intro i _
  by_cases hi : i ∈ S
  · simpa only [positionCount, activeProbability, if_pos hi] using shuffledCount_le_fair seed hB (rows i)
  · simp only [positionCount, activeProbability, if_neg hi, one_mul, shuffle_empty_eq_iid_zero]
    exact le_rfl

theorem tupleCount_le_max_fairRows {L k : ℕ} (seed : Seed k) (S : Finset (Fin L))
    (rows : Fin L → Finset (Fin (k * 24))) :
    tupleCount seed S rows ≤ ((2 : ℝ)^(k * 24) * maxShellRatio seed)^S.card *
      (fairRows S (k * 24)).p rows :=
  tupleCount_le_fairRows seed S (spectrum_le_max seed) rows

/-- The comparison transfers every nonnegative row statistic. -/
theorem activeMessages_expect_le {L k : ℕ} (seed : Seed k) (S : Finset (Fin L)) {B : ℝ}
    (hB : ∀ w ≤ k * 24, (spectrum seed w : ℝ) ≤ B * ((k * 24).choose w : ℝ))
    (F : (Fin L → Finset (Fin (k * 24))) → ℝ) (hF : ∀ rows, 0 ≤ F rows) :
    (∑ x ∈ activeMessages S, (ConcreteRoute.rowLaw (rowSupports seed x)).expect F) ≤
      ((2 : ℝ)^(k * 24) * B)^S.card * (fairRows S (k * 24)).expect F := by
  simp only [FinPMF.expect]
  rw [Finset.sum_comm]
  simp only [← Finset.sum_mul]
  change (∑ rows, tupleCount seed S rows * F rows) ≤ _
  rw [Finset.mul_sum]
  apply Finset.sum_le_sum
  intro rows _
  simpa only [mul_assoc] using mul_le_mul_of_nonneg_right
    (tupleCount_le_fairRows seed S hB rows) (hF rows)

end Spin.Structured.ConcreteOuter


