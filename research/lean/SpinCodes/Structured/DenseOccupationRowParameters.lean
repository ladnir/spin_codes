import SpinCodes.Structured.ConcreteRouteDomination

noncomputable section
namespace Spin.Structured.DenseOccupationFixed
open Finset ConcreteRoute

/-- The actual active row mean supplies every density parameter used by the dense boxes. -/
theorem row_parameters {L b : ℕ} (hL : 0 < L) (hb : 0 < b)
    (S : Finset (Fin L)) (hS : S.Nonempty) (rows : Fin L → Finset (Fin b))
    (hz : ∀ i, i ∉ S → rows i = ∅)
    (hr : ∀ i ∈ S, ((rows i).card:ℝ)/b ∈ Set.Icc (13/125) (112/125)) :
    ∃ x : ℝ, x ∈ Set.Icc (13/125) (112/125) ∧
      (∑ i ∈ S, ((rows i).card:ℝ)/b) = (S.card:ℝ)*x ∧
      density rows = ((S.card:ℝ)/L)*x ∧ activeRows rows = S.card ∧
      0 < totalWeight rows ∧ totalWeight rows < L*b := by
  have hqN : 0 < S.card := Finset.card_pos.mpr hS
  have hq : (0:ℝ) < S.card := by exact_mod_cast hqN
  have hLR : (0:ℝ) < L := by exact_mod_cast hL
  have hbR : (0:ℝ) < b := by exact_mod_cast hb
  let x : ℝ := (∑ i ∈ S, ((rows i).card:ℝ)/b)/S.card
  have hm : (∑ i ∈ S, ((rows i).card:ℝ)/b) = (S.card:ℝ)*x := by
    dsimp [x]
    field_simp
  have hlow : (S.card:ℝ)*(13/125) ≤ ∑ i ∈ S, ((rows i).card:ℝ)/b := by
    simpa only [Finset.sum_const, nsmul_eq_mul] using
      Finset.sum_le_sum (fun i hi => (hr i hi).1)
  have hhigh : (∑ i ∈ S, ((rows i).card:ℝ)/b) ≤ (S.card:ℝ)*(112/125) := by
    simpa only [Finset.sum_const, nsmul_eq_mul] using
      Finset.sum_le_sum (fun i hi => (hr i hi).2)
  have hx : x ∈ Set.Icc (13/125) (112/125) := by
    rw [hm] at hlow hhigh
    exact ⟨(mul_le_mul_iff_right₀ hq).mp hlow, (mul_le_mul_iff_right₀ hq).mp hhigh⟩
  have hs : (∑ i ∈ S, ((rows i).card:ℝ)/b) = ∑ i, ((rows i).card:ℝ)/b := by
    apply Finset.sum_subset (Finset.subset_univ S)
    intro i _ hi
    simp [hz i hi]
  have hd : density rows = ((S.card:ℝ)/L)*x := by
    rw [← average_rowDensity, ← hs, hm]
    ring
  have hactive : activeRows rows = S.card := by
    unfold activeRows
    congr 1
    ext i
    simp only [Finset.mem_filter, Finset.mem_univ, true_and]
    constructor
    · intro hi
      by_contra hn
      exact hi (hz i hn)
    · intro hi hempty
      have := (hr i hi).1
      simp only [hempty, Finset.card_empty, Nat.cast_zero, zero_div] at this
      norm_num at this
  have hQL : (S.card:ℝ) ≤ L := by
    have hh : S.card ≤ L := by simpa using Finset.card_le_univ S
    exact_mod_cast hh
  have hd0 : 0 < density rows := by rw [hd]; exact mul_pos (div_pos hq hLR) (by linarith [hx.1])
  have hd1 : density rows < 1 := by
    rw [hd]
    have hq1 : (S.card:ℝ)/L ≤ 1 := (div_le_one hLR).mpr hQL
    nlinarith [hx.1,hx.2, div_nonneg hq.le hLR.le]
  have ht0 : 0 < totalWeight rows := by
    have := (div_pos_iff_of_pos_right (mul_pos hLR hbR)).mp hd0
    exact_mod_cast this
  have ht1 : totalWeight rows < L*b := by
    have := (div_lt_one (mul_pos hLR hbR)).mp hd1
    exact_mod_cast this
  exact ⟨x,hx,hm,hd,hactive,ht0,ht1⟩

#print axioms row_parameters
end Spin.Structured.DenseOccupationFixed
