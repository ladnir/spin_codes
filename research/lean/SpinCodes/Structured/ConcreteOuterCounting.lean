import SpinCodes.Structured.ConcreteOuter
import SpinCodes.Structured.ConcreteMarkedFairRows
import SpinCodes.Selection

noncomputable section
namespace Spin.Structured.ConcreteOuter
open Finset Routing ConcreteMarked

/-- Counting nonzero messages, followed by one uniform coordinate shuffle. -/
def shuffledCount {k : ℕ} (seed : Seed k) (S : Finset (Fin (k * 24))) : ℝ :=
  ∑ x ∈ univ.filter (fun x : LocalMessage k => x ≠ 0),
    (shuffleLaw (support (encode seed x))).p S

theorem shuffledCount_nonneg {k : ℕ} (seed : Seed k) (S : Finset (Fin (k * 24))) :
    0 ≤ shuffledCount seed S :=
  Finset.sum_nonneg (fun x _ => (shuffleLaw _).nonneg S)

/-- Each support of weight w receives exactly A_seed(w)/choose(b,w) counting mass. -/
theorem shuffledCount_eq {k : ℕ} (seed : Seed k) (S : Finset (Fin (k * 24))) :
    shuffledCount seed S = (spectrum seed S.card : ℝ) / ((k * 24).choose S.card : ℝ) := by
  rw [spectrum_eq_sum, Finset.sum_div]
  apply Finset.sum_congr rfl
  intro x _
  rw [shuffleLaw_apply, support_card]
  by_cases h : wtF (encode seed x) = S.card
  · simp [h]
  · simp [h, Ne.symm h]

theorem fairBits_apply (b : ℕ) (S : Finset (Fin b)) :
    (iidBits b (1/2) (by norm_num) (by norm_num)).p S = 1 / (2 : ℝ)^b := by
  rw [iidBits, poissonBinom_const_apply]
  norm_num only [show (1 : ℝ) - 1/2 = 1/2 by norm_num]
  rw [← pow_add, Nat.add_sub_of_le (Routing.card_le_width S), div_pow, one_pow]

/-- A realized spectral bound gives an exact fair-row counting comparison. -/
theorem shuffledCount_le_fair {k : ℕ} (seed : Seed k) {B : ℝ}
    (hB : ∀ w ≤ k * 24, (spectrum seed w : ℝ) ≤ B * ((k * 24).choose w : ℝ))
    (S : Finset (Fin (k * 24))) :
    shuffledCount seed S ≤ ((2 : ℝ)^(k * 24) * B) *
      (iidBits (k * 24) (1/2) (by norm_num) (by norm_num)).p S := by
  rw [shuffledCount_eq, fairBits_apply]
  have hc := Routing.choose_card_pos S
  have hh := (div_le_iff₀ hc).mpr (hB S.card (Routing.card_le_width S))
  convert hh using 1 <;> field_simp

/-- Maximum realized count per support among the b+1 weight layers. -/
def maxShellRatio {k : ℕ} (seed : Seed k) : ℝ :=
  Finset.univ.sup' Finset.univ_nonempty
    (fun w : Fin (k * 24 + 1) => (spectrum seed w : ℝ) / ((k * 24).choose w : ℝ))

theorem shellRatio_le_max {k : ℕ} (seed : Seed k) (w : ℕ) (hw : w ≤ k * 24) :
    (spectrum seed w : ℝ) / ((k * 24).choose w : ℝ) ≤ maxShellRatio seed := by
  let j : Fin (k * 24 + 1) := ⟨w, by omega⟩
  unfold maxShellRatio
  exact Finset.le_sup' (fun v : Fin (k * 24 + 1) =>
    (spectrum seed v.val : ℝ) / ((k * 24).choose v.val : ℝ)) (Finset.mem_univ j)

theorem maxShellRatio_nonneg {k : ℕ} (seed : Seed k) : 0 ≤ maxShellRatio seed :=
  (by positivity : (0 : ℝ) ≤ (spectrum seed 0 : ℝ) / ((k * 24).choose 0 : ℝ)).trans
    (shellRatio_le_max seed 0 (by omega))

theorem spectrum_le_max {k : ℕ} (seed : Seed k) (w : ℕ) (hw : w ≤ k * 24) :
    (spectrum seed w : ℝ) ≤ maxShellRatio seed * ((k * 24).choose w : ℝ) :=
  (div_le_iff₀ (by exact_mod_cast Nat.choose_pos hw)).mp (shellRatio_le_max seed w hw)

theorem shuffledCount_le_max_fair {k : ℕ} (seed : Seed k) (S : Finset (Fin (k * 24))) :
    shuffledCount seed S ≤ ((2 : ℝ)^(k * 24) * maxShellRatio seed) *
      (iidBits (k * 24) (1/2) (by norm_num) (by norm_num)).p S :=
  shuffledCount_le_fair seed (spectrum_le_max seed) S

/-- The selected-outer event supplies a finite realized ratio bound from expected ratios. -/
theorem good_spectrum_bound {k : ℕ} (seed : Seed k) (W : Finset ℕ) {B : ℝ}
    (hB0 : 0 ≤ B)
    (hgood : Spin.Good (seedLaw k) spectrum (k * 24) W seed)
    (hB : ∀ w ∈ W, Spin.Abar (seedLaw k) spectrum w ≤
      B * ((k * 24).choose w : ℝ)) :
    ∀ w ≤ k * 24, (spectrum seed w : ℝ) ≤
      (((k * 24 : ℕ) : ℝ)^2 * B) * ((k * 24).choose w : ℝ) := by
  intro w hw
  by_cases hg : w ∈ W
  · calc (spectrum seed w : ℝ) ≤ ((k * 24 : ℕ) : ℝ)^2 * Spin.Abar (seedLaw k) spectrum w :=
        hgood.2 w hg
      _ ≤ ((k * 24 : ℕ) : ℝ)^2 * (B * ((k * 24).choose w : ℝ)) :=
        mul_le_mul_of_nonneg_left (hB w hg) (sq_nonneg _)
      _ = _ := by ring
  · have hz := hgood.1 w (by simp [Spin.badWeights, hw, hg])
    rw [hz, Nat.cast_zero]
    positivity

end Spin.Structured.ConcreteOuter


