import SpinCodes.Structured.ConcreteShufflePermutation

/-! Exact inclusion probabilities for the actual uniformly shuffled region support. -/
noncomputable section
namespace Spin.Structured.Placement
open Finset Routing

def placementLayer (n k : Nat) : Finset (Finset (Fin n)) := univ.powersetCard k

theorem shuffle_prob_eq_count {n : Nat} (S : Finset (Fin n))
    (E : Finset (Fin n) → Prop) [DecidablePred E] :
    (shuffleLaw S).prob E =
      (((placementLayer n S.card).filter E).card : ℝ) / (n.choose S.card : ℝ) := by
  classical
  unfold Spin.FinPMF.prob
  simp only [shuffleLaw_apply]
  rw [← Finset.sum_filter]
  have he : (univ.filter E).filter (fun T => T.card = S.card) =
      (placementLayer n S.card).filter E := by
    ext T
    simp [placementLayer, and_comm]
  rw [he, Finset.sum_const, nsmul_eq_mul]
  ring

theorem shuffle_prob_subset {n : Nat} (S U : Finset (Fin n)) (hU : U.card ≤ S.card) :
    (shuffleLaw S).prob (fun T => U ⊆ T) =
      ((n - U.card).choose (S.card - U.card) : ℝ) / (n.choose S.card : ℝ) := by
  rw [shuffle_prob_eq_count]
  unfold placementLayer
  rw [Finset.card_filter_powersetCard_subset U univ S.card (Finset.subset_univ _) hU,
    Finset.card_univ, Fintype.card_fin]

theorem shuffle_prob_subset_large {n : Nat} (S U : Finset (Fin n)) (hU : S.card < U.card) :
    (shuffleLaw S).prob (fun T => U ⊆ T) = 0 := by
  rw [shuffle_prob_eq_count]
  have he : (placementLayer n S.card).filter (fun T => U ⊆ T) = ∅ := by
    apply Finset.filter_eq_empty_iff.mpr
    intro T hT hUT
    have ht := (Finset.mem_powersetCard.mp hT).2
    have hu := Finset.card_le_card hUT
    omega
  rw [he, Finset.card_empty, Nat.cast_zero, zero_div]

theorem choose_ratio {n k s : Nat} (hsk : s ≤ k) (hkn : k ≤ n) :
    ((n - s).choose (k - s) : ℝ) / (n.choose k : ℝ) =
      (k.choose s : ℝ) / (n.choose s : ℝ) := by
  have hn : (n.choose k : ℝ) ≠ 0 := by exact_mod_cast Nat.choose_ne_zero hkn
  have hs : (n.choose s : ℝ) ≠ 0 := by exact_mod_cast Nat.choose_ne_zero (hsk.trans hkn)
  have hi : (n.choose k : ℝ) * (k.choose s : ℝ) =
      (n.choose s : ℝ) * ((n - s).choose (k - s) : ℝ) := by
    exact_mod_cast Nat.choose_mul (n := n) hsk
  apply (div_eq_div_iff hn hs).mpr
  nlinarith

theorem shuffle_prob_mem {n : Nat} (S : Finset (Fin n)) (x : Fin n) :
    (shuffleLaw S).prob (fun T => x ∈ T) = (S.card : ℝ) / n := by
  have he : (shuffleLaw S).prob (fun T => x ∈ T) =
      (shuffleLaw S).prob (fun T => ({x} : Finset _) ⊆ T) :=
    Spin.FinPMF.prob_congr _ (fun _ => Finset.singleton_subset_iff.symm)
  rw [he]
  by_cases hk : 1 ≤ S.card
  · rw [shuffle_prob_subset S {x} (by simpa using hk), Finset.card_singleton,
      choose_ratio hk (card_le_width S), Nat.choose_one_right, Nat.choose_one_right]
  · have hk0 : S.card = 0 := by omega
    rw [shuffle_prob_subset_large S {x} (by simp [hk0]), hk0, Nat.cast_zero, zero_div]

theorem shuffle_prob_pair {n : Nat} (S : Finset (Fin n)) (x y : Fin n) (hxy : x ≠ y) :
    (shuffleLaw S).prob (fun T => x ∈ T ∧ y ∈ T) =
      ((S.card : ℝ) * (S.card - 1)) / ((n : ℝ) * (n - 1)) := by
  have he : (shuffleLaw S).prob (fun T => x ∈ T ∧ y ∈ T) =
      (shuffleLaw S).prob (fun T => ({x, y} : Finset _) ⊆ T) := by
    apply Spin.FinPMF.prob_congr
    intro T
    simp only [Finset.insert_subset_iff, Finset.singleton_subset_iff]
  have hp : ({x, y} : Finset (Fin n)).card = 2 := by simp [hxy]
  rw [he]
  by_cases hk : 2 ≤ S.card
  · rw [shuffle_prob_subset S {x, y} (by simpa [hp] using hk), hp,
      choose_ratio hk (card_le_width S), Nat.cast_choose_two, Nat.cast_choose_two]
    exact div_div_div_cancel_right₀ (by norm_num : (2 : ℝ) ≠ 0) _ _
  · rw [shuffle_prob_subset_large S {x, y} (by omega)]
    have hk01 : S.card = 0 ∨ S.card = 1 := by omega
    rcases hk01 with h | h <;> simp [h]

/-- The formulas apply to the actual uniform region permutation, not a
replacement iid location law. -/
theorem permutation_prob_pair {n : Nat} (S : Finset (Fin n)) (x y : Fin n) (hxy : x ≠ y) :
    (Spin.FinPMF.uniform (Equiv.Perm (Fin n))).prob
      (fun σ => x ∈ shuffleSupport σ S ∧ y ∈ shuffleSupport σ S) =
      ((S.card : ℝ) * (S.card - 1)) / ((n : ℝ) * (n - 1)) := by
  rw [← Spin.FinPMF.map_prob (Spin.FinPMF.uniform (Equiv.Perm (Fin n)))
    (fun σ => shuffleSupport σ S) (fun T => x ∈ T ∧ y ∈ T), uniform_permutation_shuffleLaw]
  exact shuffle_prob_pair S x y hxy

end Spin.Structured.Placement
