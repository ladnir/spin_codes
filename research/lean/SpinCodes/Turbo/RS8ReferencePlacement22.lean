import SpinCodes.Turbo.RS8ReferenceRegion22
import SpinCodes.Structured.ConcretePlacementLaw

/-!
# Complete-value placement law in one Turbo region

The old `ConcretePlacementLaw` counts uniform permutation images of a
selected set. The Turbo reference additionally carries a uniform four-bit
value at every selected position, including zero. We keep the value law
pointwise before using the old subset-inclusion calculation.
-/

noncomputable section

namespace Spin.Turbo

open Finset

def regionalPacketSupport22 {n : ℕ} (x : Fin n → Packet4) :
    Finset (Fin n) :=
  Finset.univ.filter (fun i => x i ≠ 0)

theorem uniformPacket4_apply (x : Packet4) :
    (FinPMF.uniform Packet4).p x = 1 / (16 : ℝ) := by
  change 1 / (Fintype.card Packet4 : ℝ) = _
  simp [Packet4]

/-- Before regional shuffling, an active position has a uniform four-bit
value and an inactive position is zero. This is a statement about each
complete regional word, not merely its support. -/
theorem turboReferenceRegionPre_apply
    {R : ℕ} (S : Finset (Fin (R * 16)))
    (x : Fin (R * 16) → Packet4) :
    (Spin.piPMF (fun i : Fin (R * 16) =>
      turboReferenceCellLaw22 S i)).p x =
      if regionalPacketSupport22 x ⊆ S then
        1 / (16 : ℝ) ^ S.card else 0 := by
  rw [Spin.piPMF_apply]
  by_cases hS : regionalPacketSupport22 x ⊆ S
  · rw [if_pos hS]
    have hterm (i : Fin (R * 16)) :
        (turboReferenceCellLaw22 S i).p (x i) =
          if i ∈ S then 1 / (16 : ℝ) else 1 := by
      by_cases hi : i ∈ S
      · simp only [turboReferenceCellLaw22, hi, ite_true]
        exact uniformPacket4_apply (x i)
      · have hx : x i = 0 := by
          by_contra hn
          exact hi (hS (Finset.mem_filter.mpr
            ⟨Finset.mem_univ _, hn⟩))
        simp [turboReferenceCellLaw22, hi,
          zeroPacketLaw22_apply, hx]
    simp_rw [hterm]
    rw [Finset.prod_ite]
    simp [div_eq_mul_inv]
  · rw [if_neg hS]
    have hw : ∃ i : Fin (R * 16),
        i ∈ regionalPacketSupport22 x ∧ i ∉ S := by
      by_contra hn
      apply hS
      intro i hi
      by_contra hnot
      exact hn ⟨i, hi, hnot⟩
    obtain ⟨i, hi, hnot⟩ := hw
    apply Finset.prod_eq_zero (Finset.mem_univ i)
    have hx : x i ≠ 0 := (Finset.mem_filter.mp hi).2
    simp [turboReferenceCellLaw22, hnot,
      zeroPacketLaw22_apply, hx]

def regionalValueShuffleEquiv22 {R : ℕ}
    (σ : Equiv.Perm (Fin (R * 16))) :
    (Fin (R * 16) → Packet4) ≃
      (Fin (R * 16) → Packet4) where
  toFun x i := x (σ i)
  invFun y i := y (σ.symm i)
  left_inv := by intro x; funext i; simp
  right_inv := by intro y; funext i; simp

theorem turboReferenceOneRegionLaw22_eq_perm_mixture
    (R : ℕ) (S : Finset (Fin (R * 16))) :
    turboReferenceOneRegionLaw22 R S =
      (FinPMF.uniform (Equiv.Perm (Fin (R * 16)))).bind
        (fun σ => (Spin.piPMF (fun i : Fin (R * 16) =>
          turboReferenceCellLaw22 S i)).map
            (regionalValueShuffleEquiv22 σ)) := by
  let P := Spin.piPMF (fun i : Fin (R * 16) =>
    turboReferenceCellLaw22 S i)
  let U := FinPMF.uniform (Equiv.Perm (Fin (R * 16)))
  apply FinPMF.eq_of_expect_eq
  intro f
  rw [turboReferenceOneRegionLaw22, FinPMF.expect_map,
    FinPMF.expect_prod, FinPMF.expect_bind]
  simp only [FinPMF.expect_map]
  unfold FinPMF.expect
  change (∑ x, P.p x * (∑ σ, U.p σ * f (regionalValueShuffleEquiv22 σ x))) =
    ∑ σ, U.p σ * (∑ x, P.p x * f (regionalValueShuffleEquiv22 σ x))
  calc
    _ = ∑ x, ∑ σ,
          P.p x * (U.p σ * f (regionalValueShuffleEquiv22 σ x)) := by
      simp only [Finset.mul_sum]
    _ = ∑ σ, ∑ x,
          P.p x * (U.p σ * f (regionalValueShuffleEquiv22 σ x)) :=
      Finset.sum_comm
    _ = _ := by
      apply Finset.sum_congr rfl
      intro σ _
      rw [Finset.mul_sum]
      apply Finset.sum_congr rfl
      intro x _
      ring

theorem turboReferenceOneRegionLaw22_apply_inverse
    (R : ℕ) (S : Finset (Fin (R * 16)))
    (y : Fin (R * 16) → Packet4) :
    (turboReferenceOneRegionLaw22 R S).p y =
      (FinPMF.uniform (Equiv.Perm (Fin (R * 16)))).expect
        (fun σ =>
          (Spin.piPMF (fun i : Fin (R * 16) =>
            turboReferenceCellLaw22 S i)).p
              ((regionalValueShuffleEquiv22 σ).symm y)) := by
  rw [turboReferenceOneRegionLaw22_eq_perm_mixture, FinPMF.bind_p]
  simp_rw [FinPMF.map_equiv_apply]
  rfl

theorem regionalValueShuffle_subset_iff {R : ℕ}
    (S : Finset (Fin (R * 16)))
    (σ : Equiv.Perm (Fin (R * 16)))
    (y : Fin (R * 16) → Packet4) :
    regionalPacketSupport22 ((regionalValueShuffleEquiv22 σ).symm y) ⊆ S ↔
      regionalPacketSupport22 y ⊆
        Spin.Structured.Routing.shuffleSupport σ S := by
  constructor
  · intro h i hi
    apply (Spin.Structured.Routing.mem_shuffleSupport σ S i).mpr
    apply h
    apply Finset.mem_filter.mpr
    refine ⟨Finset.mem_univ _, ?_⟩
    change y (σ.symm (σ i)) ≠ 0
    simpa using (Finset.mem_filter.mp hi).2
  · intro h i hi
    have hj : σ.symm i ∈ regionalPacketSupport22 y := by
      apply Finset.mem_filter.mpr
      refine ⟨Finset.mem_univ _, ?_⟩
      exact (Finset.mem_filter.mp hi).2
    have hs := (Spin.Structured.Routing.mem_shuffleSupport σ S
      (σ.symm i)).mp (h hj)
    simpa using hs

theorem finite_expect_if_const
    {α : Type*} [Fintype α] (P : FinPMF α)
    (E : α → Prop) [DecidablePred E] (c : ℝ) :
    P.expect (fun x => if E x then c else 0) = P.prob E * c := by
  rw [FinPMF.prob_eq_expect_indicator]
  unfold FinPMF.expect
  rw [Finset.sum_mul]
  apply Finset.sum_congr rfl
  intro x _
  by_cases h : E x <;> simp [h]

/-- The complete regional word is obtained by uniformly placing the
`|S|` selected group positions, then giving each selected position a
uniform four-bit value (zero included). This uses the actual uniform
regional permutation, not an independent slot model. -/
theorem turboReferenceOneRegionLaw22_apply_subset
    (R : ℕ) (S : Finset (Fin (R * 16)))
    (y : Fin (R * 16) → Packet4) :
    (turboReferenceOneRegionLaw22 R S).p y =
      (Spin.Structured.Routing.shuffleLaw S).prob
        (fun T => regionalPacketSupport22 y ⊆ T) *
          (1 / (16 : ℝ) ^ S.card) := by
  rw [turboReferenceOneRegionLaw22_apply_inverse]
  simp_rw [turboReferenceRegionPre_apply]
  simp_rw [regionalValueShuffle_subset_iff]
  rw [finite_expect_if_const]
  rw [← FinPMF.map_prob
    (FinPMF.uniform (Equiv.Perm (Fin (R * 16))))
    (fun σ => Spin.Structured.Routing.shuffleSupport σ S)
    (fun T => regionalPacketSupport22 y ⊆ T)]
  rw [Spin.Structured.Routing.uniform_permutation_shuffleLaw]

theorem turboReferenceOneRegionLaw22_apply
    (R : ℕ) (S : Finset (Fin (R * 16)))
    (y : Fin (R * 16) → Packet4) :
    (turboReferenceOneRegionLaw22 R S).p y =
      if (regionalPacketSupport22 y).card ≤ S.card then
        (Nat.choose (R * 16 - (regionalPacketSupport22 y).card)
          (S.card - (regionalPacketSupport22 y).card) : ℝ) /
          ((Nat.choose (R * 16) S.card : ℝ) * (16 : ℝ) ^ S.card)
      else 0 := by
  rw [turboReferenceOneRegionLaw22_apply_subset]
  by_cases hk : (regionalPacketSupport22 y).card ≤ S.card
  · rw [if_pos hk,
      Spin.Structured.Placement.shuffle_prob_subset S
        (regionalPacketSupport22 y) hk]
    ring
  · have hlt : S.card < (regionalPacketSupport22 y).card := by omega
    rw [if_neg hk,
      Spin.Structured.Placement.shuffle_prob_subset_large S
        (regionalPacketSupport22 y) hlt]
    ring

#print axioms turboReferenceRegionPre_apply
#print axioms turboReferenceOneRegionLaw22_apply_inverse
#print axioms turboReferenceOneRegionLaw22_apply_subset
#print axioms turboReferenceOneRegionLaw22_apply

end Spin.Turbo
