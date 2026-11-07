import SpinCodes.UniformMixer.Map22SubsetEquiv

/-!
# Submask parity for six-bit indices

Transport the Boolean-lattice interval parity through the checked
index-to-subset equivalence. This gives the exact cancellation identity
needed by the production ANF inversion, without a large finite sum
being reduced by the kernel.
-/

noncomputable section

namespace Spin.UniformMixer

theorem six_bit_submask_parity (j k : Fin 64) :
    (∑ m : Fin 64,
      if isSubmask k m ∧ isSubmask m j then (1 : ZMod 2) else 0) =
    if k = j then 1 else 0 := by
  have htransport :
      (∑ m : Fin 64,
        if isSubmask k m ∧ isSubmask m j then (1 : ZMod 2) else 0) =
      ∑ s : Finset (Fin 6),
        if bitPositions k ⊆ s ∧ s ⊆ bitPositions j then
          (1 : ZMod 2) else 0 := by
    apply Fintype.sum_equiv bitPositionsEquiv
    intro m
    change (if isSubmask k m ∧ isSubmask m j then (1 : ZMod 2) else 0) =
      if bitPositions k ⊆ bitPositions m ∧
        bitPositions m ⊆ bitPositions j then 1 else 0
    simp only [isSubmask_iff_subset]
  rw [htransport]
  have hfilter :
      (Finset.univ.filter (fun s : Finset (Fin 6) =>
        bitPositions k ⊆ s ∧ s ⊆ bitPositions j)) =
      Finset.Icc (bitPositions k) (bitPositions j) := by
    ext s
    simp only [Finset.mem_filter, Finset.mem_univ, true_and,
      Finset.mem_Icc]
  calc
    (∑ s : Finset (Fin 6),
        if bitPositions k ⊆ s ∧ s ⊆ bitPositions j then
          (1 : ZMod 2) else 0) =
      ∑ s ∈ Finset.Icc (bitPositions k) (bitPositions j),
        (1 : ZMod 2) := by
      simp only [Finset.sum_ite, Finset.sum_const_zero, add_zero]
      rw [hfilter]
    _ = if bitPositions k = bitPositions j then 1 else 0 :=
      finset_interval_parity (bitPositions k) (bitPositions j)
    _ = if k = j then 1 else 0 := by
      by_cases h : k = j
      · simp [h]
      · have hne : bitPositions k ≠ bitPositions j :=
          fun he => h (bitPositions_injective he)
        simp [h, hne]

end Spin.UniformMixer
