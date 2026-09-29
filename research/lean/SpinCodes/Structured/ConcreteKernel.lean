import SpinCodes.Structured.ConcreteFourier

/-! Kernel parity and distance from the actual feedback basis columns.
These facts justify the constant-weight packing bound without assuming a
kernel-distance certificate. -/

noncomputable section
namespace Spin.Structured.ConcreteMaps
open Finset PackedMap
open scoped symmDiff

theorem Cset_singleton (i : Fin 128) :
    Cset {i} = supportEquiv 19 (C (inputBasis i)) := by
  rw [Cset, supportEquiv_symm_singleton]
  rfl

theorem Cset_singleton_card (i : Fin 128) : (Cset {i}).card = 5 := by
  rw [Cset_singleton]
  change (support 19 (C (inputBasis i))).card = 5
  rw [← weight_eq_card_support, C_inputBasis_weight]

theorem Cset_singleton_injective : Function.Injective (fun i : Fin 128 => Cset {i}) := by
  intro i j h
  apply C_inputBasis_distinct
  apply (supportEquiv 19).injective
  simpa only [← Cset_singleton] using h

theorem CtransposeSet_univ : CtransposeSet univ = univ := by
  apply Finset.eq_univ_of_forall
  intro i
  have h := C_character_pairing univ {i}
  rw [Spin.chi_singleton_right] at h
  have hr : Spin.chi univ (Cset {i}) = -1 := by
    simp only [Spin.chi, Finset.univ_inter, Cset_singleton_card]
    norm_num
  rw [hr] at h
  by_contra hi
  simp only [if_neg hi] at h
  norm_num at h

theorem kernel_even {x : Finset (Fin 128)} (hx : Cset x = ∅) : Even x.card := by
  have h := C_character_pairing univ x
  rw [CtransposeSet_univ, hx] at h
  apply (neg_one_pow_eq_one_iff_even (by norm_num : (-1 : ℤ) ≠ 1)).mp
  simpa only [Spin.chi, Finset.univ_inter, Finset.inter_empty,
    Finset.card_empty, pow_zero] using h

theorem kernel_distance {x : Finset (Fin 128)} (hx : Cset x = ∅) (hne : x ≠ ∅) :
    4 ≤ x.card := by
  have hpos : 0 < x.card := Finset.card_pos.mpr (Finset.nonempty_iff_ne_empty.mpr hne)
  obtain ⟨k, hk⟩ := kernel_even hx
  by_contra hsmall
  have htwo : x.card = 2 := by omega
  obtain ⟨i, j, hij, rfl⟩ := Finset.card_eq_two.mp htwo
  rw [Spin.insert_eq_singleton_symmDiff i {j} (by simpa using hij),
    Cset_xor, Finset.symmDiff_eq_empty] at hx
  exact hij (Cset_singleton_injective hx)

theorem fiber_packing (j : ℕ) (q : Finset (Fin 19)) :
    (syndromeFiber j q).card * j.choose (j - 1) ≤ Nat.choose 128 (j - 1) := by
  obtain ⟨x, rfl⟩ := Cset_surjective q
  rw [syndromeFiber_eq]
  exact Spin.card_fiber_mul_choose_le dual j (d := 4) (by decide)
    (fun _ hu hn => kernel_distance ((orthTo_iff_kernel _).mp hu) hn) x

theorem fiber_packing_compl (j : ℕ) (q : Finset (Fin 19)) :
    (syndromeFiber j q).card * (128 - j).choose ((128 - j) - 1) ≤
      Nat.choose 128 ((128 - j) - 1) := by
  obtain ⟨x, rfl⟩ := Cset_surjective q
  rw [syndromeFiber_eq]
  exact Spin.card_fiber_mul_choose_le_compl dual j (d := 4) (by decide)
    (fun _ hu hn => kernel_distance ((orthTo_iff_kernel _).mp hu) hn) x

end Spin.Structured.ConcreteMaps
