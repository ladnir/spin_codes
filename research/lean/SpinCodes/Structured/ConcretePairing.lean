import SpinCodes.Structured.ConcreteMaps
import SpinCodes.Structured.Parseval

/-! The actual feedback map and its listed transpose are adjoint for the
binary character pairing. This identifies the concrete kernel and syndrome
fibers with the dual-code definitions used by Fourier inversion. -/

namespace Spin.Structured.PackedMap

theorem pack_testBit (bs : List Bool) (i : ℕ) : (pack bs).testBit i = bs.getD i false := by
  induction bs generalizing i with
  | nil => simp [pack]
  | cons b bs ih =>
    cases i with
    | zero => cases b <;> simp [pack, Nat.testBit_zero]
    | succ i =>
      cases b <;> simp [pack, Nat.testBit_succ, Nat.add_div, ih]

theorem transpose_bit (rs : List ℕ) {width i j : ℕ} (hi : i < width) (hj : j < rs.length) :
    ((transpose width rs).getD i 0).testBit j = (rs.getD j 0).testBit i := by
  rw [transpose, List.getD_eq_getElem _ _ (by simpa using hi)]
  simp only [List.getElem_map, List.getElem_range, pack_testBit]
  rw [List.getD_eq_getElem _ _ (by simpa using hj), List.getD_eq_getElem _ _ hj]
  simp only [List.getElem_map]

theorem support_onehot {n : ℕ} (i : Fin n) : support n (2 ^ i.val) = {i} := by
  ext j
  simp [support, Nat.testBit_two_pow, Fin.ext_iff, eq_comm]

theorem supportEquiv_symm_singleton {n : ℕ} (i : Fin n) :
    (supportEquiv n).symm {i} =
      ⟨2 ^ i.val, Nat.pow_lt_pow_right (by decide) i.isLt⟩ := by
  apply (supportEquiv n).injective
  rw [Equiv.apply_symm_apply]
  exact (support_onehot i).symm

end Spin.Structured.PackedMap

namespace Spin
open Finset
open scoped symmDiff

theorem chi_singleton_right {n : ℕ} (x : Finset (Fin n)) (i : Fin n) :
    chi x {i} = if i ∈ x then -1 else 1 := by
  by_cases h : i ∈ x <;> simp [chi, h]

theorem chi_singleton_left {n : ℕ} (i : Fin n) (x : Finset (Fin n)) :
    chi {i} x = if i ∈ x then -1 else 1 := by
  unfold chi
  rw [Finset.inter_comm]
  exact chi_singleton_right x i

theorem insert_eq_singleton_symmDiff {n : ℕ} (i : Fin n) (x : Finset (Fin n)) (hi : i ∉ x) :
    insert i x = {i} ∆ x := by
  ext j
  simp only [Finset.mem_insert, Finset.mem_symmDiff, Finset.mem_singleton]
  by_cases h : j = i
  · subst j
    simp [hi]
  · simp [h]

/-- Binary linear maps are adjoint when their basis incidences are transposes. -/
theorem chi_adjoint_of_basis {s t : ℕ}
    (F : Finset (Fin s) → Finset (Fin t)) (G : Finset (Fin t) → Finset (Fin s))
    (hF0 : F ∅ = ∅) (hG0 : G ∅ = ∅)
    (hF : ∀ q p, F (q ∆ p) = F q ∆ F p)
    (hG : ∀ x y, G (x ∆ y) = G x ∆ G y)
    (hb : ∀ (j : Fin s) (i : Fin t), i ∈ F {j} ↔ j ∈ G {i})
    (q : Finset (Fin s)) (x : Finset (Fin t)) : chi (F q) x = chi q (G x) := by
  have hsingle (j : Fin s) (x : Finset (Fin t)) : chi (F {j}) x = chi {j} (G x) := by
    induction x using Finset.induction_on with
    | empty => simp [hG0, chi]
    | @insert i x hi ih =>
      rw [insert_eq_singleton_symmDiff i x hi, hG, chi_symmDiff_right, chi_symmDiff_right, ih]
      simp only [chi_singleton_right, chi_singleton_left, hb j i]
  induction q using Finset.induction_on with
  | empty => simp [hF0, chi]
  | @insert j q hj ih =>
    rw [insert_eq_singleton_symmDiff j q hj, hF, chi_symmDiff_left, chi_symmDiff_left,
      hsingle, ih]

end Spin

namespace Spin.Structured.ConcreteMaps
open Finset PackedMap
open scoped symmDiff

theorem C_basis_pairing (j : Fin 19) (i : Fin 128) :
    i ∈ CtransposeSet {j} ↔ j ∈ Cset {i} := by
  simp only [CtransposeSet, Cset, supportEquiv_symm_singleton]
  change i ∈ support 128 (eval cTransposeRows (2 ^ j.val)) ↔
    j ∈ support 19 (eval cRows (2 ^ i.val))
  simp only [support, Finset.mem_filter, Finset.mem_univ, true_and, eval_onehot]
  rw [show cRows = transpose 128 cTransposeRows from rfl,
    transpose_bit cTransposeRows i.isLt (by simpa only [cTranspose_length] using j.isLt)]

theorem C_character_pairing (q : Finset (Fin 19)) (x : Finset (Fin 128)) :
    Spin.chi (CtransposeSet q) x = Spin.chi q (Cset x) :=
  Spin.chi_adjoint_of_basis CtransposeSet Cset CtransposeSet_empty Cset_empty
    CtransposeSet_xor Cset_xor C_basis_pairing q x

/-- The actual dual code, as a set of words rather than an assumed spectrum. -/
noncomputable def dual : Finset (Finset (Fin 128)) := univ.image CtransposeSet

theorem dual_card : dual.card = 2 ^ 19 := by
  rw [dual, Finset.card_image_of_injective _ CtransposeSet_injective]
  simp

theorem dual_xor_closed {u v : Finset (Fin 128)} (hu : u ∈ dual) (hv : v ∈ dual) :
    u ∆ v ∈ dual := by
  obtain ⟨q, _, rfl⟩ := Finset.mem_image.mp hu
  obtain ⟨p, _, rfl⟩ := Finset.mem_image.mp hv
  exact Finset.mem_image.mpr ⟨q ∆ p, Finset.mem_univ _, CtransposeSet_xor q p⟩

theorem orthTo_iff_kernel (x : Finset (Fin 128)) : Spin.OrthTo dual x ↔ Cset x = ∅ := by
  constructor
  · intro h
    apply Finset.eq_empty_iff_forall_notMem.mpr
    intro i hi
    have hm : CtransposeSet {i} ∈ dual := Finset.mem_image.mpr ⟨{i}, Finset.mem_univ _, rfl⟩
    have he : Spin.chi (CtransposeSet {i}) x = 1 :=
      (neg_one_pow_eq_one_iff_even (by norm_num : (-1 : ℤ) ≠ 1)).mpr (h _ hm)
    rw [C_character_pairing, Spin.chi_singleton_left, if_pos hi] at he
    norm_num at he
  · intro h v hv
    obtain ⟨q, _, rfl⟩ := Finset.mem_image.mp hv
    apply (neg_one_pow_eq_one_iff_even (by norm_num : (-1 : ℤ) ≠ 1)).mp
    change Spin.chi (CtransposeSet q) x = 1
    rw [C_character_pairing, h]
    simp [Spin.chi]

theorem actual_fiber (j : ℕ) (x : Finset (Fin 128)) :
    Spin.fiber dual j x = (Spin.layer 128 j).filter (fun y => Cset y = Cset x) := by
  ext y
  simp only [Spin.fiber, Finset.mem_filter, orthTo_iff_kernel, Cset_xor, Finset.symmDiff_eq_empty]
  exact and_congr_right fun _ => eq_comm

end Spin.Structured.ConcreteMaps
