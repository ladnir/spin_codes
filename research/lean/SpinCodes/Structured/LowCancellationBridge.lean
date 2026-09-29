import SpinCodes.Structured.LowCancellationDefs
import SpinCodes.Structured.ConcreteCancellation
import SpinCodes.Structured.LowCancellationSort

noncomputable section
namespace Spin.Structured.LowCancellation
open PackedMap ConcreteMaps
open scoped symmDiff

theorem sorted_nodup (xs : List Nat)
    (h : (certSort xs).IsChain (· < ·)) : xs.Nodup := by
  have hn : (certSort xs).Nodup :=
    (List.isChain_iff_pairwise.mp h).imp (fun h => Nat.ne_of_lt h)
  exact (certSort_perm xs).nodup_iff.mp hn

theorem support_inj_bounded {n x y : Nat} (hx : x < 2 ^ n) (hy : y < 2 ^ n)
    (h : support n x = support n y) : x = y :=
  congrArg Fin.val (support_injective n (a₁ := ⟨x, hx⟩) (a₂ := ⟨y, hy⟩) h)

theorem valid_input {j : Nat} {gs : List Group} (hv : ∀ g ∈ gs, valid j g)
    {x : Nat} (hx : x ∈ inputs gs) : x < 2 ^ 128 ∧ weight 128 x = j := by
  obtain ⟨g, hg, hx⟩ := List.mem_flatMap.mp hx
  obtain ⟨r, hr, rfl⟩ := List.mem_map.mp hx
  exact ⟨(hv g hg).2.2.2 r hr |>.1, (hv g hg).2.2.2 r hr |>.2.1⟩

theorem inputs_support_nodup {j : Nat} {gs : List Group}
    (hv : ∀ g ∈ gs, valid j g) (hn : (inputs gs).Nodup) :
    ((inputs gs).map (support 128)).Nodup := by
  apply (List.nodup_map_iff_inj_on hn).mpr
  intro x hx y hy h
  exact support_inj_bounded (valid_input hv hx).1 (valid_input hv hy).1 h

theorem inputs_layer {j : Nat} {gs : List Group}
    (hv : ∀ g ∈ gs, valid j g) (hn : (inputs gs).Nodup)
    (hl : (inputs gs).length = Nat.choose 128 j) :
    ((inputs gs).map (support 128)).toFinset = Spin.layer 128 j := by
  apply Finset.eq_of_subset_of_card_le
  · intro x hx
    obtain ⟨a, ha, rfl⟩ := List.mem_map.mp (List.mem_toFinset.mp hx)
    apply Finset.mem_powersetCard.mpr
    exact ⟨Finset.subset_univ _, (weight_eq_card_support 128 a).symm.trans (valid_input hv ha).2⟩
  · rw [List.toFinset_card_of_nodup (inputs_support_nodup hv hn), List.length_map,
      hl, Spin.card_layer]

theorem layer_sum {j : Nat} {gs : List Group}
    (hv : ∀ g ∈ gs, valid j g) (hn : (inputs gs).Nodup)
    (hl : (inputs gs).length = Nat.choose 128 j) (f : Finset (Fin 128) → ℝ) :
    (∑ x ∈ Spin.layer 128 j, f x) =
      (gs.map fun g => (g.entries.map fun r => f (support 128 r.1)).sum).sum := by
  rw [← inputs_layer hv hn hl, List.sum_toFinset _ (inputs_support_nodup hv hn)]
  simp only [inputs, List.map_flatMap, List.map_map, Function.comp_def]
  clear hv hn hl j
  induction gs with
  | nil => rfl
  | cons g gs ih =>
    simp only [List.flatMap_cons, List.sum_append, List.map_cons, List.sum_cons]
    congr 1

theorem Cset_support {x : Nat} (hx : x < 2 ^ 128) :
    Cset (support 128 x) = support 19 (eval cRows x) := by
  change supportEquiv 19 (C ((supportEquiv 128).symm (supportEquiv 128 ⟨x, hx⟩))) = _
  rw [Equiv.symm_apply_apply]
  rfl

theorem Aset_support {q : Nat} (hq : q < 2 ^ 19) :
    Aset (support 19 q) = support 128 (eval aRows q) := by
  change supportEquiv 128 (A ((supportEquiv 19).symm (supportEquiv 19 ⟨q, hq⟩))) = _
  rw [Equiv.symm_apply_apply]
  rfl

theorem valid_emitted {j : Nat} {g : Group} (hv : valid j g)
    {r : Nat × Nat} (hr : r ∈ g.entries) (z : ℝ) :
    emitted z (support 19 g.syndrome) (support 128 r.1) = z ^ r.2 := by
  rw [emitted, Aset_support hv.1, ← support_xor, ← weight_eq_card_support,
    (hv.2.2.2 r hr).2.2.2]

theorem valid_syndrome {j : Nat} {g : Group} (hv : valid j g)
    {r : Nat × Nat} (hr : r ∈ g.entries) :
    Cset (support 128 r.1) = support 19 g.syndrome := by
  rw [Cset_support (hv.2.2.2 r hr).1, (hv.2.2.2 r hr).2.2.1]

end Spin.Structured.LowCancellation
