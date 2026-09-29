import SpinCodes.Structured.LowCancellationBridge

noncomputable section
namespace Spin.Structured.LowCancellation
open PackedMap ConcreteMaps

def groupMoment (z : ℝ) (g : Group) : ℝ := (g.entries.map fun r => z ^ r.2).sum

theorem cancellation_eq {j : Nat} {gs : List Group}
    (hv : ∀ g ∈ gs, valid j g) (hn : (inputs gs).Nodup)
    (hl : (inputs gs).length = Nat.choose 128 j) (z : ℝ) (q : Finset (Fin 19)) :
    cancellationMoment j z q =
      (gs.map fun g => if support 19 g.syndrome = q then groupMoment z g else 0).sum /
        Nat.choose 128 j := by
  unfold cancellationMoment syndromeFiber
  rw [Finset.sum_filter, layer_sum hv hn hl]
  congr 1
  apply congrArg List.sum
  apply List.map_congr_left
  intro g hg
  have he : (g.entries.map (fun r =>
      if Cset (support 128 r.1) = q then emitted z q (support 128 r.1) else 0)) =
      g.entries.map (fun r => if support 19 g.syndrome = q then z ^ r.2 else 0) := by
    apply List.map_congr_left
    intro r hr
    rw [valid_syndrome (hv g hg) hr]
    split_ifs with h
    · rw [← h, valid_emitted (hv g hg) hr]
    · rfl
  rw [he]
  by_cases h : support 19 g.syndrome = q <;> simp [h, groupMoment]

theorem sum_key_le {α κ : Type*} [DecidableEq κ] (xs : List α) (key : α → κ)
    (f : α → ℝ) (k : κ) (b : ℝ) (hb : 0 ≤ b)
    (hn : (xs.map key).Nodup) (hf : ∀ x ∈ xs, f x ≤ b) :
    (xs.map fun x => if key x = k then f x else 0).sum ≤ b := by
  induction xs with
  | nil => simpa using hb
  | cons a xs ih =>
    simp only [List.map_cons, List.nodup_cons] at hn
    simp only [List.map_cons, List.sum_cons]
    by_cases h : key a = k
    · have hz : (xs.map fun x => if key x = k then f x else 0).sum = 0 := by
        apply List.sum_eq_zero
        intro v hv
        obtain ⟨x, hx, rfl⟩ := List.mem_map.mp hv
        have hne : key x ≠ k := by
          intro he
          apply hn.1
          exact List.mem_map.mpr ⟨x, hx, he.trans h.symm⟩
        simp [hne]
      simpa [h, hz] using hf a (by simp)
    · simp only [if_neg h, zero_add]
      exact ih hn.2 (fun x hx => hf x (by simp [hx]))

theorem expand_sum (p : List (Nat × Nat)) (z : ℝ) :
    ((expand p).map (fun e => z ^ e)).sum =
      (p.map fun ec => (ec.2 : ℝ) * z ^ ec.1).sum := by
  induction p with
  | nil => rfl
  | cons ec p ih =>
    change ((List.replicate ec.2 ec.1 ++ expand p).map (fun e => z ^ e)).sum = _
    simp [List.map_append, List.sum_append, ih, nsmul_eq_mul]

theorem list_sum_div (xs : List ℝ) (d : ℝ) : xs.sum / d = (xs.map (· / d)).sum := by
  induction xs with
  | nil => simp
  | cons a xs ih => simp [add_div, ih]

theorem le_foldr_max {x : ℝ} {xs : List ℝ} (h : x ∈ xs) : x ≤ xs.foldr max 0 := by
  induction xs with
  | nil => simp at h
  | cons y ys ih =>
    rcases List.mem_cons.mp h with rfl | h
    · exact le_max_left _ _
    · exact (ih h).trans (le_max_right _ _)

theorem foldr_max_nonneg (xs : List ℝ) : 0 ≤ xs.foldr max 0 := by
  induction xs with
  | nil => exact le_rfl
  | cons x xs ih => exact ih.trans (le_max_right _ _)

theorem cancellation_le_patterns {j : Nat} {gs : List Group} {ps : List (List (Nat × Nat))}
    (hv : ∀ g ∈ gs, valid j g) (hn : (inputs gs).Nodup)
    (hl : (inputs gs).length = Nat.choose 128 j)
    (hk : (gs.map Group.syndrome).Nodup)
    (hp : ∀ g ∈ gs, g.entries.map Prod.snd ∈ ps.map expand) (z : ℝ) (q) :
    cancellationMoment j z q ≤
      (ps.map fun p => Spin.Imt.Occupation.Sparse.pattern j p z).foldr max 0 := by
  rw [cancellation_eq hv hn hl]
  rw [list_sum_div]
  have he : (gs.map (fun g => if support 19 g.syndrome = q then groupMoment z g else 0)).map
      (fun a => a / (Nat.choose 128 j : ℝ)) =
      gs.map (fun g => if support 19 g.syndrome = q then groupMoment z g / Nat.choose 128 j else 0) := by
    simp only [List.map_map]
    apply List.map_congr_left
    intro g _
    by_cases h : support 19 g.syndrome = q <;> simp [h]
  rw [he]
  apply sum_key_le gs (fun g : Group => support 19 g.syndrome)
    (fun g => groupMoment z g / (Nat.choose 128 j : ℝ)) q _ (foldr_max_nonneg _)
  · change (gs.map ((support 19) ∘ Group.syndrome)).Nodup
    rw [← List.map_map]
    apply (List.nodup_map_iff_inj_on hk).mpr
    intro a ha b hb hab
    obtain ⟨g, hg, rfl⟩ := List.mem_map.mp ha
    obtain ⟨g', hg', rfl⟩ := List.mem_map.mp hb
    exact support_inj_bounded (hv g hg).1 (hv g' hg').1 hab
  · intro g hg
    obtain ⟨p, hp', he⟩ := List.mem_map.mp (hp g hg)
    have hm : groupMoment z g / Nat.choose 128 j = Spin.Imt.Occupation.Sparse.pattern j p z := by
      unfold groupMoment Spin.Imt.Occupation.Sparse.pattern
      rw [← expand_sum, he, List.map_map]
      rfl
    rw [hm]
    exact le_foldr_max (List.mem_map.mpr ⟨p, hp', rfl⟩)

end Spin.Structured.LowCancellation
