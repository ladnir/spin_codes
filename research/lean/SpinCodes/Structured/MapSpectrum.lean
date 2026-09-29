import SpinCodes.Structured.MapSpectrumDefs
import SpinCodes.Structured.PackedMap

/-! A kernel-checked histogram counts the actual map's input vectors. Splitting
the input range into blocks preserves every input exactly once. -/

namespace Spin.Structured.MapSpectrum

theorem length_bump (w : ℕ) (h : List ℕ) : (bump w h).length = h.length := by
  induction h generalizing w with
  | nil => cases w <;> rfl
  | cons a h ih => cases w <;> simp [bump, ih]

theorem getD_bump (w i : ℕ) (h : List ℕ) (hi : i < h.length) :
    (bump w h).getD i 0 = h.getD i 0 + if w = i then 1 else 0 := by
  induction h generalizing w i with
  | nil => simp at hi
  | cons a h ih =>
    cases w with
    | zero => cases i <;> simp [bump]
    | succ w =>
      cases i with
      | zero => simp [bump]
      | succ i => simpa [bump] using ih w i (by simpa using hi)

theorem foldl_bump_count {α : Type*} (f : α → ℕ) (xs : List α) (h : List ℕ)
    (i : ℕ) (hi : i < h.length) :
    (xs.foldl (fun h x => bump (f x) h) h).getD i 0 =
      h.getD i 0 + xs.countP (fun x => f x == i) := by
  induction xs generalizing h with
  | nil => simp
  | cons x xs ih =>
    simp only [List.foldl_cons, List.countP_cons]
    rw [ih _ (by simpa only [length_bump] using hi), getD_bump _ _ _ hi]
    by_cases h : f x = i <;> simp [h] <;> omega

theorem histogram_count (rows inputs : List ℕ) {i : ℕ} (hi : i < 129) :
    (histogram rows inputs).getD i 0 =
      inputs.countP (fun q => PackedMap.weight 128 (PackedMap.eval rows q) == i) := by
  rw [histogram, foldl_bump_count _ _ _ _ (by simpa using hi)]
  simp only [List.getD_replicate 0 hi, Nat.zero_add]

theorem countP_range_eq_card (p : ℕ → Bool) (n : ℕ) :
    (List.range n).countP p = (Finset.univ.filter fun q : Fin n => p q).card := by
  rw [Finset.card_filter]
  change (List.range n).countP p = ∑ q : Fin n, (fun i : ℕ => if p i then 1 else 0) q
  rw [Fin.sum_univ_eq_sum_range (fun i : ℕ => if p i then (1 : ℕ) else 0) n]
  induction n with
  | zero => simp
  | succ n ih =>
    simp only [List.range_succ, List.countP_append, List.countP_cons, List.countP_nil,
      Nat.zero_add, Finset.sum_range_succ, ih]

theorem histogram_range_card (rows : List ℕ) (n : ℕ) {i : ℕ} (hi : i < 129) :
    (histogram rows (List.range n)).getD i 0 =
      (Finset.univ.filter fun q : Fin n => PackedMap.weight 128 (PackedMap.eval rows q) = i).card := by
  rw [histogram_count _ _ hi, countP_range_eq_card]
  simp

theorem range_blocks (size count : ℕ) :
    (List.range count).flatMap (fun b => List.range' (b * size) size) =
      List.range' 0 (count * size) := by
  induction count with
  | zero => simp
  | succ count ih =>
    rw [List.range_succ, List.flatMap_append, ih]
    simp only [List.flatMap_cons, List.flatMap_nil, List.append_nil]
    simpa only [Nat.zero_add, Nat.succ_mul] using
      (List.range'_append_1 (s := 0) (m := count * size) (n := size))

theorem histogram_blocks (rows : List ℕ) (size count : ℕ) {i : ℕ} (hi : i < 129) :
    (histogram rows (List.range (count * size))).getD i 0 =
      ((List.range count).map fun b => (block rows (b * size) size).getD i 0).sum := by
  rw [histogram_count _ _ hi, List.range_eq_range', ← range_blocks, List.countP_flatMap]
  apply congrArg List.sum
  apply List.map_congr_left
  intro b _
  exact (histogram_count _ _ hi).symm

theorem block_split_weights (rows : List ℕ) (k high : ℕ) (hk : k ≤ rows.length)
    (lowImages : List ℕ) (highImage : ℕ) (ws : List ℕ)
    (hlow : (List.range (2 ^ k)).map (PackedMap.eval (rows.take k)) = lowImages)
    (hhigh : PackedMap.eval (rows.drop k) high = highImage)
    (hw : lowImages.map (fun q => PackedMap.weight 128 (q ^^^ highImage)) = ws) :
    block rows (high * 2 ^ k) (2 ^ k) =
      ws.foldl (fun h w => bump w h) (List.replicate 129 0) := by
  apply block_of_weights
  rw [← hw, ← hlow, List.map_map, List.range'_eq_map_range, List.map_map]
  apply List.map_congr_left
  intro low hl
  simp only [Function.comp_def]
  rw [PackedMap.eval_split rows hk (List.mem_range.mp hl) high, hhigh]

end Spin.Structured.MapSpectrum
