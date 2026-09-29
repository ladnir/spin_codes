/-
The fourth integer fibre bound: the variance bound.

`app:imt-finite-transfers` bounds a nonzero fibre by

    ⌊(a_j + ⌈√((M-1)(M(V_j - k_j²) - a_j²))⌉) / M⌋,

obtained "by minimizing the other M-1 squares for a fixed maximum fibre".
That minimisation is Cauchy-Schwarz, and solving the resulting quadratic for
the fibre size is what produces the square root.  Here the quadratic itself is
proved,

    M·t² + a_j² ≤ 2·a_j·t + (M-1)·(V_j - k_j²),

which is the same statement without the square root, and therefore without
any rounding.  The appendix's displayed formula is this inequality solved for
`t` and rounded outward.

Unlike the other three bounds, this one needs the *whole family* of fibre
sizes at once, so most of the file is the partition of the non-kernel layer
into fibres.
-/
import SpinCodes.Structured.FiberBounds

namespace Spin

open Finset
open scoped symmDiff

/-! ## Cauchy-Schwarz on all but one term

Stated for an arbitrary finite family, entirely in `ℕ`: writing the totals as
`a = A + t` and `S = B + t²` removes every subtraction, and what is left is
exactly `A² ≤ m·B`. -/

theorem sq_max_le {ι : Type*} [DecidableEq ι] (F : Finset ι) (g : ι → ℕ)
    {i₀ : ι} (hi₀ : i₀ ∈ F) {m a S : ℕ}
    (hm : F.card ≤ m + 1) (ha : ∑ i ∈ F, g i = a) (hS : ∑ i ∈ F, (g i) ^ 2 = S) :
    (m + 1) * (g i₀) ^ 2 + a ^ 2 ≤ 2 * a * g i₀ + m * S := by
  classical
  set A := ∑ i ∈ F.erase i₀, g i with hA
  set B := ∑ i ∈ F.erase i₀, (g i) ^ 2 with hB
  have haA : a = A + g i₀ := by
    rw [← ha, hA, Finset.sum_erase_add F g hi₀]
  have hSB : S = B + (g i₀) ^ 2 := by
    rw [← hS, hB, Finset.sum_erase_add F (fun i => (g i) ^ 2) hi₀]
  have hcs : A ^ 2 ≤ (F.erase i₀).card * B := sq_sum_le_card_mul_sum_sq
  have hcard : (F.erase i₀).card ≤ m := by
    rw [Finset.card_erase_of_mem hi₀]
    omega
  have hkey : A ^ 2 ≤ m * B := le_trans hcs (Nat.mul_le_mul_right _ hcard)
  subst haA
  subst hSB
  nlinarith [hkey]

/-! ## Symmetric-difference algebra

`(Finset α, ∆)` is a group, so these are one-line consequences of
associativity and `a ∆ a = ⊥`.  Proving them by `ext; simp; tauto` instead
makes the elaborator time out. -/

lemma symmDiff_cancel_mid {n : ℕ} (a b c : Finset (Fin n)) :
    (a ∆ c) ∆ (b ∆ c) = a ∆ b := by
  rw [symmDiff_symmDiff_symmDiff_comm, symmDiff_self, symmDiff_bot]

lemma symmDiff_cancel_outer {n : ℕ} (a b c : Finset (Fin n)) :
    (a ∆ b) ∆ (a ∆ c) = b ∆ c := by
  rw [symmDiff_symmDiff_symmDiff_comm, symmDiff_self, bot_symmDiff]

lemma symmDiff_cancel_tail {n : ℕ} (a b : Finset (Fin n)) : (a ∆ b) ∆ b = a := by
  rw [symmDiff_assoc, symmDiff_self, symmDiff_bot]

/-! ## The kernel, and fibres as cosets -/

lemma orthTo_empty {n : ℕ} (D : Finset (Finset (Fin n))) : OrthTo D (∅ : Finset (Fin n)) := by
  intro v _
  simp

lemma self_mem_fiber {n : ℕ} (D : Finset (Finset (Fin n))) {j : ℕ}
    {y : Finset (Fin n)} (hy : y ∈ layer n j) : y ∈ fiber D j y := by
  refine Finset.mem_filter.mpr ⟨hy, ?_⟩
  have hyy : y ∆ y = (∅ : Finset (Fin n)) := by
    rw [symmDiff_self]; rfl
  rw [hyy]
  exact orthTo_empty D

lemma mem_fiber_iff {n : ℕ} (D : Finset (Finset (Fin n))) (j : ℕ)
    (x y : Finset (Fin n)) : y ∈ fiber D j x ↔ y ∈ layer n j ∧ OrthTo D (x ∆ y) :=
  Finset.mem_filter

/-- Fibres through words differing by a kernel element coincide. -/
lemma fiber_eq_of_orthTo {n : ℕ} (D : Finset (Finset (Fin n))) (j : ℕ)
    {x x' : Finset (Fin n)} (h : OrthTo D (x ∆ x')) :
    fiber D j x = fiber D j x' := by
  ext y
  rw [mem_fiber_iff, mem_fiber_iff]
  constructor
  · rintro ⟨hy, hxy⟩
    refine ⟨hy, ?_⟩
    have hrw : (x ∆ x') ∆ (x ∆ y) = x' ∆ y := symmDiff_cancel_outer x x' y
    have hres := orthTo_symmDiff h hxy
    rwa [hrw] at hres
  · rintro ⟨hy, hxy⟩
    refine ⟨hy, ?_⟩
    have hrw : (x ∆ x') ∆ (x' ∆ y) = x ∆ y := by
      rw [symmDiff_assoc, symmDiff_symmDiff_cancel_left]
    have hres := orthTo_symmDiff h hxy
    rwa [hrw] at hres

/-- Distinct fibres are disjoint. -/
lemma fiber_disjoint {n : ℕ} (D : Finset (Finset (Fin n))) (j : ℕ)
    {x x' : Finset (Fin n)} (h : fiber D j x ≠ fiber D j x') :
    Disjoint (fiber D j x) (fiber D j x') := by
  refine Finset.disjoint_left.mpr fun z hz hz' => ?_
  have h1 : OrthTo D (x ∆ z) := (mem_fiber_iff D j x z).mp hz |>.2
  have h2 : OrthTo D (x' ∆ z) := (mem_fiber_iff D j x' z).mp hz' |>.2
  have hrw : (x ∆ z) ∆ (x' ∆ z) = x ∆ x' := symmDiff_cancel_mid x x' z
  have hres := orthTo_symmDiff h1 h2
  rw [hrw] at hres
  exact h (fiber_eq_of_orthTo D j hres)

/-! ## The partition of the non-kernel layer -/

/-- The weight-`j` words with nonzero syndrome. -/
def nonKernelLayer {n : ℕ} (D : Finset (Finset (Fin n))) (j : ℕ) :
    Finset (Finset (Fin n)) :=
  (layer n j).filter (fun y => ¬ OrthTo D y)

/-- The weight-`j` kernel words. -/
def kernelLayer {n : ℕ} (D : Finset (Finset (Fin n))) (j : ℕ) :
    Finset (Finset (Fin n)) :=
  (layer n j).filter (fun y => OrthTo D y)

/-- The distinct nonzero fibres. -/
def fibers {n : ℕ} (D : Finset (Finset (Fin n))) (j : ℕ) :
    Finset (Finset (Finset (Fin n))) :=
  (nonKernelLayer D j).image (fun y => fiber D j y)

lemma fiber_subset_nonKernel {n : ℕ} (D : Finset (Finset (Fin n))) (j : ℕ)
    {x : Finset (Fin n)} (hx : ¬ OrthTo D x) :
    fiber D j x ⊆ nonKernelLayer D j := by
  intro y hy
  rw [mem_fiber_iff] at hy
  refine Finset.mem_filter.mpr ⟨hy.1, fun hyk => hx ?_⟩
  have hrw : (x ∆ y) ∆ y = x := symmDiff_cancel_tail x y
  have hres := orthTo_symmDiff hy.2 hyk
  rwa [hrw] at hres

lemma pairwiseDisjoint_fibers {n : ℕ} (D : Finset (Finset (Fin n))) (j : ℕ) :
    (↑(fibers D j) : Set (Finset (Finset (Fin n)))).PairwiseDisjoint (fun Φ => Φ) := by
  intro Φ hΦ Φ' hΦ' hne
  simp only [fibers, Finset.coe_image, Set.mem_image, Finset.mem_coe] at hΦ hΦ'
  obtain ⟨x, _, rfl⟩ := hΦ
  obtain ⟨x', _, rfl⟩ := hΦ'
  exact fiber_disjoint D j hne

lemma biUnion_fibers {n : ℕ} (D : Finset (Finset (Fin n))) (j : ℕ) :
    (fibers D j).biUnion (fun Φ => Φ) = nonKernelLayer D j := by
  ext y
  rw [Finset.mem_biUnion]
  constructor
  · rintro ⟨Φ, hΦ, hyΦ⟩
    obtain ⟨x, hx, rfl⟩ := Finset.mem_image.mp hΦ
    exact fiber_subset_nonKernel D j (Finset.mem_filter.mp hx).2 hyΦ
  · intro hy
    exact ⟨fiber D j y, Finset.mem_image_of_mem _ hy,
      self_mem_fiber D (Finset.mem_filter.mp hy).1⟩

/-- Every word in a fibre names that fibre. -/
lemma fiber_of_mem {n : ℕ} (D : Finset (Finset (Fin n))) (j : ℕ)
    {x y : Finset (Fin n)} (hy : y ∈ fiber D j x) : fiber D j y = fiber D j x := by
  refine fiber_eq_of_orthTo D j ?_
  have hxy : OrthTo D (x ∆ y) := (mem_fiber_iff D j x y |>.mp hy).2
  have hrw : y ∆ x = x ∆ y := symmDiff_comm y x
  rwa [hrw]

/-- `a_j` is the sum of the fibre sizes. -/
theorem sum_card_fibers {n : ℕ} (D : Finset (Finset (Fin n))) (j : ℕ) :
    ∑ Φ ∈ fibers D j, Φ.card = (nonKernelLayer D j).card := by
  rw [← biUnion_fibers D j, Finset.card_biUnion]
  intro Φ hΦ Φ' hΦ' hne
  exact pairwiseDisjoint_fibers D j (Finset.mem_coe.mpr hΦ) (Finset.mem_coe.mpr hΦ') hne

/-- `V_j - k_j²` is the sum of squares of the fibre sizes — stated additively. -/
theorem card_orth_pairs_split {n : ℕ} (D : Finset (Finset (Fin n))) (j : ℕ) :
    (((layer n j) ×ˢ (layer n j)).filter (fun p => OrthTo D (p.1 ∆ p.2))).card
      = (kernelLayer D j).card ^ 2 + ∑ Φ ∈ fibers D j, Φ.card ^ 2 := by
  classical
  -- count the product set row by row
  have hrow : (((layer n j) ×ˢ (layer n j)).filter
      (fun p => OrthTo D (p.1 ∆ p.2))).card
      = ∑ x ∈ layer n j, (fiber D j x).card := by
    rw [Finset.card_filter, Finset.sum_product]
    exact Finset.sum_congr rfl fun x _ => by
      rw [fiber, Finset.card_filter]
  rw [hrow, ← Finset.sum_filter_add_sum_filter_not (layer n j) (fun y => OrthTo D y)]
  congr 1
  · -- kernel rows: each fibre is the whole kernel layer
    have hk : ∀ x ∈ (layer n j).filter (fun y => OrthTo D y),
        (fiber D j x).card = ((layer n j).filter (fun y => OrthTo D y)).card := by
      intro x hx
      have hxk : OrthTo D x := (Finset.mem_filter.mp hx).2
      congr 1
      ext y
      rw [mem_fiber_iff]
      simp only [Finset.mem_filter]
      constructor
      · rintro ⟨hy, hxy⟩
        refine ⟨hy, ?_⟩
        have hrw : x ∆ (x ∆ y) = y := symmDiff_symmDiff_cancel_left x y
        have hres := orthTo_symmDiff hxk hxy
        rwa [hrw] at hres
      · rintro ⟨hy, hyk⟩
        exact ⟨hy, orthTo_symmDiff hxk hyk⟩
    rw [Finset.sum_congr rfl hk, Finset.sum_const, smul_eq_mul]
    show _ = ((layer n j).filter (fun y => OrthTo D y)).card ^ 2
    rw [sq]
  · -- non-kernel rows: group them into fibres
    show ∑ x ∈ nonKernelLayer D j, (fiber D j x).card = _
    rw [← biUnion_fibers D j, Finset.sum_biUnion (pairwiseDisjoint_fibers D j)]
    refine Finset.sum_congr rfl fun Φ hΦ => ?_
    obtain ⟨x, _, rfl⟩ := Finset.mem_image.mp hΦ
    rw [Finset.sum_congr rfl (fun y hy => by rw [fiber_of_mem D j hy]),
      Finset.sum_const, smul_eq_mul, sq]

/-! ## The bound -/

/-- **The variance bound.**  `m + 1` bounds the number of nonzero fibres and
`t` is the size of any one of them.  This is the appendix's bound before it is
solved for `t` and rounded outward; the sum of squares is `V_j - k_j²` by
`card_orth_pairs_split`. -/
theorem variance_bound {n : ℕ} (D : Finset (Finset (Fin n))) (j : ℕ)
    {x : Finset (Fin n)} (hxl : x ∈ layer n j) (hx : ¬ OrthTo D x) {m : ℕ}
    (hm : (fibers D j).card ≤ m + 1) :
    (m + 1) * (fiber D j x).card ^ 2 + (nonKernelLayer D j).card ^ 2
      ≤ 2 * (nonKernelLayer D j).card * (fiber D j x).card
        + m * ∑ Φ ∈ fibers D j, Φ.card ^ 2 := by
  classical
  have hmem : fiber D j x ∈ fibers D j :=
    Finset.mem_image_of_mem _ (Finset.mem_filter.mpr ⟨hxl, hx⟩)
  exact sq_max_le (fibers D j) (fun Φ => Φ.card) hmem hm (sum_card_fibers D j) rfl

/-- The same bound written with the pair count `V_j` and the kernel count
`k_j`, additively so that no truncated subtraction appears. -/
theorem variance_bound_pairs {n : ℕ} (D : Finset (Finset (Fin n))) (j : ℕ)
    {x : Finset (Fin n)} (hxl : x ∈ layer n j) (hx : ¬ OrthTo D x) {m : ℕ}
    (hm : (fibers D j).card ≤ m + 1) :
    (m + 1) * (fiber D j x).card ^ 2 + (nonKernelLayer D j).card ^ 2
        + m * (kernelLayer D j).card ^ 2
      ≤ 2 * (nonKernelLayer D j).card * (fiber D j x).card
        + m * (((layer n j) ×ˢ (layer n j)).filter
              (fun p => OrthTo D (p.1 ∆ p.2))).card := by
  have hsplit := card_orth_pairs_split D j
  have hb := variance_bound D j hxl hx hm
  rw [hsplit]
  nlinarith [hb]

end Spin
