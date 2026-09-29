/-
Integer upper bounds on a nonzero syndrome fibre.

`app:imt-finite-transfers` takes the minimum of four bounds on `f_j`, the size
of any nonzero weight-`j` syndrome fibre.  Three of them are proved here:

*  `a_j = n_j - k_j`, because a fibre with nonzero syndrome misses the kernel;
*  `⌊|D|⁻¹ ∑_w b_w |K_j(w)|⌋`, the triangle inequality on Fourier inversion;
*  the constant-weight packing bound, from the kernel distance.

The fourth (the variance bound, from Cauchy-Schwarz across all `M` nonzero
syndromes) needs the whole family of fibre sizes at once and is separate.

As in `Parseval.lean`, nothing here mentions a syndrome map: a fibre is the
set of weight-`j` words whose symmetric difference with `x` lies in the kernel.
-/
import SpinCodes.Structured.Parseval

namespace Spin

open Finset
open scoped symmDiff

/-! ## Two closure facts -/

/-- The kernel is closed under symmetric difference. -/
lemma orthTo_symmDiff {n : ℕ} {D : Finset (Finset (Fin n))} {a b : Finset (Fin n)}
    (ha : OrthTo D a) (hb : OrthTo D b) : OrthTo D (a ∆ b) := by
  intro v hv
  obtain ⟨p, hp⟩ := ha v hv
  obtain ⟨q, hq⟩ := hb v hv
  have hd : v ∩ (a ∆ b) = (v ∩ a) ∆ (v ∩ b) := inter_symmDiff_right v a b
  have hc := card_symmDiff_add_two_mul_card_inter (v ∩ a) (v ∩ b)
  rw [hd]
  exact ⟨p + q - ((v ∩ a) ∩ (v ∩ b)).card, by omega⟩

/-- Two members of one fibre differ by a kernel word. -/
lemma symmDiff_mem_kernel {n : ℕ} {D : Finset (Finset (Fin n))} {j : ℕ}
    {x y y' : Finset (Fin n)} (hy : y ∈ fiber D j x) (hy' : y' ∈ fiber D j x) :
    OrthTo D (y ∆ y') := by
  have h1 : OrthTo D (x ∆ y) := (Finset.mem_filter.mp hy).2
  have h2 : OrthTo D (x ∆ y') := (Finset.mem_filter.mp hy').2
  have hcancel : (x ∆ y) ∆ (x ∆ y') = y ∆ y' := by
    ext a; simp [Finset.mem_symmDiff]; tauto
  rw [← hcancel]
  exact orthTo_symmDiff h1 h2

lemma card_layer (n j : ℕ) : (layer n j).card = n.choose j := by
  rw [layer, Finset.card_powersetCard, Finset.card_univ, Fintype.card_fin]

lemma fiber_subset_layer {n : ℕ} (D : Finset (Finset (Fin n))) (j : ℕ)
    (x : Finset (Fin n)) : fiber D j x ⊆ layer n j := Finset.filter_subset _ _

lemma card_eq_of_mem_fiber {n : ℕ} {D : Finset (Finset (Fin n))} {j : ℕ}
    {x y : Finset (Fin n)} (hy : y ∈ fiber D j x) : y.card = j :=
  (Finset.mem_powersetCard.mp (fiber_subset_layer D j x hy)).2

/-! ## Bound 1: a nonzero fibre misses the kernel -/

/-- **`f_j ≤ a_j`.**  Stated additively to avoid truncated subtraction: the
fibre and the kernel layer are disjoint subsets of the weight-`j` layer. -/
theorem card_fiber_add_card_kernel_le {n : ℕ} (D : Finset (Finset (Fin n))) (j : ℕ)
    {x : Finset (Fin n)} (hx : ¬ OrthTo D x) :
    (fiber D j x).card + ((layer n j).filter (fun y => OrthTo D y)).card
      ≤ n.choose j := by
  classical
  have hdisj : Disjoint (fiber D j x) ((layer n j).filter (fun y => OrthTo D y)) := by
    refine Finset.disjoint_left.mpr fun y hy hy' => ?_
    have h1 : OrthTo D (x ∆ y) := (Finset.mem_filter.mp hy).2
    have h2 : OrthTo D y := (Finset.mem_filter.mp hy').2
    have hcancel : (x ∆ y) ∆ y = x := by
      ext a; simp
    exact hx (hcancel ▸ orthTo_symmDiff h1 h2)
  have hsub : fiber D j x ∪ ((layer n j).filter (fun y => OrthTo D y)) ⊆ layer n j :=
    Finset.union_subset (fiber_subset_layer D j x) (Finset.filter_subset _ _)
  have hle := Finset.card_le_card hsub
  rwa [Finset.card_union_of_disjoint hdisj, card_layer] at hle

/-! ## Bound 2: the triangle inequality on Fourier inversion -/

lemma abs_chi {n : ℕ} (u x : Finset (Fin n)) : |chi u x| = 1 := by
  unfold chi
  rw [abs_pow, abs_neg, abs_one, one_pow]

/-- **`|D| · f_j ≤ ∑_{v ∈ D} |K_j(|v|)|`.**  Splitting at `v = ∅` gives the
appendix form `n_j + ∑_w b_w |K_j(w)|`. -/
theorem card_fiber_le_abs_sum {n : ℕ} (D : Finset (Finset (Fin n)))
    (hD : ∀ u ∈ D, ∀ v ∈ D, u ∆ v ∈ D) (j : ℕ) (x : Finset (Fin n)) :
    (D.card : ℤ) * ((fiber D j x).card : ℤ)
      ≤ ∑ v ∈ D, |krawtchouk n j v.card| := by
  rw [card_fiber D hD j x]
  refine Finset.sum_le_sum fun v _ => ?_
  calc chi v x * krawtchouk n j v.card
      ≤ |chi v x * krawtchouk n j v.card| := le_abs_self _
    _ = |krawtchouk n j v.card| := by rw [abs_mul, abs_chi, one_mul]

/-! ## Bound 3: constant-weight packing

Stated for an arbitrary constant-weight set with a pairwise distance
guarantee, so that it applies both to a fibre and to its complement. -/

/-- If distinct words of weight `h` are pairwise at distance more than `2r`,
then no `(h-r)`-subset lies in two of them, and counting `(h-r)`-subsets
bounds the number of words. -/
theorem card_mul_choose_le_of_dist {n : ℕ} (W : Finset (Finset (Fin n)))
    {h r d : ℕ} (hw : ∀ y ∈ W, y.card = h) (hr : 2 * r < d)
    (hdist : ∀ y ∈ W, ∀ y' ∈ W, y ≠ y' → d ≤ (y ∆ y').card) :
    W.card * h.choose (h - r) ≤ n.choose (h - r) := by
  classical
  have hdisj : ∀ y ∈ W, ∀ y' ∈ W, y ≠ y' →
      Disjoint (y.powersetCard (h - r)) (y'.powersetCard (h - r)) := by
    intro y hy y' hy' hne
    refine Finset.disjoint_left.mpr fun T hT hT' => ?_
    rw [Finset.mem_powersetCard] at hT hT'
    have hTsub : T ⊆ y ∩ y' := Finset.subset_inter hT.1 hT'.1
    have hTle : h - r ≤ (y ∩ y').card := hT.2 ▸ Finset.card_le_card hTsub
    have hcard := card_symmDiff_add_two_mul_card_inter y y'
    rw [hw y hy, hw y' hy'] at hcard
    have hd := hdist y hy y' hy' hne
    omega
  have hbi : (W.biUnion (fun y => y.powersetCard (h - r))).card
      = W.card * h.choose (h - r) := by
    rw [Finset.card_biUnion hdisj,
      Finset.sum_congr rfl (fun y hy => by
        rw [Finset.card_powersetCard, hw y hy]),
      Finset.sum_const, smul_eq_mul]
  have hsub : W.biUnion (fun y => y.powersetCard (h - r))
      ⊆ (univ : Finset (Fin n)).powersetCard (h - r) := by
    intro T hT
    rw [Finset.mem_biUnion] at hT
    obtain ⟨y, _, hTy⟩ := hT
    rw [Finset.mem_powersetCard] at hTy ⊢
    exact ⟨Finset.subset_univ _, hTy.2⟩
  have hle := Finset.card_le_card hsub
  rwa [hbi, Finset.card_powersetCard, Finset.card_univ, Fintype.card_fin] at hle

/-- **The packing bound on a fibre.**  `hdist` is the kernel distance. -/
theorem card_fiber_mul_choose_le {n : ℕ} (D : Finset (Finset (Fin n))) (j : ℕ)
    {r d : ℕ} (hr : 2 * r < d)
    (hdist : ∀ u : Finset (Fin n), OrthTo D u → u ≠ ∅ → d ≤ u.card)
    (x : Finset (Fin n)) :
    (fiber D j x).card * j.choose (j - r) ≤ n.choose (j - r) := by
  refine card_mul_choose_le_of_dist _ (fun y hy => card_eq_of_mem_fiber hy) hr ?_
  intro y hy y' hy' hne
  refine hdist _ (symmDiff_mem_kernel hy hy') ?_
  simpa [symmDiff_eq_bot, Finset.bot_eq_empty] using hne

/-- The same bound applied to the complements, which is the stronger one when
`j > n/2` -- the appendix step of complementing when needed. -/
theorem card_fiber_mul_choose_le_compl {n : ℕ} (D : Finset (Finset (Fin n))) (j : ℕ)
    {r d : ℕ} (hr : 2 * r < d)
    (hdist : ∀ u : Finset (Fin n), OrthTo D u → u ≠ ∅ → d ≤ u.card)
    (x : Finset (Fin n)) :
    (fiber D j x).card * (n - j).choose ((n - j) - r) ≤ n.choose ((n - j) - r) := by
  classical
  have hinj : Set.InjOn (fun y : Finset (Fin n) => yᶜ) (fiber D j x) := by
    intro a _ b _ hab
    simpa using congrArg (fun s : Finset (Fin n) => sᶜ) hab
  have hcard : ((fiber D j x).image (fun y => yᶜ)).card = (fiber D j x).card :=
    Finset.card_image_of_injOn hinj
  rw [← hcard]
  refine card_mul_choose_le_of_dist _ ?_ hr ?_
  · intro y hy
    obtain ⟨z, hz, rfl⟩ := Finset.mem_image.mp hy
    rw [Finset.card_compl, Fintype.card_fin, card_eq_of_mem_fiber hz]
  · intro y hy y' hy' hne
    obtain ⟨z, hz, rfl⟩ := Finset.mem_image.mp hy
    obtain ⟨z', hz', rfl⟩ := Finset.mem_image.mp hy'
    have hcompl : zᶜ ∆ z'ᶜ = z ∆ z' := by
      ext a; simp [Finset.mem_symmDiff]
    rw [hcompl]
    refine hdist _ (symmDiff_mem_kernel hz hz') ?_
    have hzne : z ≠ z' := fun hq => hne (by rw [hq])
    simpa [symmDiff_eq_bot, Finset.bot_eq_empty] using hzne

end Spin
