/-
Fourier inversion and Parseval for a binary dual code.

`app:imt-finite-transfers` obtains its integer fibre bounds from

    k_j = 2^{-19} (n_j + ∑_w b_w K_j(w)),
    V_j = 2^{-19} (n_j² + ∑_w b_w K_j(w)²),

where `k_j` counts weight-`j` kernel words and `V_j` is the sum of squares of
the weight-`j` syndrome fibre sizes.

Both are proved here *without* the syndrome map.  For `k_j` that is standard.
For `V_j` it is the observation that two words share a syndrome exactly when
their symmetric difference lies in the kernel, so

    ∑_s f_j(s)² = #{(x,y) : |x| = |y| = j, x ∆ y ∈ ker},

which mentions only the kernel.  Since the kernel is itself determined by the
dual code, the entire development needs nothing but a subgroup `D` of
`(Finset (Fin n), ∆)` — no matrix, no syndrome space, no quotient.

The character is `v ↦ (-1)^|v ∩ x|`, and the one fact making it a character is
that `|·|` is additive mod 2 under `∆`.
-/
import SpinCodes.Structured.Krawtchouk

namespace Spin

open Finset
open scoped symmDiff

/-! ## The character -/

/-- `(-1)^{⟨u,x⟩}`, the nontrivial character of `𝔽₂` applied to `|u ∩ x|`. -/
def chi {n : ℕ} (u x : Finset (Fin n)) : ℤ := (-1 : ℤ) ^ (u ∩ x).card

/-- Weight of a symmetric difference, without truncated subtraction. -/
lemma card_symmDiff_add_two_mul_card_inter {n : ℕ} (a b : Finset (Fin n)) :
    (a ∆ b).card + 2 * (a ∩ b).card = a.card + b.card := by
  have h1 : (a ∆ b).card = (a \ b).card + (b \ a).card := by
    rw [Finset.symmDiff_def]
    exact Finset.card_union_of_disjoint disjoint_sdiff_sdiff
  have h2 : (a \ b).card + (a ∩ b).card = a.card := Finset.card_sdiff_add_card_inter a b
  have h3 : (b \ a).card + (b ∩ a).card = b.card := Finset.card_sdiff_add_card_inter b a
  have h4 : b ∩ a = a ∩ b := Finset.inter_comm b a
  rw [h4] at h3
  omega

/-- Weight is additive mod 2 under symmetric difference. -/
lemma neg_one_pow_card_symmDiff {n : ℕ} (a b : Finset (Fin n)) :
    (-1 : ℤ) ^ (a ∆ b).card = (-1 : ℤ) ^ a.card * (-1 : ℤ) ^ b.card := by
  have key := card_symmDiff_add_two_mul_card_inter a b
  have hsq : ((-1 : ℤ) ^ 2) ^ (a ∩ b).card = 1 := by norm_num
  have : (-1 : ℤ) ^ (a ∆ b).card * ((-1 : ℤ) ^ 2) ^ (a ∩ b).card
      = (-1 : ℤ) ^ a.card * (-1 : ℤ) ^ b.card := by
    rw [← pow_mul, ← pow_add, key, pow_add]
  rwa [hsq, mul_one] at this

lemma symmDiff_inter_left {n : ℕ} (u v x : Finset (Fin n)) :
    (u ∆ v) ∩ x = (u ∩ x) ∆ (v ∩ x) := by
  ext a; simp [Finset.mem_symmDiff]; tauto

lemma inter_symmDiff_right {n : ℕ} (u x y : Finset (Fin n)) :
    u ∩ (x ∆ y) = (u ∩ x) ∆ (u ∩ y) := by
  ext a; simp [Finset.mem_symmDiff]; tauto

/-- The character is multiplicative in its first argument. -/
lemma chi_symmDiff_left {n : ℕ} (u v x : Finset (Fin n)) :
    chi (u ∆ v) x = chi u x * chi v x := by
  unfold chi
  rw [symmDiff_inter_left u v x]
  exact neg_one_pow_card_symmDiff _ _

/-- And in its second — which is what turns Parseval into a statement about
`x ∆ y`. -/
lemma chi_symmDiff_right {n : ℕ} (u x y : Finset (Fin n)) :
    chi u (x ∆ y) = chi u x * chi u y := by
  unfold chi
  rw [inter_symmDiff_right u x y]
  exact neg_one_pow_card_symmDiff _ _

/-! ## Orthogonality

A character of a finite group is either trivial or sums to zero.  Here the
group is `(D, ∆)` and the character is `chi · x`. -/

/-- `x` is orthogonal to every element of `D`. -/
def OrthTo {n : ℕ} (D : Finset (Finset (Fin n))) (x : Finset (Fin n)) : Prop :=
  ∀ v ∈ D, Even (v ∩ x).card

instance {n : ℕ} (D : Finset (Finset (Fin n))) (x : Finset (Fin n)) :
    Decidable (OrthTo D x) := by unfold OrthTo; infer_instance

theorem sum_chi {n : ℕ} (D : Finset (Finset (Fin n)))
    (hD : ∀ u ∈ D, ∀ v ∈ D, u ∆ v ∈ D) (x : Finset (Fin n)) :
    ∑ v ∈ D, chi v x = if OrthTo D x then (D.card : ℤ) else 0 := by
  classical
  by_cases horth : OrthTo D x
  · rw [if_pos horth]
    rw [Finset.sum_congr rfl (fun v hv => ?_), Finset.sum_const, nsmul_eq_mul, mul_one]
    exact (horth v hv).neg_one_pow
  · rw [if_neg horth]
    obtain ⟨u, huD, hu⟩ : ∃ u ∈ D, ¬ Even (u ∩ x).card := by
      by_contra hc
      exact horth fun v hv => by
        by_contra hev
        exact hc ⟨v, hv, hev⟩
    have hchiu : chi u x = -1 := by
      unfold chi
      exact Odd.neg_one_pow (Nat.not_even_iff_odd.mp hu)
    have hreindex : ∑ v ∈ D, chi v x = ∑ v ∈ D, chi (u ∆ v) x := by
      refine (Finset.sum_nbij' (fun v => u ∆ v) (fun v => u ∆ v) ?_ ?_ ?_ ?_ ?_).symm
      · intro v hv
        exact hD u huD v hv
      · intro v hv
        exact hD u huD v hv
      · intro v _
        exact symmDiff_symmDiff_cancel_left u v
      · intro v _
        exact symmDiff_symmDiff_cancel_left u v
      · intro v _
        rfl
    have hsplit : ∑ v ∈ D, chi (u ∆ v) x = chi u x * ∑ v ∈ D, chi v x := by
      rw [Finset.mul_sum]
      exact Finset.sum_congr rfl fun v _ => chi_symmDiff_left u v x
    rw [hsplit, hchiu] at hreindex
    linarith [hreindex]

/-! ## The two identities -/

/-- The weight-`j` layer. -/
abbrev layer (n j : ℕ) : Finset (Finset (Fin n)) :=
  (univ : Finset (Fin n)).powersetCard j

lemma sum_chi_layer {n : ℕ} (v : Finset (Fin n)) (j : ℕ) :
    ∑ x ∈ layer n j, chi v x = krawtchouk n j v.card := by
  rw [← krawtchouk_char_sum v j]
  exact Finset.sum_congr rfl fun x _ => by unfold chi; rw [Finset.inter_comm]

/-- The weight-`j` words sharing `x`'s syndrome: those `y` whose symmetric
difference with `x` lies in the kernel. -/
def fiber {n : ℕ} (D : Finset (Finset (Fin n))) (j : ℕ) (x : Finset (Fin n)) :
    Finset (Finset (Fin n)) :=
  (layer n j).filter (fun y => OrthTo D (x ∆ y))

/-- **Fourier inversion at an arbitrary syndrome**:

    |D| · |fibre through x| = ∑_{v ∈ D} χ(v,x) · K_j(|v|).

At `x = ∅` every character is `1` and this is the kernel count. -/
theorem card_fiber {n : ℕ} (D : Finset (Finset (Fin n)))
    (hD : ∀ u ∈ D, ∀ v ∈ D, u ∆ v ∈ D) (j : ℕ) (x : Finset (Fin n)) :
    (D.card : ℤ) * ((fiber D j x).card : ℤ)
      = ∑ v ∈ D, chi v x * krawtchouk n j v.card := by
  classical
  calc (D.card : ℤ) * ((fiber D j x).card : ℤ)
      = ∑ y ∈ layer n j, (if OrthTo D (x ∆ y) then (D.card : ℤ) else 0) := by
        unfold fiber
        rw [Finset.sum_ite, Finset.sum_const, Finset.sum_const]
        simp [mul_comm]
    _ = ∑ y ∈ layer n j, ∑ v ∈ D, chi v (x ∆ y) :=
        Finset.sum_congr rfl fun y _ => (sum_chi D hD _).symm
    _ = ∑ v ∈ D, ∑ y ∈ layer n j, chi v (x ∆ y) := Finset.sum_comm
    _ = ∑ v ∈ D, chi v x * krawtchouk n j v.card := by
        refine Finset.sum_congr rfl fun v _ => ?_
        rw [← sum_chi_layer v j, Finset.mul_sum]
        exact Finset.sum_congr rfl fun y _ => chi_symmDiff_right v x y

/-- **Fourier inversion**: `|D| · k_j = ∑_{v ∈ D} K_j(|v|)`.

Splitting the right-hand sum at `v = ∅` (where `K_j(0) = n_j`) and grouping
the rest by weight gives the appendix's `n_j + ∑_w b_w K_j(w)`. -/
theorem card_orth_layer {n : ℕ} (D : Finset (Finset (Fin n)))
    (hD : ∀ u ∈ D, ∀ v ∈ D, u ∆ v ∈ D) (j : ℕ) :
    (D.card : ℤ) * (((layer n j).filter (fun x => OrthTo D x)).card : ℤ)
      = ∑ v ∈ D, krawtchouk n j v.card := by
  have h := card_fiber D hD j ∅
  have hemp : ∀ y : Finset (Fin n), (∅ : Finset (Fin n)) ∆ y = y := by
    intro y; ext a; simp [Finset.mem_symmDiff]
  simpa [fiber, chi, hemp] using h

/-- **Parseval**: `|D| · V_j = ∑_{v ∈ D} K_j(|v|)²`, where `V_j` is written in
its coset form — two weight-`j` words share a syndrome exactly when their
symmetric difference is orthogonal to the whole dual code. -/
theorem card_orth_pairs {n : ℕ} (D : Finset (Finset (Fin n)))
    (hD : ∀ u ∈ D, ∀ v ∈ D, u ∆ v ∈ D) (j : ℕ) :
    (D.card : ℤ)
        * ((((layer n j) ×ˢ (layer n j)).filter
            (fun p => OrthTo D (p.1 ∆ p.2))).card : ℤ)
      = ∑ v ∈ D, (krawtchouk n j v.card) ^ 2 := by
  classical
  calc (D.card : ℤ)
        * ((((layer n j) ×ˢ (layer n j)).filter
            (fun p => OrthTo D (p.1 ∆ p.2))).card : ℤ)
      = ∑ p ∈ (layer n j) ×ˢ (layer n j),
          (if OrthTo D (p.1 ∆ p.2) then (D.card : ℤ) else 0) := by
        rw [Finset.sum_ite, Finset.sum_const, Finset.sum_const]
        simp [mul_comm]
    _ = ∑ p ∈ (layer n j) ×ˢ (layer n j), ∑ v ∈ D, chi v (p.1 ∆ p.2) :=
        Finset.sum_congr rfl fun p _ => (sum_chi D hD _).symm
    _ = ∑ v ∈ D, ∑ p ∈ (layer n j) ×ˢ (layer n j), chi v (p.1 ∆ p.2) :=
        Finset.sum_comm
    _ = ∑ v ∈ D, (krawtchouk n j v.card) ^ 2 := by
        refine Finset.sum_congr rfl fun v _ => ?_
        rw [Finset.sum_product]
        calc ∑ x ∈ layer n j, ∑ y ∈ layer n j, chi v (x ∆ y)
            = ∑ x ∈ layer n j, ∑ y ∈ layer n j, chi v x * chi v y :=
              Finset.sum_congr rfl fun x _ =>
                Finset.sum_congr rfl fun y _ => chi_symmDiff_right v x y
          _ = (∑ x ∈ layer n j, chi v x) * (∑ y ∈ layer n j, chi v y) := by
              rw [Finset.sum_mul_sum]
          _ = (krawtchouk n j v.card) ^ 2 := by
              rw [sum_chi_layer]; ring

end Spin
