/-
The direct-sum weight enumerator.

The scalable outer takes `H_{a,b}`, the direct sum of `b/a` copies of a base
map `H_a`, and the paper writes its weight distribution as

    G_b(a) := [u^a] G(u)^{b/24}

for the Golay instantiation.  The general fact behind that notation is
independent of Golay and of coding theory: the number of `k`-tuples of block
messages whose total weight is `a` equals the `a`-th coefficient of the `k`-th
power of the per-block enumerator.

Stated over an arbitrary finite block-message type with an arbitrary weight
function, which is also how the paper sets it up ("an even base length `a` and
an injective rate-one-half map `H_a` ... its weight enumerator `W_a(z)`").
-/
import Mathlib

set_option linter.unusedSectionVars false

namespace Spin.Structured

open Finset Polynomial

variable {σ : Type*} [Fintype σ] [DecidableEq σ]

/-- The weight enumerator of a block: `W(u) = Σ_s u^{w(s)}`. -/
noncomputable def enumerator (w : σ → ℕ) : Polynomial ℕ :=
  ∑ s : σ, (X : Polynomial ℕ) ^ (w s)

lemma coeff_enumerator (w : σ → ℕ) (i : ℕ) :
    (enumerator w).coeff i = (univ.filter (fun s : σ => w s = i)).card := by
  classical
  unfold enumerator
  rw [finsetSum_coeff, Finset.card_filter]
  refine Finset.sum_congr rfl fun s _ => ?_
  rw [Polynomial.coeff_X_pow]
  by_cases h : w s = i
  · simp [h]
  · rw [if_neg (fun he => h he.symm), if_neg h]

/-- Splitting a tuple on its first coordinate. -/
lemma card_filter_cons_split (k : ℕ) (P : (Fin (k + 1) → σ) → Prop) [DecidablePred P] :
    (univ.filter P).card
      = ∑ s : σ, (univ.filter (fun g : Fin k → σ => P (Fin.cons s g))).card := by
  classical
  simp only [Finset.card_filter]
  rw [← Fintype.sum_equiv (Fin.consEquiv (fun _ : Fin (k + 1) => σ))
    (fun p : σ × (Fin k → σ) => if P (Fin.cons p.1 p.2) then 1 else 0)
    (fun f => if P f then 1 else 0) (fun p => rfl)]
  rw [Fintype.sum_prod_type]

/-- **The direct-sum enumerator.**  The number of `k`-tuples of block messages
of total weight `a` is `[u^a] W(u)^k`. -/
theorem card_tuples_weight (w : σ → ℕ) (k a : ℕ) :
    (univ.filter (fun f : Fin k → σ => ∑ i, w (f i) = a)).card
      = ((enumerator w) ^ k).coeff a := by
  classical
  induction k generalizing a with
  | zero =>
      simp only [pow_zero, Polynomial.coeff_one]
      have hsum : ∀ f : Fin 0 → σ, ∑ i, w (f i) = 0 := by intro f; simp
      rcases Nat.eq_zero_or_pos a with rfl | ha
      · rw [Finset.filter_true_of_mem (fun f _ => hsum f)]
        simp
      · rw [Finset.filter_false_of_mem (fun f _ => by rw [hsum f]; omega),
          if_neg (by omega)]
        simp
  | succ k ih =>
      rw [card_filter_cons_split, pow_succ', Polynomial.coeff_mul]
      -- left side: group the tuples by their first block
      have hL : ∀ s : σ,
          (univ.filter (fun g : Fin k → σ =>
              ∑ i, w ((Fin.cons s g : Fin (k + 1) → σ) i) = a)).card
            = if w s ≤ a then ((enumerator w) ^ k).coeff (a - w s) else 0 := by
        intro s
        have hcons : ∀ g : Fin k → σ,
            ∑ i, w ((Fin.cons s g : Fin (k + 1) → σ) i) = w s + ∑ i, w (g i) := by
          intro g
          rw [Fin.sum_univ_succ]
          simp
        by_cases hs : w s ≤ a
        · rw [if_pos hs, ← ih (a - w s)]
          refine congrArg Finset.card (Finset.filter_congr fun g _ => ?_)
          rw [hcons g]
          constructor
          · intro h; omega
          · intro h; omega
        · rw [if_neg hs]
          refine Finset.card_eq_zero.mpr (Finset.filter_eq_empty_iff.mpr fun g _ => ?_)
          rw [hcons g]
          omega
      rw [Finset.sum_congr rfl (fun s _ => hL s)]
      -- right side: group the antidiagonal by the weight of the first block
      rw [Finset.Nat.sum_antidiagonal_eq_sum_range_succ
        (fun i j => (enumerator w).coeff i * ((enumerator w) ^ k).coeff j)]
      have hR : ∀ i ∈ range (a + 1),
          (enumerator w).coeff i * ((enumerator w) ^ k).coeff (a - i)
            = ∑ s ∈ univ.filter (fun s : σ => w s = i),
                ((enumerator w) ^ k).coeff (a - w s) := by
        intro i _
        rw [coeff_enumerator]
        rw [Finset.sum_congr rfl (fun s hs => by
          simp only [Finset.mem_filter] at hs
          rw [hs.2])]
        rw [Finset.sum_const, smul_eq_mul]
      rw [Finset.sum_congr rfl hR]
      -- both sides are now the sum over blocks of weight at most `a`
      have hLHS : (∑ s : σ, if w s ≤ a then ((enumerator w) ^ k).coeff (a - w s) else 0)
          = ∑ s ∈ univ.filter (fun s : σ => w s ≤ a),
              ((enumerator w) ^ k).coeff (a - w s) := by
        rw [Finset.sum_filter]
      have hmaps : ∀ s ∈ univ.filter (fun s : σ => w s ≤ a), w s ∈ range (a + 1) := by
        intro s hs
        simp only [Finset.mem_filter, Finset.mem_univ, true_and] at hs
        simp only [Finset.mem_range]
        omega
      have hRHS : (∑ i ∈ range (a + 1), ∑ s ∈ univ.filter (fun s : σ => w s = i),
            ((enumerator w) ^ k).coeff (a - w s))
          = ∑ s ∈ univ.filter (fun s : σ => w s ≤ a),
              ((enumerator w) ^ k).coeff (a - w s) := by
        rw [← Finset.sum_fiberwise_of_maps_to hmaps
          (fun s => ((enumerator w) ^ k).coeff (a - w s))]
        refine Finset.sum_congr rfl fun i hi => ?_
        simp only [Finset.mem_range] at hi
        refine Finset.sum_congr ?_ fun _ _ => rfl
        ext s
        simp only [Finset.mem_filter, Finset.mem_univ, true_and]
        constructor
        · intro h; exact ⟨by omega, h⟩
        · rintro ⟨-, h⟩; exact h
      rw [hLHS, hRHS]


/-! ## Sanity check

For a single bit, the enumerator is `1 + u`, so `card_tuples_weight` specialises
to "the number of weight-`a` binary words of length `k` is `[u^a](1+u)^k`".
This pins the definition, which the general theorem alone does not. -/

example : enumerator (fun b : Bool => if b then 1 else 0) = 1 + X := by
  unfold enumerator
  rw [Fintype.sum_bool]
  simp [add_comm]

end Spin.Structured
