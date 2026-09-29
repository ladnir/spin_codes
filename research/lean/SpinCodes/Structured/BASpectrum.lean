/-
The exact expected BA spectrum (`eq:structured-exact-ba-spectrum`):

    Ā_b(w) = Σ_{a=1}^{b} Σ_{c=0}^{b} G_b(a) P_b(a,c) P_b(c,w).

Everything it rests on is already proved:

* `prob_two_stage` — the two interleaved accumulators compose by
  Chapman-Kolmogorov, giving `Σ_c P(a,c)P(c,w)` for a fixed input of weight `a`;
* `card_tuples_weight` — the direct-sum block count is `[u^a] W(u)^k`, which is
  the paper's `G_b(a)`;
* `Golay.enumerator_eq` — what `W(u)` is for the Golay instantiation.

What is left is the bookkeeping: concatenate the blocks, group the nonzero
messages by their outer weight, and check that the sum really starts at `a = 1`.
The last point is the only one with content — it needs the base map to send
only the zero message to weight zero, which is the paper's injectivity
assumption on `H_a`.
-/
import SpinCodes.Structured.Composition
import SpinCodes.Structured.Enumerator

set_option linter.unusedSectionVars false

namespace Spin.Structured

open Finset

/-! ## Concatenating blocks -/

/-- The direct sum of `k` blocks of length `n`, as a word of length `k*n`. -/
def concatBlocks {k n : ℕ} (f : Fin k → (Fin n → Bool)) : Fin (k * n) → Bool :=
  fun j => f (finProdFinEquiv.symm j).1 (finProdFinEquiv.symm j).2

lemma wtF_concatBlocks {k n : ℕ} (f : Fin k → (Fin n → Bool)) :
    wtF (concatBlocks f) = ∑ i, wtF (f i) := by
  classical
  rw [wtF_eq_sum]
  rw [← Fintype.sum_equiv finProdFinEquiv
    (fun p : Fin k × Fin n => if f p.1 p.2 then 1 else 0)
    (fun j : Fin (k * n) => if concatBlocks f j then 1 else 0)
    (fun p => by simp [concatBlocks])]
  rw [Fintype.sum_prod_type]
  exact Finset.sum_congr rfl fun i _ => (wtF_eq_sum (f i)).symm

/-! ## The expected spectrum -/

variable {k n : ℕ} {σ : Type} [Fintype σ] [DecidableEq σ] [Zero σ]

/-- The weight of one encoded block. -/
def blockW (blockEnc : σ → (Fin n → Bool)) (s : σ) : ℕ := wtF (blockEnc s)

/-- The outer word of a direct-sum message. -/
def outerWord (blockEnc : σ → (Fin n → Bool)) (x : Fin k → σ) : Fin (k * n) → Bool :=
  concatBlocks (fun i => blockEnc (x i))

lemma wtF_outerWord (blockEnc : σ → (Fin n → Bool)) (x : Fin k → σ) :
    wtF (outerWord blockEnc x) = ∑ i, blockW blockEnc (x i) :=
  wtF_concatBlocks _

/-- **`eq:structured-exact-ba-spectrum`.**

The expected number of nonzero constituent words of weight `Wt`, over the two
BA permutations, is `Σ_{a≥1} G_b(a) Σ_c P_b(a,c) P_b(c,Wt)` with
`G_b(a) = [u^a] W(u)^k`.

The hypothesis `hzero` is the paper's injectivity assumption on the base map
`H_a`, and it is exactly what makes the outer sum start at `a = 1`. -/
theorem expected_spectrum (blockEnc : σ → (Fin n → Bool))
    (hzero : ∀ s : σ, blockW blockEnc s = 0 ↔ s = 0) (Wt : ℕ) :
    (∑ x ∈ univ.filter (fun x : Fin k → σ => x ≠ 0),
        ((FinPMF.uniform (Equiv.Perm (Fin (k * n)))).prod
            (FinPMF.uniform (Equiv.Perm (Fin (k * n))))).prob
          (fun τ => accWtF ((accF (outerWord blockEnc x ∘ τ.1)) ∘ τ.2) = Wt))
      = ∑ a ∈ Icc 1 (k * n),
          (((enumerator (blockW blockEnc)) ^ k).coeff a : ℝ)
            * ∑ c ∈ range (k * n + 1), Pt (k * n) a c * Pt (k * n) c Wt := by
  classical
  set b := k * n with hb
  set g : ℕ → ℝ := fun a => ∑ c ∈ range (b + 1), Pt b a c * Pt b c Wt with hg
  -- each message contributes according to the weight of its outer word
  have hterm : ∀ x : Fin k → σ,
      ((FinPMF.uniform (Equiv.Perm (Fin b))).prod
          (FinPMF.uniform (Equiv.Perm (Fin b)))).prob
        (fun τ => accWtF ((accF (outerWord blockEnc x ∘ τ.1)) ∘ τ.2) = Wt)
        = g (wtF (outerWord blockEnc x)) :=
    fun x => prob_two_stage (outerWord blockEnc x) rfl
  rw [Finset.sum_congr rfl (fun x _ => hterm x)]
  -- the outer weight of a nonzero message lies in `[1, b]`
  have hmaps : ∀ x ∈ univ.filter (fun x : Fin k → σ => x ≠ 0),
      wtF (outerWord blockEnc x) ∈ Icc 1 b := by
    intro x hx
    simp only [Finset.mem_filter, Finset.mem_univ, true_and] at hx
    refine Finset.mem_Icc.mpr ⟨?_, wtF_le _⟩
    by_contra hlt
    push_neg at hlt
    have hz : wtF (outerWord blockEnc x) = 0 := by omega
    rw [wtF_outerWord] at hz
    refine hx (funext fun i => ?_)
    have hi : blockW blockEnc (x i) = 0 :=
      Nat.eq_zero_of_le_zero (hz ▸ Finset.single_le_sum
        (f := fun i => blockW blockEnc (x i)) (fun j _ => Nat.zero_le _) (mem_univ i))
    exact (hzero (x i)).mp hi
  rw [Finset.sum_fiberwise_of_maps_to hmaps
      (fun x => g (wtF (outerWord blockEnc x))) |>.symm]
  refine Finset.sum_congr rfl fun a ha => ?_
  simp only [Finset.mem_Icc] at ha
  -- inside a weight class the summand is constant
  have hconst : ∀ x ∈ (univ.filter (fun x : Fin k → σ => x ≠ 0)).filter
      (fun x => wtF (outerWord blockEnc x) = a),
      g (wtF (outerWord blockEnc x)) = g a := by
    intro x hx
    simp only [Finset.mem_filter] at hx
    rw [hx.2]
  rw [Finset.sum_congr rfl hconst, Finset.sum_const, nsmul_eq_mul]
  congr 1
  -- for `a ≥ 1` the nonzero restriction is vacuous, and the count is `[u^a] W^k`
  have hset : (univ.filter (fun x : Fin k → σ => x ≠ 0)).filter
      (fun x => wtF (outerWord blockEnc x) = a)
      = univ.filter (fun x : Fin k → σ => ∑ i, blockW blockEnc (x i) = a) := by
    ext x
    simp only [Finset.mem_filter, Finset.mem_univ, true_and, wtF_outerWord]
    constructor
    · rintro ⟨-, h2⟩; exact h2
    · intro h2
      refine ⟨?_, h2⟩
      rintro rfl
      have : ∑ i, blockW blockEnc ((0 : Fin k → σ) i) = 0 := by
        refine Finset.sum_eq_zero fun i _ => ?_
        exact (hzero _).mpr rfl
      omega
  rw [hset, card_tuples_weight (blockW blockEnc) k a]

end Spin.Structured
