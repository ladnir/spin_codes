/-
One-sample selection of the shared outer constituent
(paper: `lem:structured-one-sample-selection`).

Setup samples the BA-3 constituent `𝒪_b` *once* and reuses that realization
in every outer position, so the whole proof must pay for a single good event
`𝒢_b` determined by the two BA permutations.  On `𝒢_b`:

* no nonzero constituent word has relative weight outside `𝒲`;
* every remaining weight satisfies `A_{𝒪_b}(w) ≤ b² · Ā_b(w)`.

The paper states `Pr[𝒢_b] = 1 - o(1)`.  We prove the quantitative statement
it is derived from — the `o(1)` is then immediate and is recorded as a
corollary.  Keeping the finite-`b` inequality explicit is what lets the
downstream union be checked.
-/
import SpinCodes.Prob

set_option linter.unusedSectionVars false

namespace Spin

open Finset

attribute [local instance] Classical.propDecidable

variable {Ω : Type*} [Fintype Ω] (P : FinPMF Ω) (A : Ω → ℕ → ℕ)

/-- `Ā_b(w) = 𝔼_{τ₁,τ₂}[A_{𝒪_b}(w)]`, the expected spectrum. -/
def Abar (w : ℕ) : ℝ := P.expect (fun ω => (A ω w : ℝ))

lemma Abar_nonneg (w : ℕ) : 0 ≤ Abar P A w :=
  FinPMF.expect_nonneg _ fun _ => by positivity

/-- The weights outside the certified window `𝒲`. -/
def badWeights (b : ℕ) (Wgood : Finset ℕ) : Finset ℕ := (range (b + 1)) \ Wgood

/-- The good event `𝒢_b`. -/
def Good (b : ℕ) (Wgood : Finset ℕ) (ω : Ω) : Prop :=
  (∀ w ∈ badWeights b Wgood, A ω w = 0) ∧
  (∀ w ∈ Wgood, (A ω w : ℝ) ≤ (b : ℝ) ^ 2 * Abar P A w)

/-- Per-weight Markov bound for the selected spectrum.  The degenerate case
`Ā_b(w) = 0` is handled by the counting form of Markov, not by dividing by
zero. -/
lemma prob_spectrum_fail_le (b : ℕ) (hb : 0 < b) (w : ℕ) :
    P.prob (fun ω => (b : ℝ) ^ 2 * Abar P A w < (A ω w : ℝ)) ≤ 1 / (b : ℝ) ^ 2 := by
  have hbR : (0 : ℝ) < (b : ℝ) ^ 2 := by positivity
  rcases eq_or_lt_of_le (Abar_nonneg P A w) with hz | hpos
  · -- expected multiplicity zero: no word of this weight survives
    have h1 : P.prob (fun ω => (b : ℝ) ^ 2 * Abar P A w < (A ω w : ℝ))
        ≤ P.prob (fun ω => 1 ≤ A ω w) := by
      refine FinPMF.prob_mono P fun ω hω => ?_
      rw [← hz, mul_zero] at hω
      exact_mod_cast Nat.one_le_iff_ne_zero.mpr (by exact_mod_cast hω.ne')
    have h2 := FinPMF.markov_one P (fun ω => A ω w)
    have : P.prob (fun ω => 1 ≤ A ω w) ≤ 0 :=
      le_trans h2 (le_of_eq hz.symm)
    calc P.prob (fun ω => (b : ℝ) ^ 2 * Abar P A w < (A ω w : ℝ)) ≤ 0 := le_trans h1 this
      _ ≤ 1 / (b : ℝ) ^ 2 := by positivity
  · -- positive expectation: ordinary Markov at threshold `b² · Ā_b(w)`
    have hthr : 0 < (b : ℝ) ^ 2 * Abar P A w := by positivity
    have hmk := FinPMF.markov P (f := fun ω => (A ω w : ℝ)) (fun ω => by positivity) hthr
    have hmono : P.prob (fun ω => (b : ℝ) ^ 2 * Abar P A w < (A ω w : ℝ))
        ≤ P.prob (fun ω => (b : ℝ) ^ 2 * Abar P A w ≤ (A ω w : ℝ)) :=
      FinPMF.prob_mono P fun ω hω => le_of_lt hω
    refine le_trans hmono (le_trans hmk (le_of_eq ?_))
    show Abar P A w / ((b : ℝ) ^ 2 * Abar P A w) = 1 / (b : ℝ) ^ 2
    field_simp

/-- **One-sample selection, quantitative form.**

`Pr[¬𝒢_b] ≤ (tail mass of the expected spectrum) + (b+1)/b²`.

The first term is Markov applied to the expected number of nonzero tail
words; the second is a union bound over at most `b+1` weights, each costing
`b⁻²`. -/
theorem prob_not_good_le (b : ℕ) (hb : 0 < b) (Wgood : Finset ℕ)
    (hsub : Wgood ⊆ range (b + 1)) :
    P.prob (fun ω => ¬ Good P A b Wgood ω)
      ≤ (∑ w ∈ badWeights b Wgood, Abar P A w) + (b + 1) / (b : ℝ) ^ 2 := by
  classical
  -- `¬𝒢` splits into a tail failure and a spectrum failure
  have hsplit : ∀ ω, ¬ Good P A b Wgood ω →
      (∃ w ∈ badWeights b Wgood, 1 ≤ A ω w) ∨
      (∃ w ∈ Wgood, (b : ℝ) ^ 2 * Abar P A w < (A ω w : ℝ)) := by
    intro ω hω
    unfold Good at hω
    rw [not_and_or] at hω
    rcases hω with h | h
    · left
      push_neg at h
      obtain ⟨w, hw, hne⟩ := h
      exact ⟨w, hw, Nat.one_le_iff_ne_zero.mpr hne⟩
    · right
      push_neg at h
      obtain ⟨w, hw, hlt⟩ := h
      exact ⟨w, hw, hlt⟩
  refine le_trans (FinPMF.prob_mono P hsplit) (le_trans (FinPMF.prob_or_le P _ _) ?_)
  refine add_le_add ?_ ?_
  · -- tail: Markov on the expected number of nonzero tail words
    refine le_trans (FinPMF.prob_biUnion_le P _ (fun w ω => 1 ≤ A ω w)) ?_
    exact Finset.sum_le_sum fun w _ => FinPMF.markov_one P (fun ω => A ω w)
  · -- spectrum: union bound over the certified window
    refine le_trans (FinPMF.prob_biUnion_le P _
      (fun w ω => (b : ℝ) ^ 2 * Abar P A w < (A ω w : ℝ))) ?_
    refine le_trans (Finset.sum_le_sum fun w _ => prob_spectrum_fail_le P A b hb w) ?_
    rw [Finset.sum_const, nsmul_eq_mul]
    have hcard : (Wgood.card : ℝ) ≤ (b : ℝ) + 1 := by
      have := Finset.card_le_card hsub
      rw [Finset.card_range] at this
      exact_mod_cast this
    have hinv : (0:ℝ) ≤ ((b : ℝ) ^ 2)⁻¹ := by positivity
    simpa [div_eq_mul_inv] using mul_le_mul_of_nonneg_right hcard hinv

/-- The paper's form: if the certified tail mass vanishes and `b → ∞`, then
`Pr[𝒢_b] = 1 - o(1)`. -/
theorem prob_not_good_tendsto_zero
    {Ωm : ℕ → Type*} [∀ m, Fintype (Ωm m)]
    (Pm : ∀ m, FinPMF (Ωm m)) (Am : ∀ m, Ωm m → ℕ → ℕ)
    (bm : ℕ → ℕ) (Wg : ℕ → Finset ℕ)
    (hbpos : ∀ m, 0 < bm m) (hsub : ∀ m, Wg m ⊆ range (bm m + 1))
    (hb : Filter.Tendsto bm Filter.atTop Filter.atTop)
    (htail : Filter.Tendsto
      (fun m => ∑ w ∈ badWeights (bm m) (Wg m), Abar (Pm m) (Am m) w)
      Filter.atTop (nhds 0)) :
    Filter.Tendsto
      (fun m => (Pm m).prob (fun ω => ¬ Good (Pm m) (Am m) (bm m) (Wg m) ω))
      Filter.atTop (nhds 0) := by
  have hnonneg : ∀ m, 0 ≤ (Pm m).prob (fun ω => ¬ Good (Pm m) (Am m) (bm m) (Wg m) ω) :=
    fun m => FinPMF.prob_nonneg _ _
  have hub : ∀ m, (Pm m).prob (fun ω => ¬ Good (Pm m) (Am m) (bm m) (Wg m) ω)
      ≤ (∑ w ∈ badWeights (bm m) (Wg m), Abar (Pm m) (Am m) w)
        + ((bm m : ℝ) + 1) / (bm m : ℝ) ^ 2 :=
    fun m => prob_not_good_le (Pm m) (Am m) (bm m) (hbpos m) (Wg m) (hsub m)
  have hbne : ∀ m, ((bm m : ℝ)) ≠ 0 := fun m => Nat.cast_ne_zero.mpr (hbpos m).ne'
  have hcast : Filter.Tendsto (fun m => (bm m : ℝ)) Filter.atTop Filter.atTop :=
    tendsto_natCast_atTop_atTop.comp hb
  have hinv : Filter.Tendsto (fun m => ((bm m : ℝ))⁻¹) Filter.atTop (nhds 0) :=
    hcast.inv_tendsto_atTop
  have hsecond : Filter.Tendsto (fun m => ((bm m : ℝ) + 1) / (bm m : ℝ) ^ 2)
      Filter.atTop (nhds 0) := by
    have h0 : Filter.Tendsto
        (fun m => ((bm m : ℝ))⁻¹ + ((bm m : ℝ))⁻¹ * ((bm m : ℝ))⁻¹)
        Filter.atTop (nhds 0) := by simpa using hinv.add (hinv.mul hinv)
    refine h0.congr fun m => ?_
    have h := hbne m
    field_simp
  have hsum := htail.add hsecond
  rw [add_zero] at hsum
  exact tendsto_of_tendsto_of_tendsto_of_le_of_le tendsto_const_nhds hsum hnonneg hub

end Spin
