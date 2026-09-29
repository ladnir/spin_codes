/-
`lem:structured-exact-ba-spectrum`, instantiated at Golay.

    Ā_b(w) = Σ_{a=1}^{b} Σ_{c=0}^{b} G_b(a) P_b(a,c) P_b(c,w),
    G_b(a) = [u^a] G(u)^{b/24},
    G(u)   = 1 + 759u^8 + 2576u^12 + 759u^16 + u^24.

Every ingredient is proved rather than assumed: `G(u)` is computed from the
4096 codewords, `G_b(a)` is the direct-sum convolution, `P_b(a,c)` is the
uniform-interleaver transition law, and the two stages compose by
Chapman-Kolmogorov.
-/
import SpinCodes.Structured.BASpectrum
import SpinCodes.Structured.GolayEnum

set_option linter.unusedSectionVars false

namespace Spin.Structured.Golay

open Finset Polynomial Spin.Structured

/-- The zero message encodes to the zero word. -/
lemma W_zero : W (0 : Fin 12 → Bool) = 0 := by decide

/-- The base map is injective at zero: only the zero message has weight zero.
This is the content of the paper's injectivity assumption on `H_a`, and here it
is a consequence of the computed distribution (`A 0 = 1`) rather than a
hypothesis. -/
lemma W_eq_zero_iff (m : Fin 12 → Bool) : W m = 0 ↔ m = 0 := by
  classical
  constructor
  · intro h
    have hcard : (univ.filter (fun m' : Fin 12 → Bool => W m' = 0)).card = 1 := by
      rw [card_weight 0]
      rfl
    have hm : m ∈ univ.filter (fun m' : Fin 12 → Bool => W m' = 0) := by
      simp [h]
    have hz : (0 : Fin 12 → Bool) ∈ univ.filter (fun m' : Fin 12 → Bool => W m' = 0) := by
      simp [W_zero]
    exact Finset.card_le_one.mp (le_of_eq hcard) m hm 0 hz
  · rintro rfl
    exact W_zero

lemma blockW_enc : blockW enc = W := rfl

/-- **`lem:structured-exact-ba-spectrum` for the Golay-BA-3 constituent.**

The expected number of nonzero constituent words of weight `Wt`, over the two
BA permutations, with the outer enumerator computed rather than assumed. -/
theorem ba_expected_spectrum (k Wt : ℕ) :
    (∑ x ∈ univ.filter (fun x : Fin k → (Fin 12 → Bool) => x ≠ 0),
        ((FinPMF.uniform (Equiv.Perm (Fin (k * 24)))).prod
            (FinPMF.uniform (Equiv.Perm (Fin (k * 24))))).prob
          (fun τ => accWtF ((accF (outerWord enc x ∘ τ.1)) ∘ τ.2) = Wt))
      = ∑ a ∈ Icc 1 (k * 24),
          ((((1 : Polynomial ℕ) + C 759 * X ^ 8 + C 2576 * X ^ 12
              + C 759 * X ^ 16 + X ^ 24) ^ k).coeff a : ℝ)
            * ∑ c ∈ range (k * 24 + 1), Pt (k * 24) a c * Pt (k * 24) c Wt := by
  rw [← enumerator_eq]
  exact expected_spectrum enc (fun s => by rw [blockW_enc]; exact W_eq_zero_iff s) Wt

end Spin.Structured.Golay
