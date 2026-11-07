import SpinCodes.UniformMixer.ScalableSpectrumPolynomial
import SpinCodes.UniformMixer.ScalableTransferCollatz

/-!
# Paper-oriented transfer and its Collatz interface

The paper evolves row vectors, whereas the Lean envelope evolves column
vectors. Its displayed transfer is the transpose of `scalableThreeTransfer`.
The sparse witness therefore acts on the transposed matrix and starts from
the all-one terminal vector.
-/

noncomputable section

namespace Spin.UniformMixer

def paperSpectrumTransfer (β z : ℝ) : Matrix (Fin 3) (Fin 3) ℝ :=
  (spectrumTransfer β z).transpose

theorem paperSpectrumTransfer_eq_actual (β z : ℝ)
    (hβ0 : 0 ≤ β) (hβ1 : β ≤ 1) :
    paperSpectrumTransfer β z =
      (scalableThreeTransfer β z hβ0 hβ1).transpose := by
  rw [scalableThreeTransfer_eq_spectrum]
  rfl

theorem transfer_total_eq_paper_first
    (T : Matrix (Fin 3) (Fin 3) ℝ) (n : ℕ) :
    (∑ i : Fin 3, (T ^ n).mulVec ![1, 0, 0] i) =
      (((T.transpose) ^ n).mulVec (fun _ => 1) 0) := by
  rw [← Matrix.transpose_pow]
  simp [Matrix.mulVec, dotProduct, Fin.sum_univ_three]

theorem transfer_paper_collatz_total
    (T : Matrix (Fin 3) (Fin 3) ℝ)
    (hT : ∀ i j, 0 ≤ T i j)
    (v : Fin 3 → ℝ) (c ρ : ℝ)
    (hc : 0 ≤ c) (hρ : 0 ≤ ρ)
    (hones : ∀ i, (1 : ℝ) ≤ c * v i)
    (hv : ∀ i, T.transpose.mulVec v i ≤ ρ * v i)
    (n : ℕ) :
    ∑ i : Fin 3, (T ^ n).mulVec ![1, 0, 0] i ≤
      c * ρ ^ n * v 0 := by
  rw [transfer_total_eq_paper_first]
  exact transfer_collatz_pointwise T.transpose
    (fun i j => hT j i) v (fun _ => 1) c ρ hc hρ hones hv n 0

theorem scalable_bernoulli_paper_collatz_bound (β z : ℝ)
    (hβ0 : 0 ≤ β) (hβ1 : β ≤ 1)
    (hz0 : 0 ≤ z) (hz1 : z ≤ 1)
    (v : Fin 3 → ℝ) (c ρ : ℝ)
    (hc : 0 ≤ c) (hρ : 0 ≤ ρ)
    (hones : ∀ i, (1 : ℝ) ≤ c * v i)
    (hv : ∀ i,
      (paperSpectrumTransfer β z).mulVec v i ≤ ρ * v i)
    (n : ℕ) :
    ∑ r : State22, (bernoulliStepKernel β z hβ0 hβ1 ^ n) 0 r ≤
      c * ρ ^ n * v 0 := by
  calc
    _ ≤ scalableTransferMoment β z hβ0 hβ1 n :=
      scalable_bernoulli_transfer_bound β z hβ0 hβ1 hz0 hz1 n
    _ ≤ c * ρ ^ n * v 0 := by
      unfold scalableTransferMoment
      apply transfer_paper_collatz_total
        (scalableThreeTransfer β z hβ0 hβ1)
        (scalableThreeTransfer_nonneg β z hβ0 hβ1 hz0 hz1)
        v c ρ hc hρ hones
      intro i
      rw [← paperSpectrumTransfer_eq_actual β z hβ0 hβ1]
      exact hv i

end Spin.UniformMixer
