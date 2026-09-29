import SpinCodes.Structured.SparseRateNative
import SpinCodes.Structured.ConcreteOuterCountingSelected

noncomputable section
namespace Spin.Structured.SparseRate
open Finset ConcreteOuter

/-- Uniform sparse first moment for each fixed shared outer satisfying its stated row envelope. -/
theorem outer_occupation_sparse_rate {L k R d Q : ℕ} (seed : Seed k)
    (e : (Fin (k*24) × Fin L) ≃ (Fin R × Fin 128))
    (hQ : 4096 ≤ Q) (hQL : Q ≤ L) (hb : 1000 ≤ k*24)
    {B α ε : ℝ} (hB0 : 0 ≤ B)
    (hB : ∀ w ≤ k*24, (spectrum seed w:ℝ) ≤ B*((k*24).choose w:ℝ))
    (hrow : (2:ℝ)^(k*24)*B ≤
      Real.exp ((k*24:ℕ)*(Real.log 2/2+1281/100000+ε)))
    (hα0 : 0 < α) (hα1 : α ≤ 1/10000) (hα : (L:ℝ)*α = Q)
    (hrounds : 128*R = L*(k*24)) (hd : (d:ℝ) ≤ (11/100)*(L:ℝ)*(k*24:ℕ))
    (hschedule : Real.log (L:ℝ) ≤ (4*Real.log 2/39)*(k*24:ℕ)) (hε : ε ≤ 1/500) :
    (∑ x ∈ (Spin.nonzeroMsgs (Fin L → LocalMessage k)).filter (fun x => occupation x = Q),
      failureProbability e (rowSupports seed x) d) ≤
        Real.exp (-(3/500)*(Q:ℝ)*(k*24:ℕ)) := by
  refine (nonzero_occupation_sparse_failure seed Q e d hB hα0 hα1).trans ?_
  apply countedMassBound_le hQ hQL hb hα0 hα1 hα hrounds hd hschedule hε
  · positivity
  · have hh := pow_le_pow_left₀ (by positivity : 0 ≤ (2:ℝ)^(k*24)*B) hrow Q
    rw [← Real.exp_nat_mul] at hh
    convert hh using 1 <;> ring

end Spin.Structured.SparseRate
