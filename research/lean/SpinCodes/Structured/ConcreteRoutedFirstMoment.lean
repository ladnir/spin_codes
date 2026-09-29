import SpinCodes.FiniteMoment
import SpinCodes.Structured.ConcreteRoutedMoment

noncomputable section
attribute [local instance] Classical.propDecidable
namespace Spin.Structured.ConcreteRoutedFirstMoment
open Finset ConcreteRoute ConcreteRoutedEncoder Spin.Imt

abbrev InnerSeeds (L b R : ℕ) := Seeds L b × (Fin R → ConcreteEncoder.Transvection)

/-- Attach the proved route and IMT experiment to an arbitrary supplied outer encoder. -/
def setup {M Ω : Type*} [Fintype M] [Fintype Ω] {L b R : ℕ}
    (Pout : FinPMF Ω) (outer : Ω → M → (Fin L → Finset (Fin b)))
    (e : (Fin b × Fin L) ≃ (Fin R × Fin 128)) :
    Setup M (Fin L → Finset (Fin b)) Ω (InnerSeeds L b R) where
  Pout := Pout
  Pin := ConcreteRoutedEncoder.experimentLaw L b R
  outer := outer
  innerWt seed rows := weight e rows seed

def occupation {M Ω : Type*} {L b : ℕ}
    (outer : Ω → M → (Fin L → Finset (Fin b))) (ω : Ω) (x : M) : ℕ :=
  activeRows (outer ω x)

def matrixBound {L b : ℕ} (R d : ℕ) (z : ℝ) (rows : Fin L → Finset (Fin b)) : ℝ :=
  ((((b : ℝ) + 1) ^ activeRows rows * ((L : ℝ) + 1) ^ b) *
    (((Occupation.Sparse.numericalMatrix (density rows) z).apply)^[R]
      (Coords.eZ 5)).total) / z ^ d

theorem low_weight_probability {L b R : ℕ} (hL : 0 < L) (hb : 0 < b)
    (e : (Fin b × Fin L) ≃ (Fin R × Fin 128))
    (rows : Fin L → Finset (Fin b)) (hq0 : 0 < density rows) (hq1 : density rows < 1)
    (d : ℕ) {z : ℝ} (hz0 : 0 < z) (hz1 : z ≤ 1) :
    (ConcreteRoutedEncoder.experimentLaw L b R).prob (fun ω => weight e rows ω ≤ d) ≤
      matrixBound R d z rows :=
  FinPMF.prob_weight_le_of_moment _ _ d hz0 hz1
    (moment_bound hL hb e rows hq0 hq1 hz0.le hz1)

/-- A total envelope; the exceptional all-zero/all-one densities use probability at most one. -/
def envelope {L b : ℕ} (R d : ℕ) (z : ℝ) (rows : Fin L → Finset (Fin b)) : ℝ :=
  if 0 < density rows ∧ density rows < 1 then matrixBound R d z rows else 1

theorem qd_le_envelope {M Ω : Type*} [Fintype M] [Fintype Ω] {L b R : ℕ}
    (Pout : FinPMF Ω) (outer : Ω → M → (Fin L → Finset (Fin b)))
    (hL : 0 < L) (hb : 0 < b) (e : (Fin b × Fin L) ≃ (Fin R × Fin 128))
    (rows : Fin L → Finset (Fin b)) (d : ℕ) {z : ℝ} (hz0 : 0 < z) (hz1 : z ≤ 1) :
    (setup Pout outer e).qd d rows ≤ envelope R d z rows := by
  unfold envelope
  split_ifs with h
  · exact low_weight_probability hL hb e rows h.1 h.2 d hz0 hz1
  · exact FinPMF.prob_le_one _ _

/-- The finite first moment sums the actual inner failure bounds over the supplied outer words. -/
theorem occupation_first_moment {M Ω : Type*} [Fintype M] [DecidableEq M] [Zero M]
    [Fintype Ω] {L b R : ℕ}
    (Pout : FinPMF Ω) (outer : Ω → M → (Fin L → Finset (Fin b)))
    (hL : 0 < L) (hb : 0 < b) (e : (Fin b × Fin L) ≃ (Fin R × Fin 128))
    (d Q : ℕ) {z : ℝ} (hz0 : 0 < z) (hz1 : z ≤ 1) :
    (setup Pout outer e).joint.expect
        (fun ω => ((setup Pout outer e).ZQ d (occupation outer) Q ω : ℝ)) ≤
      Pout.expect (fun ω => ∑ x ∈ (nonzeroMsgs M).filter
        (fun x => activeRows (outer ω x) = Q), envelope R d z (outer ω x)) := by
  exact Setup.expect_ZQ_le (setup Pout outer e) d (occupation outer) Q
    (fun ω x => envelope R d z (outer ω x))
    (fun ω x _ _ => qd_le_envelope Pout outer hL hb e (outer ω x) d hz0 hz1)

end Spin.Structured.ConcreteRoutedFirstMoment
