import SpinCodes.Structured.ConcreteRoutedMoment
import SpinCodes.Structured.Collatz
import SpinCodes.FiniteMoment

noncomputable section
namespace Spin.Structured.PositiveFinite
open Spin.Imt ConcreteRoute ConcreteRoutedEncoder

/-- A checked positive column contracts the concrete seven-state occupation matrix. -/
theorem occupation_iterate_le {q z lam wmin : ℝ} {w : Coords 5}
    (hq0 : 0 ≤ q) (hq1 : q ≤ 1) (hz0 : 0 ≤ z)
    (hwmin : 0 < wmin) (hwZ : wmin ≤ w.Z) (hwD : wmin ≤ w.D)
    (hwS : ∀ i, wmin ≤ w.S i)
    (hw : ((Occupation.Sparse.numericalMatrix q z).applyCol w).le (Coords.smul lam w))
    (R : ℕ) :
    (((Occupation.Sparse.numericalMatrix q z).apply)^[R] (Coords.eZ 5)).total ≤
      lam ^ R * w.Z / wmin := by
  apply (le_div_iff₀ hwmin).mpr
  exact collatz_total (Occupation.Sparse.numericalMatrix_nonneg hq0 hq1 hz0)
    hw hwmin hwZ hwD hwS R

/-- The actual independent permutations and transvections obey the finite Collatz bound. -/
theorem routed_moment_le {L b R : ℕ} (hL : 0 < L) (hb : 0 < b)
    (e : (Fin b × Fin L) ≃ (Fin R × Fin 128))
    (rows : Fin L → Finset (Fin b))
    (hq0 : 0 < density rows) (hq1 : density rows < 1)
    {z lam wmin : ℝ} {w : Coords 5} (hz0 : 0 ≤ z) (hz1 : z ≤ 1)
    (hwmin : 0 < wmin) (hwZ : wmin ≤ w.Z) (hwD : wmin ≤ w.D)
    (hwS : ∀ i, wmin ≤ w.S i)
    (hw : ((Occupation.Sparse.numericalMatrix (density rows) z).applyCol w).le
      (Coords.smul lam w)) :
    (experimentLaw L b R).expect (fun ω => z ^ weight e rows ω) ≤
      (((b:ℝ)+1)^activeRows rows * ((L:ℝ)+1)^b) * (lam^R * w.Z / wmin) :=
  (moment_bound hL hb e rows hq0 hq1 hz0 hz1).trans
    (mul_le_mul_of_nonneg_left
      (occupation_iterate_le hq0.le hq1.le hz0 hwmin hwZ hwD hwS hw R) (by positivity))

/-- Finite low-weight probability, retaining every route and eigenvector prefactor. -/
theorem routed_probability_le {L b R : ℕ} (hL : 0 < L) (hb : 0 < b)
    (e : (Fin b × Fin L) ≃ (Fin R × Fin 128))
    (rows : Fin L → Finset (Fin b))
    (hq0 : 0 < density rows) (hq1 : density rows < 1)
    {z lam wmin : ℝ} {w : Coords 5} (hz0 : 0 < z) (hz1 : z ≤ 1)
    (hwmin : 0 < wmin) (hwZ : wmin ≤ w.Z) (hwD : wmin ≤ w.D)
    (hwS : ∀ i, wmin ≤ w.S i)
    (hw : ((Occupation.Sparse.numericalMatrix (density rows) z).applyCol w).le
      (Coords.smul lam w)) (d : ℕ) :
    (experimentLaw L b R).prob (fun ω => weight e rows ω ≤ d) ≤
      (((b:ℝ)+1)^activeRows rows * ((L:ℝ)+1)^b) * (lam^R * w.Z / wmin) / z^d :=
  FinPMF.prob_weight_le_of_moment _ _ d hz0 hz1
    (routed_moment_le hL hb e rows hq0 hq1 hz0.le hz1 hwmin hwZ hwD hwS hw)

end Spin.Structured.PositiveFinite

