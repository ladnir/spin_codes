import SpinCodes.Structured.ConcreteRoute
import SpinCodes.Structured.ConcreteRowLaw
import SpinCodes.Structured.ConcreteShuffleMixture
import SpinCodes.FiniteProductLaw

/-! The concrete two-stage route and its Bernoulli domination. -/
noncomputable section
namespace Spin.Structured.ConcreteRoute

open Finset Routing

/-- Independently shuffle the supplied rows. -/
def rowLaw {L b : ℕ} (rows : Fin L → Finset (Fin b)) :
    FinPMF (Fin L → Finset (Fin b)) := piPMF (fun i => shuffleLaw (rows i))

/-- Independently shuffle every region after the fixed transpose. -/
def regionKernel {L b : ℕ} (regions : Fin b → Finset (Fin L)) :
    FinPMF (Fin b → Finset (Fin L)) := piPMF (fun j => shuffleLaw (regions j))

/-- The structured route: row shuffles, fixed wiring, then region shuffles. -/
def law {L b : ℕ} (rows : Fin L → Finset (Fin b)) :
    FinPMF (Fin b → Finset (Fin L)) := ((rowLaw rows).map transpose).bind regionKernel

/-- The independent-cell comparison law after the region shuffles. -/
def midLaw {L b : ℕ} (rows : Fin L → Finset (Fin b)) :
    FinPMF (Fin b → Finset (Fin L)) :=
  piPMF (fun _ : Fin b => shuffledLaw (poissonBinom (fun i => rowDensity (rows i))
    (fun i => rowDensity_nonneg (rows i)) (fun i => rowDensity_le_one (rows i))))

/-- Independent Bernoulli cells arranged as regions. -/
def iidLaw (L b : ℕ) (q : ℝ) (hq0 : 0 ≤ q) (hq1 : q ≤ 1) :
    FinPMF (Fin b → Finset (Fin L)) :=
  piPMF (fun _ : Fin b => poissonBinom (fun _ => q) (fun _ => hq0) (fun _ => hq1))

theorem rowLaw_dominates {L b : ℕ} (hb : 0 < b) (rows : Fin L → Finset (Fin b)) :
    Dominates (rowLaw rows) (piPMF (fun i => rowBernoulli (rows i)))
      (((b : ℝ) + 1) ^ activeRows rows) := by
  have h := dominates_piPMF (fun i => shuffleLaw (rows i))
    (fun i => rowBernoulli (rows i)) (fun i => rowCost (rows i))
    (fun i => shuffleLaw_dominates_rowBernoulli hb (rows i))
  simpa only [rowCost, prod_rowCost, rowLaw] using h

theorem law_dominates_mid {L b : ℕ} (hb : 0 < b) (rows : Fin L → Finset (Fin b)) :
    Dominates (law rows) (midLaw rows) (((b : ℝ) + 1) ^ activeRows rows) := by
  have h := ((rowLaw_dominates hb rows).map transpose).bind
    (regionKernel (L := L) (b := b))
  change Dominates (law rows)
    (((piPMF (fun i => rowBernoulli (rows i))).map transpose).bind regionKernel) _ at h
  simp only [rowBernoulli, transpose_bernoulli] at h
  rw [show regionKernel (L := L) (b := b) = (fun x => piPMF (fun j => shuffleLaw (x j))) from rfl] at h
  rw [piPMF_bind] at h
  exact h

theorem region_dominates {L b : ℕ} (hL : 0 < L) (rows : Fin L → Finset (Fin b))
    (hq0 : 0 < density rows) (hq1 : density rows < 1) :
    Dominates
      (shuffledLaw (poissonBinom (fun i => rowDensity (rows i))
        (fun i => rowDensity_nonneg (rows i)) (fun i => rowDensity_le_one (rows i))))
      (poissonBinom (fun _ : Fin L => density rows) (fun _ => hq0.le) (fun _ => hq1.le))
      ((L : ℝ) + 1) := by
  apply dominates_region hL _ (fun i => rowDensity (rows i))
    (fun i => rowDensity_nonneg (rows i)) (fun i => rowDensity_le_one (rows i))
    (average_rowDensity rows) hq0 hq1
  · exact shuffledLaw_fiber_uniform _
  · exact shuffledLaw_prob_card _

theorem midLaw_dominates {L b : ℕ} (hL : 0 < L) (rows : Fin L → Finset (Fin b))
    (hq0 : 0 < density rows) (hq1 : density rows < 1) :
    Dominates (midLaw rows) (iidLaw L b (density rows) hq0.le hq1.le)
      (((L : ℝ) + 1) ^ b) := by
  exact dominates_regions _ _ (fun _ => region_dominates hL rows hq0 hq1)

/-- The concrete route has the exact polynomial domination cost from the paper. -/
theorem law_dominates {L b : ℕ} (hL : 0 < L) (hb : 0 < b)
    (rows : Fin L → Finset (Fin b)) (hq0 : 0 < density rows) (hq1 : density rows < 1) :
    Dominates (law rows) (iidLaw L b (density rows) hq0.le hq1.le)
      (((b : ℝ) + 1) ^ activeRows rows * ((L : ℝ) + 1) ^ b) :=
  (law_dominates_mid hb rows).trans (midLaw_dominates hL rows hq0 hq1) (by positivity)

/-- Every nonnegative statistic transfers from the actual route to iid bits. -/
theorem expect_le {L b : ℕ} (hL : 0 < L) (hb : 0 < b)
    (rows : Fin L → Finset (Fin b)) (hq0 : 0 < density rows) (hq1 : density rows < 1)
    (F : (Fin b → Finset (Fin L)) → ℝ) (hF : ∀ x, 0 ≤ F x) :
    (law rows).expect F ≤
      (((b : ℝ) + 1) ^ activeRows rows * ((L : ℝ) + 1) ^ b) *
        (iidLaw L b (density rows) hq0.le hq1.le).expect F :=
  (law_dominates hL hb rows hq0 hq1).expect_le hF

end Spin.Structured.ConcreteRoute

