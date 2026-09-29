import SpinCodes.Structured.ConcreteReshape
import SpinCodes.Structured.ConcreteRoutePermutation
import SpinCodes.Structured.ConcreteEncoderMoment
import SpinCodes.Structured.ConcreteSerialization

noncomputable section
namespace Spin.Structured.ConcreteRoutedEncoder
open ConcreteRoute ConcreteEncoder Spin.Imt

/-- All randomness: independent row permutations, region permutations, and transvections. -/
def experimentLaw (L b R : ℕ) : FinPMF (Seeds L b × (Fin R → Transvection)) :=
  (seedLaw L b).prod (piPMF (fun _ : Fin R => transvectionLaw))

/-- Emitted weight for a fixed input matrix and a fixed wiring into IMT blocks. -/
def weight {L b R : ℕ} (e : (Fin b × Fin L) ≃ (Fin R × Fin 128))
    (rows : Fin L → Finset (Fin b)) (ω : Seeds L b × (Fin R → Transvection)) : ℕ :=
  outputWeight (reshape e (routeEval rows ω.1)) ω.2 ∅

theorem moment_eq {L b R : ℕ} (e : (Fin b × Fin L) ≃ (Fin R × Fin 128))
    (rows : Fin L → Finset (Fin b)) (z : ℝ) :
    (experimentLaw L b R).expect (fun ω => z ^ weight e rows ω) =
      (seedLaw L b).expect (fun seed => inputMoment z (reshape e (routeEval rows seed)) ∅) := by
  rw [experimentLaw, FinPMF.expect_prod]
  rfl

theorem iid_moment_eq {L b R : ℕ} (e : (Fin b × Fin L) ≃ (Fin R × Fin 128))
    (q z : ℝ) (hq0 : 0 ≤ q) (hq1 : q ≤ 1) :
    (iidLaw L b q hq0 hq1).expect (fun regions => inputMoment z (reshape e regions) ∅) =
      (iidLaw 128 R q hq0 hq1).expect (fun xs => inputMoment z xs ∅) := by
  rw [← FinPMF.expect_map (iidLaw L b q hq0 hq1) (reshape e)
    (fun xs => inputMoment z xs ∅), reshape_iid]

/-- Finite-round matrix bound for the actual routed IMT encoder. -/
theorem moment_bound {L b R : ℕ} (hL : 0 < L) (hb : 0 < b)
    (e : (Fin b × Fin L) ≃ (Fin R × Fin 128))
    (rows : Fin L → Finset (Fin b)) (hq0 : 0 < density rows) (hq1 : density rows < 1)
    {z : ℝ} (hz0 : 0 ≤ z) (hz1 : z ≤ 1) :
    (experimentLaw L b R).expect (fun ω => z ^ weight e rows ω) ≤
      (((b : ℝ) + 1) ^ activeRows rows * ((L : ℝ) + 1) ^ b) *
        (((Occupation.Sparse.numericalMatrix (density rows) z).apply)^[R]
          (Coords.eZ 5)).total := by
  rw [moment_eq]
  refine (permutation_expect_le hL hb rows hq0 hq1
    (fun regions => inputMoment z (reshape e regions) ∅)
    (fun regions => inputMoment_nonneg hz0 _ _)).trans ?_
  rw [iid_moment_eq]
  apply mul_le_mul_of_nonneg_left _ (by positivity)
  exact encoder_moment_bound
    (poissonBinom (fun _ : Fin 128 => density rows) (fun _ => hq0.le) (fun _ => hq1.le))
    (fun x => poissonBinom_const_apply _ _ x) hq0.le hq1.le hz0 hz1 R

/-- Sparse specialization, including the full route cost. -/
theorem sparse_moment_bound {L b R : ℕ} (hL : 0 < L) (hb : 0 < b)
    (e : (Fin b × Fin L) ≃ (Fin R × Fin 128))
    (rows : Fin L → Finset (Fin b)) {α : ℝ}
    (hα0 : 0 < α) (hα1 : α ≤ 1 / 10000) (hq : density rows = (4 / 5) * α) :
    (experimentLaw L b R).expect
        (fun ω => (1 - (8 / 5) * α) ^ weight e rows ω) ≤
      (((b : ℝ) + 1) ^ activeRows rows * ((L : ℝ) + 1) ^ b) *
        (2048 * (1 - 96 * α) ^ R) := by
  have hq0 : 0 < density rows := by rw [hq]; positivity
  have hq1 : density rows < 1 := by rw [hq]; linarith
  have hz0 : 0 ≤ 1 - (8 / 5) * α := by linarith
  rw [moment_eq]
  refine (permutation_expect_le hL hb rows hq0 hq1
    (fun regions => inputMoment (1 - (8 / 5) * α) (reshape e regions) ∅)
    (fun regions => inputMoment_nonneg hz0 _ _)).trans ?_
  rw [iid_moment_eq]
  apply mul_le_mul_of_nonneg_left _ (by positivity)
  apply sparse_encoder_moment
    (poissonBinom (fun _ : Fin 128 => density rows) (fun _ => hq0.le) (fun _ => hq1.le))
    _ hα0 hα1 R
  intro x
  rw [poissonBinom_const_apply, hq]

/-- The paper's region-major stream, with the number of rounds determined by its bit length. -/
theorem stream_moment_bound {L b : ℕ} (hL : 0 < L) (hb : 0 < b)
    (hdiv : 128 ∣ b * L) (rows : Fin L → Finset (Fin b))
    (hw0 : 0 < totalWeight rows) (hw1 : totalWeight rows < L * b)
    {z : ℝ} (hz0 : 0 ≤ z) (hz1 : z ≤ 1) :
    (experimentLaw L b (b * L / 128)).expect
        (fun ω => z ^ weight (streamWiring hdiv) rows ω) ≤
      (((b : ℝ) + 1) ^ activeRows rows * ((L : ℝ) + 1) ^ b) *
        (((Occupation.Sparse.numericalMatrix (density rows) z).apply)^[b * L / 128]
          (Coords.eZ 5)).total := by
  have hq0 : 0 < density rows := by
    unfold density
    exact div_pos (by exact_mod_cast hw0) (by positivity)
  have hq1 : density rows < 1 := by
    unfold density
    apply (div_lt_one (by positivity)).mpr
    exact_mod_cast hw1
  exact moment_bound hL hb (streamWiring hdiv) rows hq0 hq1 hz0 hz1

end Spin.Structured.ConcreteRoutedEncoder
