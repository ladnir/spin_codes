import SpinCodes.Structured.ConcreteMarkedConditioning
import SpinCodes.Structured.ConcreteRoutedMoment

noncomputable section
namespace Spin.Structured.ConcreteMarked
open ConcreteRoute ConcreteEncoder

/-- The fair-row reference consists of independent marked regions. -/
def regionsLaw {L : ℕ} (S : Finset (Fin L)) (b : ℕ) :
    FinPMF (Fin b → Finset (Fin L)) := piPMF (fun _ => regionLaw S)

theorem regionsLaw_dominates {L : ℕ} (S : Finset (Fin L)) (b : ℕ)
    (p : ℝ) (hp0 : 0 ≤ p) (hp1 : p ≤ 1) (hm : 0 < markMass L S.card p) :
    Dominates (regionsLaw S b) (iidLaw L b (p/2) (by positivity) (by linarith))
      ((1 / markMass L S.card p) ^ b) := by
  have hh := dominates_piPMF (fun _ : Fin b => regionLaw S)
    (fun _ => iidBits L (p/2) (by positivity) (by linarith))
    (fun _ => 1 / markMass L S.card p)
    (fun _ => regionLaw_dominates S p hp0 hp1 hm)
  simpa only [regionsLaw, iidBits, iidLaw, Finset.prod_const,
    Finset.card_univ, Fintype.card_fin] using hh

/-- Finite sparse moment after the efficient mark-count conditioning comparison. -/
theorem sparse_moment {L b R : ℕ} (S : Finset (Fin L))
    (e : (Fin b × Fin L) ≃ (Fin R × Fin 128))
    {α : ℝ} (hα0 : 0 < α) (hα1 : α ≤ 1 / 10000) :
    (regionsLaw S b).expect
      (fun regions => inputMoment (1 - (8/5)*α) (reshape e regions) ∅) ≤
        ((1 / markMass L S.card ((8/5)*α)) ^ b) *
          (2048 * (1 - 96*α) ^ R) := by
  have hp0 : 0 < (8/5 : ℝ)*α := by positivity
  have hp1 : (8/5 : ℝ)*α < 1 := by linarith
  have hz : 0 ≤ 1 - (8/5 : ℝ)*α := by linarith
  have hd := regionsLaw_dominates S b ((8/5)*α) hp0.le hp1.le
    (markMass_pos S hp0 hp1)
  refine (hd.expect_le (fun regions => inputMoment_nonneg hz _ _)).trans ?_
  rw [ConcreteRoutedEncoder.iid_moment_eq]
  apply mul_le_mul_of_nonneg_left _
    (pow_nonneg (div_nonneg zero_le_one (markMass_pos S hp0 hp1).le) _)
  apply sparse_encoder_moment (iidBits 128 (((8/5)*α)/2) (by positivity) (by linarith))
    _ hα0 hα1 R
  intro x
  unfold iidBits
  rw [poissonBinom_const_apply]
  congr 2 <;> ring

end Spin.Structured.ConcreteMarked

