import SpinCodes.Structured.ConcreteMarkedBits
import SpinCodes.Structured.ConcreteShufflePermutation

noncomputable section
namespace Spin.Structured.ConcreteMarked
open Finset Routing

def markMass (n Q : ℕ) (p : ℝ) : ℝ := (n.choose Q : ℝ) * p ^ Q * (1-p) ^ (n-Q)

theorem iidBits_prob_card (n Q : ℕ) (p : ℝ) (hp0 : 0 ≤ p) (hp1 : p ≤ 1) :
    (iidBits n p hp0 hp1).prob (fun S => S.card = Q) = markMass n Q p :=
  poissonBinom_const_prob Q _ _

theorem markMass_pos {n : ℕ} (S : Finset (Fin n)) {p : ℝ}
    (hp0 : 0 < p) (hp1 : p < 1) : 0 < markMass n S.card p := by
  unfold markMass
  have hc := choose_card_pos S
  have hpc : 0 < 1-p := by linarith
  positivity

/-- Conditioning independent marks on their count gives a uniform subset of that size. -/
theorem conditioned_marks_eq_shuffleLaw {n : ℕ} (S : Finset (Fin n))
    (p : ℝ) (hp0 : 0 ≤ p) (hp1 : p ≤ 1)
    (h : 0 < (iidBits n p hp0 hp1).prob (fun T => T.card = S.card)) :
    (iidBits n p hp0 hp1).condition (fun T => T.card = S.card) h = shuffleLaw S := by
  classical
  ext T
  rw [FinPMF.condition_p_apply, shuffleLaw_apply]
  by_cases ht : T.card = S.card
  · rw [if_pos ht, if_pos ht]
    apply (div_eq_div_iff h.ne' (choose_card_pos S).ne').mpr
    rw [iidBits_prob_card]
    unfold iidBits markMass
    rw [poissonBinom_const_apply, ht]
    ring
  · rw [if_neg ht, if_neg ht]

/-- A uniform set of marks carrying independent fair bits. -/
def regionLaw {n : ℕ} (S : Finset (Fin n)) : FinPMF (Finset (Fin n)) :=
  (shuffleLaw S).bind fairOnMarks

theorem regionLaw_eq_permutation {n : ℕ} (S : Finset (Fin n)) :
    regionLaw S = (FinPMF.uniform (Equiv.Perm (Fin n))).bind
      (fun σ => fairOnMarks (shuffleSupport σ S)) := by
  rw [regionLaw, ← uniform_permutation_shuffleLaw S, FinPMF.bind_map]

/-- The marked-region experiment is exactly iid marks conditioned on their count. -/
theorem regionLaw_eq_conditioned {n : ℕ} (S : Finset (Fin n))
    (p : ℝ) (hp0 : 0 ≤ p) (hp1 : p ≤ 1)
    (h : 0 < (iidBits n p hp0 hp1).prob (fun T => T.card = S.card)) :
    regionLaw S = ((iidBits n p hp0 hp1).condition
      (fun T => T.card = S.card) h).bind fairOnMarks := by
  rw [conditioned_marks_eq_shuffleLaw]
  rfl

/-- Removing the mark-count conditioning costs exactly its reciprocal probability. -/
theorem regionLaw_dominates {n : ℕ} (S : Finset (Fin n)) (p : ℝ)
    (hp0 : 0 ≤ p) (hp1 : p ≤ 1) (hm : 0 < markMass n S.card p) :
    Dominates (regionLaw S) (iidBits n (p/2) (by positivity) (by linarith))
      (1 / markMass n S.card p) := by
  have h : 0 < (iidBits n p hp0 hp1).prob (fun T => T.card = S.card) := by
    rwa [iidBits_prob_card]
  rw [regionLaw_eq_conditioned S p hp0 hp1 h]
  have hd := (dominates_condition (iidBits n p hp0 hp1)
    (fun T => T.card = S.card) h).bind fairOnMarks
  simpa only [iid_thinning, iidBits_prob_card] using hd

end Spin.Structured.ConcreteMarked
