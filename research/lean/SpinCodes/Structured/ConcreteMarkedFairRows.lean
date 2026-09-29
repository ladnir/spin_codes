import SpinCodes.Structured.ConcreteMarkedShuffle
import SpinCodes.Structured.ConcreteMarkedMoment

noncomputable section
namespace Spin.Structured.ConcreteMarked
open Finset Routing ConcreteRoute

theorem support_coins {n : ℕ} (p : Fin n → ℝ)
    (hp0 : ∀ i, 0 ≤ p i) (hp1 : ∀ i, p i ≤ 1) :
    (piPMF (fun i => coin (p i) (hp0 i) (hp1 i))).map support =
      poissonBinom p hp0 hp1 := by
  ext S
  change ((piPMF (fun i => coin (p i) (hp0 i) (hp1 i))).map (supportEquiv n)).p S = _
  rw [FinPMF.map_equiv_apply]
  simp only [piPMF_apply, coin, poissonBinom_eq_prod]
  change (∏ i : Fin n, if decide (i ∈ S) = true then p i else 1-p i) = _
  simp

theorem coin_mask (mark : Bool) :
    (coin (1/2) (by norm_num) (by norm_num)).map (fun bit => mark && bit) =
      coin (if mark then 1/2 else 0) (by cases mark <;> norm_num)
        (by cases mark <;> norm_num) := by
  ext bit
  cases mark <;> cases bit <;> norm_num [coin, FinPMF.map_p]

def activeProbability {L : ℕ} (S : Finset (Fin L)) (i : Fin L) : ℝ :=
  if i ∈ S then 1/2 else 0

theorem activeProbability_nonneg {L : ℕ} (S : Finset (Fin L)) (i : Fin L) :
    0 ≤ activeProbability S i := by unfold activeProbability; split_ifs <;> norm_num

theorem activeProbability_le_one {L : ℕ} (S : Finset (Fin L)) (i : Fin L) :
    activeProbability S i ≤ 1 := by unfold activeProbability; split_ifs <;> norm_num

theorem fairOnMarks_eq_poisson {L : ℕ} (S : Finset (Fin L)) :
    fairOnMarks S = poissonBinom (activeProbability S)
      (activeProbability_nonneg S) (activeProbability_le_one S) := by
  have hs : support (fun i => decide (i ∈ S)) = S := by ext i; simp
  rw [← hs, fairOnMarks_support]
  simp only [coin_mask]
  rw [support_coins]
  congr 1
  funext i
  simp [activeProbability]

theorem shuffled_iidBits {n : ℕ} (p : ℝ) (hp0 : 0 ≤ p) (hp1 : p ≤ 1) :
    shuffledLaw (iidBits n p hp0 hp1) = iidBits n p hp0 hp1 := by
  apply FinPMF.eq_of_expect_eq
  intro F
  rw [shuffledLaw_expect]
  simp only [iidBits_shuffle]
  exact FinPMF.expect_const _ _

/-- At the chosen active rows, supply independent fair bits; all other rows are zero. -/
def fairRows {L : ℕ} (S : Finset (Fin L)) (b : ℕ) :
    FinPMF (Fin L → Finset (Fin b)) :=
  piPMF (fun i => iidBits b (activeProbability S i)
    (activeProbability_nonneg S i) (activeProbability_le_one S i))

theorem fairRows_rowShuffle {L : ℕ} (S : Finset (Fin L)) (b : ℕ) :
    (fairRows S b).bind rowLaw = fairRows S b := by
  unfold fairRows rowLaw
  rw [piPMF_bind]
  have hi (p : ℝ) (hp0 : 0 ≤ p) (hp1 : p ≤ 1) :
      (iidBits b p hp0 hp1).bind shuffleLaw = iidBits b p hp0 hp1 :=
    shuffled_iidBits p hp0 hp1
  simp only [hi]

theorem fairRows_transpose {L : ℕ} (S : Finset (Fin L)) (b : ℕ) :
    (fairRows S b).map transpose = piPMF (fun _ : Fin b => fairOnMarks S) := by
  unfold fairRows iidBits
  rw [transpose_bernoulli]
  simp only [fairOnMarks_eq_poisson]

/-- Applying the actual two-stage route to fair active rows gives independent marked regions. -/
theorem fairRows_route {L : ℕ} (S : Finset (Fin L)) (b : ℕ) :
    (fairRows S b).bind law = regionsLaw S b := by
  have he : (fairRows S b).bind law =
      (((fairRows S b).bind rowLaw).map transpose).bind regionKernel := by
    rw [FinPMF.map_bind, FinPMF.bind_assoc]
    rfl
  rw [he, fairRows_rowShuffle, fairRows_transpose]
  unfold regionKernel
  rw [piPMF_bind]
  change (piPMF (fun _ : Fin b => shuffledLaw (fairOnMarks S))) = regionsLaw S b
  simp only [shuffled_fairOnMarks, regionsLaw]

end Spin.Structured.ConcreteMarked

