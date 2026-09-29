import SpinCodes.Structured.ConcreteMarkedShuffle

noncomputable section
namespace Spin.Structured.ConcreteMarked
open Routing

theorem masked_product_eq_bind {n : ℕ} (P : FinPMF (Finset (Fin n))) :
    (P.prod (iidBits n (1/2) (by norm_num) (by norm_num))).map
      (fun x => x.1 ∩ x.2) = P.bind fairOnMarks := by
  apply FinPMF.eq_of_expect_eq
  intro F
  rw [FinPMF.expect_map, FinPMF.expect_prod, FinPMF.expect_bind]
  simp only [fairOnMarks_eq_mask, FinPMF.expect_map]

/-- The explicit mark/bit experiment uses independent iid marks and fair bits. -/
def markBitLaw (n : ℕ) (p : ℝ) (hp0 : 0 ≤ p) (hp1 : p ≤ 1) :
    FinPMF (Finset (Fin n) × Finset (Fin n)) :=
  (iidBits n p hp0 hp1).prod (iidBits n (1/2) (by norm_num) (by norm_num))

theorem markBitLaw_unconditioned (n : ℕ) (p : ℝ) (hp0 : 0 ≤ p) (hp1 : p ≤ 1) :
    (markBitLaw n p hp0 hp1).map (fun x => x.1 ∩ x.2) =
      iidBits n (p/2) (by positivity) (by linarith) := by
  rw [markBitLaw, masked_product_eq_bind, iid_thinning]

theorem markBitLaw_prob_count (n Q : ℕ) (p : ℝ) (hp0 : 0 ≤ p) (hp1 : p ≤ 1) :
    (markBitLaw n p hp0 hp1).prob (fun x => x.1.card = Q) = markMass n Q p := by
  exact (FinPMF.prob_prod_left (iidBits n p hp0 hp1)
    (iidBits n (1/2) (by norm_num) (by norm_num)) (fun x => x.card = Q)).trans
      (iidBits_prob_card n Q p hp0 hp1)

/-- Conditioning only the marks leaves the independent fair bits intact. -/
theorem markBitLaw_conditioned {n : ℕ} (S : Finset (Fin n))
    (p : ℝ) (hp0 : 0 ≤ p) (hp1 : p ≤ 1)
    (h : 0 < (markBitLaw n p hp0 hp1).prob (fun x => x.1.card = S.card)) :
    ((markBitLaw n p hp0 hp1).condition (fun x => x.1.card = S.card) h).map
      (fun x => x.1 ∩ x.2) = regionLaw S := by
  have hm : 0 < (iidBits n p hp0 hp1).prob (fun x => x.card = S.card) := by
    rw [iidBits_prob_card]
    rwa [markBitLaw_prob_count] at h
  have he := FinPMF.condition_prod_left (iidBits n p hp0 hp1)
    (iidBits n (1/2) (by norm_num) (by norm_num)) (fun x => x.card = S.card) hm h
  change (((iidBits n p hp0 hp1).prod (iidBits n (1/2) (by norm_num) (by norm_num))).condition
    (fun x => x.1.card = S.card) h).map (fun x => x.1 ∩ x.2) = _
  rw [he, masked_product_eq_bind, conditioned_marks_eq_shuffleLaw]
  rfl

end Spin.Structured.ConcreteMarked

