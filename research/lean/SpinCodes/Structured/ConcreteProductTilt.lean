import SpinCodes.Structured.ConcreteFourierProduct

noncomputable section
attribute [local instance] Classical.propDecidable
namespace Spin.Structured.ConcreteFourier
open Finset ConcreteRoute

def tiltMass {n : ℕ} (p f0 f1 : Fin n → ℝ) (i : Fin n) : ℝ :=
  (1-p i)*f0 i+p i*f1 i

def tiltProb {n : ℕ} (p f0 f1 : Fin n → ℝ) (i : Fin n) : ℝ :=
  p i*f1 i/tiltMass p f0 f1 i

theorem tiltProb_nonneg {n : ℕ} (p f0 f1 : Fin n → ℝ) (hp0 : ∀ i, 0 ≤ p i)
    (hf1 : ∀ i, 0 ≤ f1 i) (hg : ∀ i, 0 < tiltMass p f0 f1 i) (i : Fin n) :
    0 ≤ tiltProb p f0 f1 i := div_nonneg (mul_nonneg (hp0 i) (hf1 i)) (hg i).le

theorem tiltProb_le_one {n : ℕ} (p f0 f1 : Fin n → ℝ) (hp1 : ∀ i, p i ≤ 1)
    (hf0 : ∀ i, 0 ≤ f0 i) (hg : ∀ i, 0 < tiltMass p f0 f1 i) (i : Fin n) :
    tiltProb p f0 f1 i ≤ 1 := by
  apply (div_le_one (hg i)).mpr
  have h := mul_nonneg (sub_nonneg.mpr (hp1 i)) (hf0 i)
  unfold tiltMass
  linarith

/-- Tilting a product law by coordinate factors keeps coordinates independent. -/
theorem poissonBinom_tilt {n : ℕ} (p f0 f1 : Fin n → ℝ)
    (hp0 : ∀ i, 0 ≤ p i) (hp1 : ∀ i, p i ≤ 1)
    (hf0 : ∀ i, 0 ≤ f0 i) (hf1 : ∀ i, 0 ≤ f1 i)
    (hg : ∀ i, 0 < tiltMass p f0 f1 i) (X : Finset (Fin n)) :
    (poissonBinom p hp0 hp1).p X * (∏ i, if i ∈ X then f1 i else f0 i) =
      (∏ i, tiltMass p f0 f1 i) *
        (poissonBinom (tiltProb p f0 f1) (tiltProb_nonneg p f0 f1 hp0 hf1 hg)
          (tiltProb_le_one p f0 f1 hp1 hf0 hg)).p X := by
  simp only [poissonBinom_eq_prod, ← prod_mul_distrib]
  apply prod_congr rfl
  intro i _
  by_cases hi : i ∈ X
  · simp only [hi, ite_true, tiltProb]
    field_simp [(hg i).ne']
  · simp only [hi, ite_false, tiltProb]
    field_simp [(hg i).ne']
    unfold tiltMass
    ring

theorem poissonBinom_tilt_event {n : ℕ} (p f0 f1 : Fin n → ℝ)
    (hp0 : ∀ i, 0 ≤ p i) (hp1 : ∀ i, p i ≤ 1)
    (hf0 : ∀ i, 0 ≤ f0 i) (hf1 : ∀ i, 0 ≤ f1 i)
    (hg : ∀ i, 0 < tiltMass p f0 f1 i) (E : Finset (Fin n) → Prop) :
    (poissonBinom p hp0 hp1).expect
      (fun X => (∏ i, if i ∈ X then f1 i else f0 i) * if E X then 1 else 0) =
      (∏ i, tiltMass p f0 f1 i) *
        (poissonBinom (tiltProb p f0 f1) (tiltProb_nonneg p f0 f1 hp0 hf1 hg)
          (tiltProb_le_one p f0 f1 hp1 hf0 hg)).prob E := by
  rw [FinPMF.prob_eq_expect_indicator]
  simp only [FinPMF.expect, ← mul_assoc, poissonBinom_tilt p f0 f1 hp0 hp1 hf0 hf1 hg,
    Finset.mul_sum]

end Spin.Structured.ConcreteFourier
