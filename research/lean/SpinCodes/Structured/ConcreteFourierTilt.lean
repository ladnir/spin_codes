import SpinCodes.Structured.ConcreteProductTilt
import SpinCodes.Structured.ConcreteFourierProbability
import SpinCodes.Structured.ConcreteScalarMoment

noncomputable section
attribute [local instance] Classical.propDecidable
namespace Spin.Structured.ConcreteFourier
open Finset ConcreteMaps ConcreteScalar
open scoped symmDiff

def emit0 {n : ℕ} (Y : Finset (Fin n)) (z : ℝ) (i : Fin n) : ℝ := if i ∈ Y then z else 1
def emit1 {n : ℕ} (Y : Finset (Fin n)) (z : ℝ) (i : Fin n) : ℝ := if i ∈ Y then 1 else z

theorem emit0_nonneg {n : ℕ} (Y : Finset (Fin n)) {z : ℝ} (hz : 0 ≤ z) (i : Fin n) :
    0 ≤ emit0 Y z i := by unfold emit0; split_ifs <;> positivity
theorem emit1_nonneg {n : ℕ} (Y : Finset (Fin n)) {z : ℝ} (hz : 0 ≤ z) (i : Fin n) :
    0 ≤ emit1 Y z i := by unfold emit1; split_ifs <;> positivity

theorem emission_mass_pos {n : ℕ} (Y : Finset (Fin n)) {β z : ℝ}
    (hb0 : 0 < β) (hb1 : β < 1) (hz : 0 < z) (i : Fin n) :
    0 < tiltMass (fun _ => β) (emit0 Y z) (emit1 Y z) i := by
  have hb : 0 < 1-β := by linarith
  unfold tiltMass emit0 emit1
  split_ifs <;> positivity

theorem pow_card_product {n : ℕ} (z : ℝ) (X : Finset (Fin n)) :
    z^X.card = ∏ i : Fin n, if i ∈ X then z else 1 := by
  rw [prod_ite]
  have h : (univ : Finset (Fin n)).filter (fun i => i ∈ X) = X := by ext; simp
  rw [h]
  simp

theorem emission_product {n : ℕ} (z : ℝ) (X Y : Finset (Fin n)) :
    z^(X ∆ Y).card = ∏ i : Fin n, if i ∈ X then emit1 Y z i else emit0 Y z i := by
  rw [pow_card_product]
  apply prod_congr rfl
  intro i _
  by_cases hx : i ∈ X <;> by_cases hy : i ∈ Y <;>
    simp [emit0, emit1, mem_symmDiff, hx, hy]

theorem emission_mass_product (β z : ℝ) (Y : Finset (Fin 128)) :
    (∏ i, tiltMass (fun _ => β) (emit0 Y z) (emit1 Y z) i) = scalarF β z Y.card := by
  have he : (∏ i, tiltMass (fun _ => β) (emit0 Y z) (emit1 Y z) i) =
      ∏ i : Fin 128, if i ∈ Y then β+(1-β)*z else 1-β+β*z := by
    apply prod_congr rfl
    intro i _
    by_cases hi : i ∈ Y <;> simp [tiltMass, emit0, emit1, hi] <;> ring
  rw [he, prod_ite]
  have hyes : (univ : Finset (Fin 128)).filter (fun i => i ∈ Y) = Y := by ext; simp
  have hno : (univ : Finset (Fin 128)).filter (fun i => i ∉ Y) = univ \ Y := by ext; simp
  rw [hyes, hno]
  simp only [prod_const, card_sdiff, inter_univ, card_univ, Fintype.card_fin, scalarF]
  ring

def emissionTilt (β z : ℝ) (Y : Finset (Fin 128)) : Fin 128 → ℝ :=
  tiltProb (fun _ => β) (emit0 Y z) (emit1 Y z)

def fourierCap (p : Fin 128 → ℝ) : ℝ :=
  (∑ q : Finset (Fin 19), ∏ i ∈ CtransposeSet q, |1-2*p i|) / 524288

theorem fourierCap_nonneg (p : Fin 128 → ℝ) : 0 ≤ fourierCap p := by
  unfold fourierCap
  positivity

theorem prob_decidability {Ω : Type*} [Fintype Ω] (P : FinPMF Ω) (E : Ω → Prop)
    (d₁ d₂ : DecidablePred E) : @FinPMF.prob Ω _ P E d₁ = @FinPMF.prob Ω _ P E d₂ :=
  congrArg (fun D => @FinPMF.prob Ω _ P E D) (Subsingleton.elim d₁ d₂)

/-- The actual emission-weighted syndrome mass has the tilted Fourier cap. -/
theorem weighted_syndrome_le {β z : ℝ} (hb0 : 0 < β) (hb1 : β < 1) (hz : 0 < z)
    (Y : Finset (Fin 128)) (s : Finset (Fin 19)) :
    (poissonBinom (fun _ : Fin 128 => β) (fun _ => hb0.le) (fun _ => hb1.le)).expect
      (fun X => z^(X ∆ Y).card * if Cset X = s then 1 else 0) ≤
      scalarF β z Y.card * fourierCap (emissionTilt β z Y) := by
  simp_rw [emission_product]
  have ht := poissonBinom_tilt_event (fun _ => β) (emit0 Y z) (emit1 Y z)
    (fun _ => hb0.le) (fun _ => hb1.le) (emit0_nonneg Y hz.le) (emit1_nonneg Y hz.le)
    (emission_mass_pos Y hb0 hb1 hz) (fun X => Cset X = s)
  have hb := mul_le_mul_of_nonneg_left
    (syndrome_probability_le _ (tiltProb_nonneg _ _ _ (fun _ => hb0.le)
      (emit1_nonneg Y hz.le) (emission_mass_pos Y hb0 hb1 hz))
      (tiltProb_le_one _ _ _ (fun _ => hb1.le) (emit0_nonneg Y hz.le)
        (emission_mass_pos Y hb0 hb1 hz)) s)
    (Finset.prod_nonneg (s := univ) (fun i _ => (emission_mass_pos Y hb0 hb1 hz i).le))
  have hmid : (∏ i, tiltMass (fun _ => β) (emit0 Y z) (emit1 Y z) i) *
      (poissonBinom (tiltProb (fun _ => β) (emit0 Y z) (emit1 Y z))
        (tiltProb_nonneg _ _ _ (fun _ => hb0.le) (emit1_nonneg Y hz.le)
          (emission_mass_pos Y hb0 hb1 hz))
        (tiltProb_le_one _ _ _ (fun _ => hb1.le) (emit0_nonneg Y hz.le)
          (emission_mass_pos Y hb0 hb1 hz))).prob (fun X => Cset X = s) ≤
      scalarF β z Y.card * fourierCap (emissionTilt β z Y) := by
    simpa only [emission_mass_product, fourierCap, emissionTilt] using hb
  rw [prob_decidability _ (fun X : Finset (Fin 128) => Cset X = s) _
    (fun X => Classical.propDecidable (Cset X = s))] at hmid
  have hf := le_trans ht.le hmid
  convert hf using 1
  apply congrArg ((poissonBinom (fun _ : Fin 128 => β) (fun _ => hb0.le) (fun _ => hb1.le)).expect)
  funext X
  by_cases hx : Cset X = s <;> simp only [hx, ite_true, ite_false, mul_one, mul_zero]

end Spin.Structured.ConcreteFourier


