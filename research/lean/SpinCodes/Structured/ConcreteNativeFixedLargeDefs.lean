import SpinCodes.Structured.ConcreteFixedNumericPin
import SpinCodes.Structured.ConcreteOuterMajorantSpectrum
import SpinCodes.Structured.ConcreteNativeFixedTwoLimit

noncomputable section
namespace Spin.Structured.ConcreteNativeFamily
open Finset Filter ConcreteOuter ConcreteFixedNumeric Placement

abbrev FixedFugacityChoice (Q : ℕ) := Fin Q → {u : ℚ // u∈ConcreteFixedNumeric.fugacities}

def fixedFugacity {Q : ℕ} (u : FixedFugacityChoice Q) (i : Fin Q) : ℝ := u i

def fixedNormCost {Q : ℕ} (u : FixedFugacityChoice Q) : ℝ :=
  Real.exp ((Q:ℝ)*gThree)*∏ i, Spin.mv (127/250) (3/1600) (1/524287) (fixedFugacity u i)

def fixedNormSlack {Q : ℕ} (u : FixedFugacityChoice Q) : ℝ :=
  fixedNormCost u*(Real.exp ((Q:ℝ)/10000)-1)

def FixedLargeNorm (Q : ℕ) : Prop := ∀ u : FixedFugacityChoice Q,
  Spin.RowNormLe (3/1600) (fixedNormCost u)
    (fugacityContinuum ((Q:ℝ)*(133/125)/(262144/524287)) (fixedFugacity u))

theorem fixedFugacity_pos {Q : ℕ} (u : FixedFugacityChoice Q) (i : Fin Q) : 0 < fixedFugacity u i := by
  have h := (fugacities_range (u i).val (u i).property).1
  have hh : (((9/10:ℚ):ℝ)) ≤ (((u i).val:ℚ):ℝ) := Rat.cast_le.mpr h
  norm_num at hh
  change 0 < (((u i).val:ℚ):ℝ)
  linarith

theorem fixedMv_pos {Q : ℕ} (u : FixedFugacityChoice Q) (i : Fin Q) :
    0 < Spin.mv (127/250) (3/1600) (1/524287) (fixedFugacity u i) := by
  unfold Spin.mv
  apply lt_of_lt_of_le _ (le_max_left _ _)
  have := fixedFugacity_pos u i
  positivity

theorem fixedNormCost_pos {Q : ℕ} (u : FixedFugacityChoice Q) : 0 < fixedNormCost u :=
  mul_pos (Real.exp_pos _) (prod_pos (fun i _ => fixedMv_pos u i))

theorem fixedNormSlack_pos {Q : ℕ} (hQ : 0 < Q) (u : FixedFugacityChoice Q) : 0 < fixedNormSlack u := by
  unfold fixedNormSlack
  apply mul_pos (fixedNormCost_pos u)
  apply sub_pos.mpr
  exact Real.one_lt_exp_iff.mpr (by positivity)

theorem fixedNormCost_add_slack {Q : ℕ} (u : FixedFugacityChoice Q) :
    fixedNormCost u+fixedNormSlack u=fixedNormCost u*Real.exp ((Q:ℝ)/10000) := by
  unfold fixedNormSlack
  ring

theorem choose_fixed_fugacities {Q : ℕ} (x : Fin Q → ℝ)
    (hx : ∀ i, x i∈Set.Icc (13/125) (112/125)) :
    ∃ u : FixedFugacityChoice Q, ∀ i, ConcreteFixedNumeric.rate (x i) (fixedFugacity u i) ≤ -(8679/10000000) := by
  classical
  choose us hm hp hr using fun i => ConcreteFixedNumeric.rate_uniform (hx i)
  exact ⟨fun i => ⟨us i,hm i⟩,hr⟩

#print axioms choose_fixed_fugacities
#print axioms fixedNormCost_pos
#print axioms fixedNormSlack_pos
end Spin.Structured.ConcreteNativeFamily

