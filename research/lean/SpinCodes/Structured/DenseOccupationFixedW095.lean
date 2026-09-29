import SpinCodes.Structured.DenseOccupationFixedCheck
import SpinCodes.Structured.DenseOccupationFixedW095Data

noncomputable section
namespace Spin.Structured.DenseOccupationFixed.W095
open Spin.Numeric Spin.Imt

def qReal : ℝ := (qn:ℝ)/qd
def zReal : ℝ := (z.lo:ℝ)/scale
def radiusReal : ℝ := (radius.lo:ℝ)/scale
def w : Coords 5 := ⟨(v.Z.lo:ℝ)/scale, (v.D.lo:ℝ)/scale, fun i => ((v.S i).lo:ℝ)/scale⟩

theorem w_mem : v.Mem w := by
  refine ⟨?_,?_,fun i => ?_⟩
  · exact Fix.sc_mem _
  · exact Fix.sc_mem _
  · fin_cases i <;> exact Fix.sc_mem _

theorem collatz : ((Occupation.Sparse.numericalMatrix qReal zReal).applyCol w).le
    (Coords.smul radiusReal w) := by
  exact checked_collatz qn qd (by decide) (powers z zPowers 208) zReal
    (powers_mem (show Fix.Mem z zReal from Fix.sc_mem _) powers_checked)
    w_mem (show Fix.Mem radius radiusReal from Fix.sc_mem _)
    column_Z_checked column_D_checked column_S_checked

theorem parameters : 0<qReal ∧ qReal<1 ∧ 0<zReal ∧ zReal<1 ∧ 0<radiusReal ∧ radiusReal<1 := by
  norm_num [qReal, qn, qd, zReal, z, radiusReal, radius, Fix.sc, scale]

theorem witness_floor : (1:ℝ)/33 ≤ w.Z ∧ (1:ℝ)/33 ≤ w.D ∧
    ∀ i, (1:ℝ)/33 ≤ w.S i := by
  refine ⟨?_,?_,fun i => ?_⟩
  · norm_num [w, v, Fix.sc, scale]
  · norm_num [w, v, Fix.sc, scale]
  · fin_cases i <;> norm_num [w, v, Fix.sc, scale]

/-- The numerical occupation witness 95 is kernel checked for every finite round count. -/
theorem iterate (R : Nat) :
    (((Occupation.Sparse.numericalMatrix qReal zReal).apply)^[R] (Coords.eZ 5)).total ≤
      33*radiusReal^R := by
  have hp := parameters
  have hw := witness_floor
  have h := PositiveFinite.occupation_iterate_le hp.1.le hp.2.1.le hp.2.2.1.le
    (by norm_num : (0:ℝ)<1/33) hw.1 hw.2.1 hw.2.2 collatz R
  simpa [w, v, Fix.sc, scale, mul_comm] using h

/-- End-to-end probability bound for the actual routed encoder at this certified witness. -/
theorem routed_probability {L b R : Nat} (hL : 0<L) (hb : 0<b)
    (e : (Fin b × Fin L) ≃ (Fin R × Fin 128))
    (rows : Fin L → Finset (Fin b)) (hw0 : 0<ConcreteRoute.totalWeight rows)
    (hw1 : ConcreteRoute.totalWeight rows<L*b) (d : Nat) :
    (ConcreteRoutedEncoder.experimentLaw L b R).prob
      (fun ω => ConcreteRoutedEncoder.weight e rows ω≤d) ≤
      ((((b:ℝ)+1)^ConcreteRoute.activeRows rows*((L:ℝ)+1)^b)*
        Real.exp (((L:ℝ)*b)*binKL (ConcreteRoute.density rows) qReal)) *
        (33*radiusReal^R)/zReal^d := by
  have hp := parameters
  have hf := witness_floor
  have h := PositiveFinite.routed_tilted_probability_le hL hb e rows hw0 hw1
    hp.1 hp.2.1 hp.2.2.1 hp.2.2.2.1.le
    (by norm_num : (0:ℝ)<1/33) hf.1 hf.2.1 hf.2.2 collatz d
  simpa [w, v, Fix.sc, scale, mul_comm] using h

#print axioms collatz
#print axioms iterate
#print axioms routed_probability

end Spin.Structured.DenseOccupationFixed.W095
