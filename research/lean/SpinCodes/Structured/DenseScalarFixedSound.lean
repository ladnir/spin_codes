import SpinCodes.Structured.DenseScalarFixedDefs
import SpinCodes.Numeric.Fixed
import SpinCodes.Structured.ConcreteScalarRouted

noncomputable section
namespace Spin.Structured.DenseScalarFixed
open Spin.Numeric ConcreteScalar ConcreteMaps Finset

theorem fEval_mem {b z : Fix} {β ζ : ℝ} (hb : Fix.Mem b β) (hz : Fix.Mem z ζ) (d : ℕ) :
    Fix.Mem (fEval b z d) (scalarF β ζ d) :=
  Fix.mul_mem
    (Fix.pow_mem (Fix.add_mem (Fix.sub_mem (by simpa using Fix.ofInt_mem 1) hb) (Fix.mul_mem hb hz)) _)
    (Fix.pow_mem (Fix.add_mem hb (Fix.mul_mem (Fix.sub_mem (by simpa using Fix.ofInt_mem 1) hb) hz)) _)

/-- Six kernel integer inequalities certify the actual scalar transfer bound. -/
theorem checked_scalar {b z : Fix} {β ζ : ℝ} (hb : Fix.Mem b β) (hz : Fix.Mem z ζ)
    {r : Int} (h0 : (fEval b z 0).hi ≤ r)
    (hs : ∀ i : Fin 5, (fEval b z (shellWeight i)).hi ≤ r) :
    scalarBound β ζ ≤ (r : ℝ)/scale := by
  apply max_le
  · exact (fEval_mem hb hz 0).2.trans
      (div_le_div_of_nonneg_right (by exact_mod_cast h0) scaleR_pos.le)
  · apply (Finset.sup'_le_iff _ _).mpr
    intro i _
    exact (fEval_mem hb hz (shellWeight i)).2.trans
      (div_le_div_of_nonneg_right (by exact_mod_cast hs i) scaleR_pos.le)

/-- An actual routed probability bound from the scalar checker, with every finite route cost. -/
theorem checked_routed_probability {L b R : ℕ} (hL : 0 < L) (hb : 0 < b)
    (e : (Fin b × Fin L) ≃ (Fin R × Fin 128)) (rows : Fin L → Finset (Fin b))
    (hw0 : 0 < ConcreteRoute.totalWeight rows) (hw1 : ConcreteRoute.totalWeight rows < L*b)
    {p ζ : ℝ} (hp0 : 0 < p) (hp1 : p < 1) (hz : 0 < ζ) (hz1 : ζ ≤ 1)
    {pi zi : Fix} (hpi : Fix.Mem pi p) (hzi : Fix.Mem zi ζ)
    {r : Int} (h0 : (fEval pi zi 0).hi ≤ r)
    (hs : ∀ i : Fin 5, (fEval pi zi (shellWeight i)).hi ≤ r) (d : ℕ) :
    (ConcreteRoutedEncoder.experimentLaw L b R).prob (fun ω => ConcreteRoutedEncoder.weight e rows ω ≤ d) ≤
      ((((b : ℝ)+1)^ConcreteRoute.activeRows rows * ((L : ℝ)+1)^b) *
        Real.exp (((L : ℝ)*b)*binKL (ConcreteRoute.density rows) p)) * ((r : ℝ)/scale)^R / ζ^d := by
  apply (routed_scalar_probability hL hb e rows hw0 hw1 hp0 hp1 hz hz1 d).trans
  apply div_le_div_of_nonneg_right _ (pow_nonneg hz.le _)
  apply mul_le_mul_of_nonneg_left _ (by positivity)
  exact pow_le_pow_left₀ (scalarBound_nonneg hp0.le hp1.le hz.le) (checked_scalar hpi hzi h0 hs) R

#print axioms checked_scalar
#print axioms checked_routed_probability

end Spin.Structured.DenseScalarFixed
