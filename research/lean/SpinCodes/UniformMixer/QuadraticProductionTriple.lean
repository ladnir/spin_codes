import SpinCodes.UniformMixer.QuadraticCoefficientTriple

/-!
# The actual production word is a quadratic coefficient triple

This bridge identifies both setup distributions and output weights:
each production state corresponds to exactly one constant, linear,
and quadratic coefficient triple.
-/

namespace Spin.UniformMixer

noncomputable def productionTriple (q : State22) : QuadraticTriple :=
  coefficientsTriple (productionQuadraticCoefficients q)

theorem productionTriple_bijective :
    Function.Bijective productionTriple := by
  constructor
  · intro q r h
    apply productionQuadraticCoefficients_bijective.1
    exact quadraticTripleEquivCoefficients.symm.injective h
  · intro t
    obtain ⟨c, hc⟩ := quadraticTripleEquivCoefficients.symm.surjective t
    obtain ⟨q, hq⟩ := productionQuadraticCoefficients_bijective.2 c
    refine ⟨q, ?_⟩
    change coefficientsTriple c = t at hc
    simpa [productionTriple, hq] using hc

set_option maxRecDepth 8192 in
theorem productionWeight_as_triple (q : State22) :
    coordinateWeight coordMap q =
      boolWeight6 (fun x =>
        (productionTriple q).1 +
        dot6Linear x (productionTriple q).2.1 +
        quadraticPart (productionTriple q).2.2 x) := by
  rw [productionWeight_as_quadraticWeight]
  apply congrArg boolWeight6
  funext x
  rw [← tripleQuadraticCoefficients_coefficientsTriple
    (productionQuadraticCoefficients q)]
  rw [tripleQuadraticEvaluation_formula]
  have hx : outputIndexState (outputIndexEquiv.symm x) = x :=
    outputIndexEquiv.apply_symm_apply x
  simp [productionTriple, hx]

end Spin.UniformMixer
