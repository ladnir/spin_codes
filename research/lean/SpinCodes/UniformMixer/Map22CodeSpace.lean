import SpinCodes.UniformMixer.Map22Polynomial
import SpinCodes.UniformMixer.Map22Rank

/-!
# Coefficient space of the production expansion

There are exactly 22 Boolean monomials of degree at most two in six
variables. The production map gives distinct coefficient vectors for
distinct states. This is the finite-dimensional identification needed
before applying Reed--Muller weight arguments.
-/

noncomputable section

namespace Spin.UniformMixer

def lowDegreeMonomials : Finset (Fin 64) :=
  Finset.univ.filter (fun m => monomialDegree m ≤ 2)

theorem lowDegreeMonomials_card : lowDegreeMonomials.card = 22 := by
  decide

theorem productionANFCoeff_injective :
    Function.Injective productionANFCoeff := by
  intro q r h
  apply expansion_injective
  funext j
  change coordMap j q = coordMap j r
  rw [coordMap_anf q j, coordMap_anf r j]
  apply Finset.sum_congr rfl
  intro m _
  have hm : productionANFCoeff q m = productionANFCoeff r m :=
    congrFun h m
  rw [hm]

end Spin.UniformMixer
