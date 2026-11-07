import SpinCodes.UniformMixer.Map22CodeSpace

/-!
# The production image is the full quadratic Boolean code

The coefficients of degree at most two form a 22-dimensional binary
space. The production coefficient map is injective between two finite
spaces of the same size, hence bijective. This establishes that the
actual 22-to-64 expansion realizes every quadratic Boolean polynomial
on six variables; the weight enumerator remains to be proved.
-/

noncomputable section

namespace Spin.UniformMixer

def QuadraticCoefficients :=
  {c : Fin 64 → ZMod 2 //
    ∀ m : Fin 64, m ∉ lowDegreeMonomials → c m = 0}

def quadraticCoefficientRestriction :
    QuadraticCoefficients ≃ (lowDegreeMonomials → ZMod 2) where
  toFun c m := c.1 m.val
  invFun f := ⟨fun m => if h : m ∈ lowDegreeMonomials then
      f ⟨m, h⟩ else 0, by
        intro m hm
        simp [hm]⟩
  left_inv c := by
    apply Subtype.ext
    funext m
    by_cases hm : m ∈ lowDegreeMonomials
    · simp [hm]
    · simp [hm, c.2 m hm]
  right_inv f := by
    funext m
    simp

noncomputable instance : Fintype QuadraticCoefficients :=
  Fintype.ofEquiv (lowDegreeMonomials → ZMod 2)
    quadraticCoefficientRestriction.symm

theorem quadraticCoefficients_card :
    Fintype.card QuadraticCoefficients = 2 ^ 22 := by
  classical
  rw [Fintype.card_congr quadraticCoefficientRestriction,
    Fintype.card_fun]
  simp [lowDegreeMonomials_card, ZMod.card]

def productionQuadraticCoefficients (q : State22) :
    QuadraticCoefficients :=
  ⟨productionANFCoeff q, by
    intro m hm
    exact productionANFCoeff_high_zero q m (by
      simp only [lowDegreeMonomials, Finset.mem_filter,
        Finset.mem_univ, true_and] at hm
      omega)⟩

theorem productionQuadraticCoefficients_val (q : State22) :
    (productionQuadraticCoefficients q).1 = productionANFCoeff q := rfl

theorem productionQuadraticCoefficients_bijective :
    Function.Bijective productionQuadraticCoefficients := by
  have hinj : Function.Injective productionQuadraticCoefficients := by
    intro q r h
    have hv := congrArg (fun c : QuadraticCoefficients => c.1) h
    rw [productionQuadraticCoefficients_val,
      productionQuadraticCoefficients_val] at hv
    exact productionANFCoeff_injective hv
  have hcard : Fintype.card State22 =
      Fintype.card QuadraticCoefficients := by
    rw [quadraticCoefficients_card, card_binary_vector]
  exact (Fintype.bijective_iff_injective_and_card
    productionQuadraticCoefficients).mpr ⟨hinj, hcard⟩

def quadraticEvaluation (c : QuadraticCoefficients) (j : Fin 64) :
    ZMod 2 :=
  ∑ m : Fin 64, if isSubmask m j then c.1 m else 0

theorem production_quadratic_evaluation (q : State22) (j : Fin 64) :
    coordMap j q = quadraticEvaluation
      (productionQuadraticCoefficients q) j := by
  rw [coordMap_anf]
  rfl

theorem every_quadratic_is_production (c : QuadraticCoefficients) :
    ∃ q : State22, ∀ j : Fin 64,
      coordMap j q = quadraticEvaluation c j := by
  obtain ⟨q, hq⟩ := productionQuadraticCoefficients_bijective.2 c
  refine ⟨q, fun j => ?_⟩
  rw [production_quadratic_evaluation]
  rw [hq]

end Spin.UniformMixer
