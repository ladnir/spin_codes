import SpinCodes.UniformMixer.QuadraticRankWalshClasses
import SpinCodes.UniformMixer.Map22QuadraticSpace
import SpinCodes.UniformMixer.QuadraticSpectrumBridge

/-!
# Canonical constant, linear, and quadratic coefficients

The 22 low-degree monomials are indexed by one constant, six variables,
and fifteen variable pairs. This explicit order matches `pairLeft` and
`pairRight` in the polar analysis.
-/

namespace Spin.UniformMixer

def canonicalMonomial : Fin 22 → Fin 64 := ![
  0, 1, 2, 4, 8, 16, 32,
  3, 5, 9, 17, 33, 6, 10, 18, 34,
  12, 20, 36, 24, 40, 48]

theorem canonicalMonomial_lowDegree (i : Fin 22) :
    canonicalMonomial i ∈ lowDegreeMonomials := by
  decide +revert

def canonicalMonomialSubtype (i : Fin 22) : lowDegreeMonomials :=
  ⟨canonicalMonomial i, canonicalMonomial_lowDegree i⟩

theorem canonicalMonomial_injective :
    Function.Injective canonicalMonomial := by
  decide

theorem canonicalMonomialSubtype_bijective :
    Function.Bijective canonicalMonomialSubtype := by
  have hinj : Function.Injective canonicalMonomialSubtype := by
    intro i j h
    exact canonicalMonomial_injective (congrArg Subtype.val h)
  have hcard : Fintype.card (Fin 22) =
      Fintype.card lowDegreeMonomials := by
    simp [lowDegreeMonomials_card]
  exact (Fintype.bijective_iff_injective_and_card _).mpr
    ⟨hinj, hcard⟩

noncomputable def canonicalMonomialEquiv : Fin 22 ≃ lowDegreeMonomials :=
  Equiv.ofBijective canonicalMonomialSubtype
    canonicalMonomialSubtype_bijective

abbrev QuadraticTriple := ZMod 2 × State6 × (Fin 15 → ZMod 2)

def tripleVector (t : QuadraticTriple) : State22 :=
  fun i =>
    if h0 : i.val = 0 then t.1
    else if h6 : i.val ≤ 6 then
      t.2.1 ⟨i.val - 1, by omega⟩
    else t.2.2 ⟨i.val - 7, by omega⟩

def vectorTriple (v : State22) : QuadraticTriple :=
  (v ⟨0, by decide⟩,
    fun i => v ⟨i.val + 1, by omega⟩,
    fun p => v ⟨p.val + 7, by omega⟩)

theorem vectorTriple_tripleVector (t : QuadraticTriple) :
    vectorTriple (tripleVector t) = t := by
  rcases t with ⟨b, a, c⟩
  apply Prod.ext
  · rfl
  · apply Prod.ext
    · funext i
      simp [vectorTriple, tripleVector]
    · funext p
      simp [vectorTriple, tripleVector]

theorem tripleVector_vectorTriple (v : State22) :
    tripleVector (vectorTriple v) = v := by
  funext i
  fin_cases i <;> rfl

def quadraticTripleEquivVector : QuadraticTriple ≃ State22 where
  toFun := tripleVector
  invFun := vectorTriple
  left_inv := vectorTriple_tripleVector
  right_inv := tripleVector_vectorTriple

noncomputable def tripleQuadraticCoefficients
    (t : QuadraticTriple) : QuadraticCoefficients :=
  quadraticCoefficientRestriction.symm
    (fun m => tripleVector t (canonicalMonomialEquiv.symm m))

noncomputable def coefficientsTriple
    (c : QuadraticCoefficients) : QuadraticTriple :=
  vectorTriple (fun i => quadraticCoefficientRestriction c
    (canonicalMonomialEquiv i))

theorem coefficientsTriple_tripleQuadraticCoefficients
    (t : QuadraticTriple) :
    coefficientsTriple (tripleQuadraticCoefficients t) = t := by
  have hv : (fun i => quadraticCoefficientRestriction
      (tripleQuadraticCoefficients t) (canonicalMonomialEquiv i)) =
      tripleVector t := by
    funext i
    simp [tripleQuadraticCoefficients]
  exact (congrArg vectorTriple hv).trans (vectorTriple_tripleVector t)

theorem tripleQuadraticCoefficients_coefficientsTriple
    (c : QuadraticCoefficients) :
    tripleQuadraticCoefficients (coefficientsTriple c) = c := by
  apply quadraticCoefficientRestriction.injective
  funext m
  simp [tripleQuadraticCoefficients, coefficientsTriple,
    tripleVector_vectorTriple]

noncomputable def quadraticTripleEquivCoefficients :
    QuadraticTriple ≃ QuadraticCoefficients where
  toFun := tripleQuadraticCoefficients
  invFun := coefficientsTriple
  left_inv := coefficientsTriple_tripleQuadraticCoefficients
  right_inv := tripleQuadraticCoefficients_coefficientsTriple

theorem quadraticEvaluation_canonical (c : QuadraticCoefficients)
    (j : Fin 64) :
    quadraticEvaluation c j =
      ∑ i : Fin 22,
        if isSubmask (canonicalMonomial i) j then
          quadraticCoefficientRestriction c (canonicalMonomialEquiv i)
        else 0 := by
  have hrestrict : quadraticEvaluation c j =
      ∑ m ∈ lowDegreeMonomials,
        if isSubmask m j then c.1 m else 0 := by
    unfold quadraticEvaluation
    change (∑ m : Fin 64, if isSubmask m j then c.1 m else 0) =
      ∑ m ∈ Finset.univ.filter (fun m => monomialDegree m ≤ 2),
        if isSubmask m j then c.1 m else 0
    rw [Finset.sum_filter]
    apply Finset.sum_congr rfl
    intro m _
    by_cases hm : m ∈ lowDegreeMonomials
    · have hd : monomialDegree m ≤ 2 := by
        simpa [lowDegreeMonomials] using hm
      simp [hd]
    · have hd : ¬ monomialDegree m ≤ 2 := by
        simpa [lowDegreeMonomials] using hm
      simp [hd, c.2 m hm]
  rw [hrestrict]
  rw [← Finset.sum_coe_sort]
  apply Fintype.sum_equiv canonicalMonomialEquiv.symm
  intro m
  have hm : canonicalMonomial (canonicalMonomialEquiv.symm m) = m.1 :=
    congrArg Subtype.val (canonicalMonomialEquiv.apply_symm_apply m)
  rw [hm]
  rw [canonicalMonomialEquiv.apply_symm_apply]
  rfl

theorem tripleQuadraticEvaluation_as_sum (t : QuadraticTriple)
    (j : Fin 64) :
    quadraticEvaluation (tripleQuadraticCoefficients t) j =
      ∑ i : Fin 22,
        if isSubmask (canonicalMonomial i) j then
          tripleVector t i else 0 := by
  rw [quadraticEvaluation_canonical]
  apply Finset.sum_congr rfl
  intro i _
  congr 1
  simp [tripleQuadraticCoefficients]

set_option maxRecDepth 100000 in
set_option maxHeartbeats 0 in
theorem tripleQuadraticEvaluation_formula (t : QuadraticTriple)
    (j : Fin 64) :
    quadraticEvaluation (tripleQuadraticCoefficients t) j =
      t.1 + dot6Linear (outputIndexState j) t.2.1 +
        quadraticPart t.2.2 (outputIndexState j) := by
  rw [tripleQuadraticEvaluation_as_sum]
  fin_cases j <;>
    simp [Fin.sum_univ_succ, canonicalMonomial, tripleVector,
      outputIndexState, dot6Linear, quadraticPart,
      isSubmask, bitScalar, pairLeft, pairRight,
      Nat.testBit_succ] <;> ring

end Spin.UniformMixer
