import SpinCodes.UniformMixer.Map22ANFReconstruct

/-!
# The production expansion has a quadratic Boolean representation

Each output bit of the actual 22-to-64 map is the evaluation of an
algebraic normal form of degree at most two. The coefficients are
defined from the checked production row masks. This identifies the
map's structure without replacing its sampled or shared setup law.
-/

noncomputable section

namespace Spin.UniformMixer

def productionANFCoeff (q : State22) (m : Fin 64) : ZMod 2 :=
  ∑ i : Fin 22, anfCoefficient (rowMask i) m * q i

theorem productionANFCoeff_high_zero (q : State22) (m : Fin 64)
    (hm : 2 < monomialDegree m) :
    productionANFCoeff q m = 0 := by
  unfold productionANFCoeff
  apply Finset.sum_eq_zero
  intro i _
  rw [rowMask_degree_le_two i m hm]
  simp

theorem coordMap_anf (q : State22) (j : Fin 64) :
    coordMap j q =
      ∑ m : Fin 64,
        if isSubmask m j then productionANFCoeff q m else 0 := by
  have hrow (i : Fin 22) :
      coordCoeff j i =
        ∑ m : Fin 64,
          if isSubmask m j then anfCoefficient (rowMask i) m else 0 := by
    change truthBit (rowMask i) j = _
    exact (anf_reconstruct (rowMask i) j).symm
  change (∑ i : Fin 22, coordCoeff j i * q i) = _
  simp_rw [hrow]
  simp only [Finset.sum_mul]
  rw [Finset.sum_comm]
  apply Finset.sum_congr rfl
  intro m _
  by_cases hm : isSubmask m j
  · simp [productionANFCoeff]
  · simp [hm]

theorem coordMap_quadratic (q : State22) (j : Fin 64) :
    coordMap j q =
      ∑ m : Fin 64,
        if isSubmask m j ∧ monomialDegree m ≤ 2 then
          productionANFCoeff q m else 0 := by
  rw [coordMap_anf]
  apply Finset.sum_congr rfl
  intro m _
  by_cases hm : monomialDegree m ≤ 2
  · simp [hm]
  · have hzero := productionANFCoeff_high_zero q m (by omega)
    simp [hm, hzero]

end Spin.UniformMixer
