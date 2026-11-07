import SpinCodes.UniformMixer.QuadraticSpectrumAggregate

/-!
# Transfer of the quadratic spectrum to the actual production map

The coefficient triple is only an analysis coordinate system. This file
reindexes the count by the proved production-state equivalence, so the
result applies to the sampled map used by the encoder.
-/

namespace Spin.UniformMixer

private def tripleEnumerationEquiv :
    QuadraticTriple ≃ (Fin 15 → ZMod 2) × State6 × ZMod 2 where
  toFun t := (t.2.2, t.2.1, t.1)
  invFun t := (t.2.2, t.2.1, t.1)
  left_inv t := by cases t with | mk b tail => cases tail; rfl
  right_inv t := by cases t with | mk c tail => cases tail; rfl

private noncomputable def productionSpectrumEquiv :
    State22 ≃ (Fin 15 → ZMod 2) × State6 × ZMod 2 :=
  (Equiv.ofBijective productionTriple productionTriple_bijective).trans
    tripleEnumerationEquiv

noncomputable def productionWeightCount (w : Nat) : Nat :=
  ∑ q : State22, if coordinateWeight coordMap q = w then 1 else 0

theorem productionWeightCount_eq_card_filter (w : Nat) :
    productionWeightCount w =
      (Finset.univ.filter (fun q : State22 =>
        coordinateWeight coordMap q = w)).card := by
  unfold productionWeightCount
  rw [← Finset.sum_filter]
  simp

theorem productionWeightCount_eq_quadratic (w : Nat) :
    productionWeightCount w = quadraticTotalWeightCount w := by
  unfold productionWeightCount quadraticTotalWeightCount
  calc
    (∑ q : State22, if coordinateWeight coordMap q = w then 1 else 0) =
        ∑ t : (Fin 15 → ZMod 2) × State6 × ZMod 2,
          if quadraticWordWeight t.2.2 t.2.1 t.1 = w then 1 else 0 := by
      apply Fintype.sum_equiv productionSpectrumEquiv
      intro q
      simp only [productionSpectrumEquiv, Equiv.trans_apply,
        tripleEnumerationEquiv, Equiv.ofBijective_apply]
      rw [productionWeight_as_triple]
      rfl
    _ = ∑ c : Fin 15 → ZMod 2, quadraticPairLowCount c w := by
      rw [Fintype.sum_prod_type]
      apply Finset.sum_congr rfl
      intro c _
      rw [Fintype.sum_prod_type]
      unfold quadraticPairLowCount
      apply Finset.sum_congr rfl
      intro a _
      have hz : (Finset.univ : Finset (ZMod 2)) = {0, 1} := by decide
      have hne : (0 : ZMod 2) ∉ ({1} : Finset (ZMod 2)) := by decide
      rw [hz, Finset.sum_insert hne, Finset.sum_singleton]

theorem productionWeightCount_0 : productionWeightCount 0 = 1 := by
  rw [productionWeightCount_eq_quadratic, quadraticTotalWeightCount_0]

theorem productionWeightCount_16 : productionWeightCount 16 = 2604 := by
  rw [productionWeightCount_eq_quadratic, quadraticTotalWeightCount_16]

theorem productionWeightCount_24 : productionWeightCount 24 = 291648 := by
  rw [productionWeightCount_eq_quadratic, quadraticTotalWeightCount_24]

theorem productionWeightCount_28 : productionWeightCount 28 = 888832 := by
  rw [productionWeightCount_eq_quadratic, quadraticTotalWeightCount_28]

theorem productionWeightCount_36 : productionWeightCount 36 = 888832 := by
  rw [productionWeightCount_eq_quadratic, quadraticTotalWeightCount_36]

theorem productionWeightCount_40 : productionWeightCount 40 = 291648 := by
  rw [productionWeightCount_eq_quadratic, quadraticTotalWeightCount_40]

theorem productionWeightCount_48 : productionWeightCount 48 = 2604 := by
  rw [productionWeightCount_eq_quadratic, quadraticTotalWeightCount_48]

theorem productionWeightCount_64 : productionWeightCount 64 = 1 := by
  rw [productionWeightCount_eq_quadratic, quadraticTotalWeightCount_64]

end Spin.UniformMixer
