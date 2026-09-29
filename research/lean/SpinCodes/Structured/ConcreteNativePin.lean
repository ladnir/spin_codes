import SpinCodes.Structured.ConcreteNativeSelection

open Spin.Structured Spin.Structured.ConcreteNativeFamily

example (m : ℕ) : threshold m = ⌊(0.11 : ℝ) * Nsched m⌋₊ := threshold_floor m

example (m : ℕ) {seed : ConcreteOuter.NativeSeed m} (hg : nativeGood m seed)
    {message : Fin (Nsched m / 2) → ZMod 2} (hm : message ≠ 0) :
    0 < ConcreteRoute.density (rows m seed message) ∧
      ConcreteRoute.density (rows m seed message) < 1 := native_good_density m hg hm

#print axioms Spin.FinPMF.prob_weight_le
#print axioms Spin.FinPMF.expect_card_filter
#print axioms Spin.Setup.expect_ZQ_eq
#print axioms Spin.Setup.expect_ZQ_condition_le
#print axioms Spin.Structured.ConcreteRoutedFirstMoment.low_weight_probability
#print axioms Spin.Structured.ConcreteRoutedFirstMoment.occupation_first_moment
#print axioms Spin.Structured.ConcreteBinaryEncoding.wordToRows_rowsToWord
#print axioms Spin.Structured.ConcreteNativeFamily.occ_mem
#print axioms Spin.Structured.ConcreteNativeFamily.setup_inner_weight
#print axioms Spin.Structured.ConcreteNativeFamily.threshold_floor
#print axioms Spin.Structured.ConcreteNativeFamily.family
#print axioms Spin.Structured.ConcreteNativeFamily.qd_le_matrixBound
#print axioms Spin.Structured.ConcreteNativeFamily.EZ_le_matrix_sum
#print axioms Spin.Structured.ConcreteNativeFamily.native_good_density
#print axioms Spin.Structured.ConcreteNativeFamily.selected_EZ_le_matrix_sum
