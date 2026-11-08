import SpinCodes.Turbo.PhysicalLiveRow22

/-!
# Sound symbolic 11-coordinate Turbo physical operator

This matrix has a complete zero-source birth row and a complete
live-source row. It is defined by exact finite expectations and maxima,
not by the retained numerical coefficients. The latter require a
separate checked domination bridge.
-/

noncomputable section

namespace Spin.Turbo

def turboSymbolicRow22 (j : Fin 17) (z : ℝ) :
    Matrix (Fin 11) (Fin 11) ℝ :=
  Matrix.of fun i k =>
    if i = 0 then turboBirthClassMass22 j z k
    else turboLiveRow22 j z i k

theorem turboSymbolicRow22_sound (j : Fin 17)
    (z : ℝ) (hz : 0 ≤ z) :
    physicalRowSound22
      (Matrix.of (turboSelectedKernel22 j z))
      (turboSymbolicRow22 j z) := by
  intro i b μ hμ
  by_cases hi : i = 0
  · subst i
    simpa only [turboSymbolicRow22, Matrix.of_apply,
      ite_true] using
      turboSelectedKernel22_zero_row_sound j z hz b μ hμ
  · simpa only [turboSymbolicRow22, Matrix.of_apply,
      hi, ite_false, turboLiveOutput22] using
      turboSelectedKernel22_live_row_sound j z hz i hi b μ hμ

#print axioms turboSymbolicRow22_sound

end Spin.Turbo
