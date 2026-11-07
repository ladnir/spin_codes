import SpinCodes.UniformMixer.QuadraticRankFinal

/-!
# Fast rank equals the actual quadratic radical

This connects the executable row elimination to the mathematical radical
of every one of the 32,768 packed quadratic forms. Subsequent finite
histogram checks may use the rank evaluator without changing their meaning.
-/

namespace Spin.UniformMixer

theorem solutionCount_polarRows_eq_radical (q : Fin (2 ^ 15)) :
    solutionCount (polarRows q.val) =
      (quadraticRadical (packedQuadratic q)).card := by
  classical
  unfold solutionCount
  congr 1
  ext x
  simp only [Finset.mem_filter, Finset.mem_univ, true_and]
  exact polarRows_vanish_iff q x

theorem packedRadicalSize_eq_fast (q : Fin (2 ^ 15)) :
    packedRadicalSize q = fastRadicalSize q.val := by
  rw [packedRadicalSize_eq]
  rw [← solutionCount_polarRows_eq_radical]
  exact polar_solutionCount_rank q.val

end Spin.UniformMixer
