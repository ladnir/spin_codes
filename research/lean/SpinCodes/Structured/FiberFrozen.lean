import SpinCodes.Structured.FiberFrozenData
import SpinCodes.Structured.FiberNumericsBridge
import SpinCodes.Structured.SparseProgramDefs

namespace Spin.Structured.FiberNumerics

theorem kernel_frozen_eq (j : Fin 129) :
    (SparsePolynomial.weightN j).kernel = Data.kernels.getD j 0 := by
  have hj : j.val < SparsePolynomial.Data.weights.length := by
    rw [SparsePolynomial.Data.weights_length]
    exact j.isLt
  rw [SparsePolynomial.weightN,
    List.getD_eq_getElem SparsePolynomial.Data.weights _ hj,
    ← Data.kernels_frozen_checked,
    List.getD_eq_getElem (SparsePolynomial.Data.weights.map (·.kernel)) _ (by simpa using hj)]
  simp only [List.getElem_map]

theorem cap_frozen_eq (j : Fin 129) :
    (SparsePolynomial.weightN j).cap = Data.caps.getD j 0 := by
  have hj : j.val < SparsePolynomial.Data.weights.length := by
    rw [SparsePolynomial.Data.weights_length]
    exact j.isLt
  rw [SparsePolynomial.weightN,
    List.getD_eq_getElem SparsePolynomial.Data.weights _ hj,
    ← Data.caps_frozen_checked,
    List.getD_eq_getElem (SparsePolynomial.Data.weights.map (·.cap)) _ (by simpa using hj)]
  simp only [List.getElem_map]

end Spin.Structured.FiberNumerics
