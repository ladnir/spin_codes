import SpinCodes.Structured.FiberNumericsData.Tables
import SpinCodes.Structured.SparseWeights

/-! The proposed kernel and cap arrays are exactly the frozen arrays already
used by the numerical transfer matrix. This is an identity of literal data;
its mathematical identification with the actual map is separate. -/

namespace Spin.Structured.FiberNumerics.Data
set_option maxRecDepth 100000

theorem kernels_frozen_checked : SparsePolynomial.Data.weights.map (·.kernel) = kernels := by
  decide +kernel

theorem caps_frozen_checked : SparsePolynomial.Data.weights.map (·.cap) = caps := by
  decide +kernel

end Spin.Structured.FiberNumerics.Data
