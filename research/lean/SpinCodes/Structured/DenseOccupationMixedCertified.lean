import SpinCodes.Structured.DenseOccupationAllCertified
import SpinCodes.Structured.DenseOccupationScalarCertified
import SpinCodes.Structured.DenseOccupationFourierCertified

namespace Spin.Structured.DenseGeometry

/-- The complete mixed numerical certificate for the actual dense routed experiment. -/
theorem certified_denseRates : DenseOccupationFixed.DenseRates (4/10000000) 24000000000000 := by
  apply denseRates_of_families
  · intro i hi
    exact (occupation_certified i hi).mono_constant (by norm_num)
  · intro i hi
    exact (scalar_certified i hi).mono_constant (by norm_num)
  · exact fourier_certified

#print axioms certified_denseRates
end Spin.Structured.DenseGeometry
