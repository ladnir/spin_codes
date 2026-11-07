import SpinCodes.UniformMixer.ScalableMomentSpectrumExact

/-!
# Exact mean-moment spectrum audit

This checks the structural equality feeding the Scalable SPIN numerical
transfer. It does not assert a numerical distance gate.
-/

#check Spin.UniformMixer.productionSpectrum_sum
#check Spin.UniformMixer.scalableMeanLiveMoment_exact_spectrum
#print axioms Spin.UniformMixer.productionSpectrum_sum
#print axioms Spin.UniformMixer.scalableMeanLiveMoment_exact_spectrum
