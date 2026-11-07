import SpinCodes.UniformMixer.ScalableLiveMaxSpectrum

/-!
# Spectrum interface for the actual three-coordinate transfer

These statements express the four input quantities F, U, H, and u of the
Scalable SPIN transfer in terms of the checked production spectrum. They
do not assert the numerical Collatz or distance gates.
-/

#check Spin.UniformMixer.scalableMaxLiveMoment_eq_spectrum
#check Spin.UniformMixer.scalableMeanLiveMoment_exact_spectrum
#check Spin.UniformMixer.scalableZeroHold_exact_spectrum
#check Spin.UniformMixer.scalableFeedbackActive_exact_spectrum
#print axioms Spin.UniformMixer.scalableMaxLiveMoment_eq_spectrum
#print axioms Spin.UniformMixer.scalableMeanLiveMoment_exact_spectrum
#print axioms Spin.UniformMixer.scalableZeroHold_exact_spectrum
#print axioms Spin.UniformMixer.scalableFeedbackActive_exact_spectrum
