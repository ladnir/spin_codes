import SpinCodes.UniformMixer.ScalableTransferPaperOrientation

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
#check Spin.UniformMixer.scalableThreeTransfer_eq_spectrum
#check Spin.UniformMixer.productionSpectrum_weighted_sum
#check Spin.UniformMixer.spectrumU_polynomial
#check Spin.UniformMixer.spectrumH_polynomial
#check Spin.UniformMixer.spectrumu_polynomial
#check Spin.UniformMixer.transfer_total_eq_paper_first
#check Spin.UniformMixer.scalable_bernoulli_paper_collatz_bound
#print axioms Spin.UniformMixer.scalableMaxLiveMoment_eq_spectrum
#print axioms Spin.UniformMixer.scalableMeanLiveMoment_exact_spectrum
#print axioms Spin.UniformMixer.scalableZeroHold_exact_spectrum
#print axioms Spin.UniformMixer.scalableFeedbackActive_exact_spectrum
#print axioms Spin.UniformMixer.scalableThreeTransfer_eq_spectrum
#print axioms Spin.UniformMixer.productionSpectrum_weighted_sum
#print axioms Spin.UniformMixer.spectrumU_polynomial
#print axioms Spin.UniformMixer.spectrumH_polynomial
#print axioms Spin.UniformMixer.spectrumu_polynomial
#print axioms Spin.UniformMixer.transfer_total_eq_paper_first
#print axioms Spin.UniformMixer.scalable_bernoulli_paper_collatz_bound
