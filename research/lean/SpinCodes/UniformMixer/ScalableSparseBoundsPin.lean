import SpinCodes.UniformMixer.ScalableSparseFeedbackUpper
import SpinCodes.UniformMixer.ScalableSparseLiveUpper

/-! Audited analytic bounds feeding the three sparse Collatz inequalities.
These are bounds on H, G, u, and F; the U bound and final arithmetic gate
remain open. -/

#check Spin.UniformMixer.scalableZeroHold_sparse_quadratic
#check Spin.UniformMixer.spectrumG_sparse_upper
#check Spin.UniformMixer.spectrumu_sparse_upper
#check Spin.UniformMixer.spectrumF_sparse_upper
#print axioms Spin.UniformMixer.scalableZeroHold_sparse_quadratic
#print axioms Spin.UniformMixer.spectrumG_sparse_upper
#print axioms Spin.UniformMixer.spectrumu_sparse_upper
#print axioms Spin.UniformMixer.spectrumF_sparse_upper
