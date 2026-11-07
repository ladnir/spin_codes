import SpinCodes.UniformMixer.QuadraticRankHistogram

/-! Kernel-checked rank histogram, quadratic coefficients 28672 through 32767. -/

namespace Spin.UniformMixer

set_option maxRecDepth 100000 in
set_option maxHeartbeats 0 in
theorem rank_histogram_block_56 :
    fastRankHistogramBlock 28672 512 = (0, 16, 304, 192) := by
  decide +kernel

set_option maxRecDepth 100000 in
set_option maxHeartbeats 0 in
theorem rank_histogram_block_57 :
    fastRankHistogramBlock 29184 512 = (0, 0, 256, 256) := by
  decide +kernel

set_option maxRecDepth 100000 in
set_option maxHeartbeats 0 in
theorem rank_histogram_block_58 :
    fastRankHistogramBlock 29696 512 = (0, 0, 256, 256) := by
  decide +kernel

set_option maxRecDepth 100000 in
set_option maxHeartbeats 0 in
theorem rank_histogram_block_59 :
    fastRankHistogramBlock 30208 512 = (0, 16, 304, 192) := by
  decide +kernel

set_option maxRecDepth 100000 in
set_option maxHeartbeats 0 in
theorem rank_histogram_block_60 :
    fastRankHistogramBlock 30720 512 = (0, 0, 256, 256) := by
  decide +kernel

set_option maxRecDepth 100000 in
set_option maxHeartbeats 0 in
theorem rank_histogram_block_61 :
    fastRankHistogramBlock 31232 512 = (0, 16, 304, 192) := by
  decide +kernel

set_option maxRecDepth 100000 in
set_option maxHeartbeats 0 in
theorem rank_histogram_block_62 :
    fastRankHistogramBlock 31744 512 = (0, 16, 304, 192) := by
  decide +kernel

set_option maxRecDepth 100000 in
set_option maxHeartbeats 0 in
theorem rank_histogram_block_63 :
    fastRankHistogramBlock 32256 512 = (0, 0, 256, 256) := by
  decide +kernel

end Spin.UniformMixer
