import SpinCodes.UniformMixer.QuadraticRankHistogram

/-! Kernel-checked rank histogram, quadratic coefficients 4096 through 8191. -/

namespace Spin.UniformMixer

set_option maxRecDepth 100000 in
set_option maxHeartbeats 0 in
theorem rank_histogram_block_8 :
    fastRankHistogramBlock 4096 512 = (0, 16, 304, 192) := by
  decide +kernel

set_option maxRecDepth 100000 in
set_option maxHeartbeats 0 in
theorem rank_histogram_block_9 :
    fastRankHistogramBlock 4608 512 = (0, 16, 304, 192) := by
  decide +kernel

set_option maxRecDepth 100000 in
set_option maxHeartbeats 0 in
theorem rank_histogram_block_10 :
    fastRankHistogramBlock 5120 512 = (0, 16, 304, 192) := by
  decide +kernel

set_option maxRecDepth 100000 in
set_option maxHeartbeats 0 in
theorem rank_histogram_block_11 :
    fastRankHistogramBlock 5632 512 = (0, 16, 304, 192) := by
  decide +kernel

set_option maxRecDepth 100000 in
set_option maxHeartbeats 0 in
theorem rank_histogram_block_12 :
    fastRankHistogramBlock 6144 512 = (0, 0, 256, 256) := by
  decide +kernel

set_option maxRecDepth 100000 in
set_option maxHeartbeats 0 in
theorem rank_histogram_block_13 :
    fastRankHistogramBlock 6656 512 = (0, 0, 256, 256) := by
  decide +kernel

set_option maxRecDepth 100000 in
set_option maxHeartbeats 0 in
theorem rank_histogram_block_14 :
    fastRankHistogramBlock 7168 512 = (0, 0, 256, 256) := by
  decide +kernel

set_option maxRecDepth 100000 in
set_option maxHeartbeats 0 in
theorem rank_histogram_block_15 :
    fastRankHistogramBlock 7680 512 = (0, 0, 256, 256) := by
  decide +kernel

end Spin.UniformMixer
