import SpinCodes.UniformMixer.QuadraticRankHistogram

/-! Kernel-checked rank histogram, quadratic coefficients 12288 through 16383. -/

namespace Spin.UniformMixer

set_option maxRecDepth 100000 in
set_option maxHeartbeats 0 in
theorem rank_histogram_block_24 :
    fastRankHistogramBlock 12288 512 = (0, 16, 304, 192) := by
  decide +kernel

set_option maxRecDepth 100000 in
set_option maxHeartbeats 0 in
theorem rank_histogram_block_25 :
    fastRankHistogramBlock 12800 512 = (0, 16, 304, 192) := by
  decide +kernel

set_option maxRecDepth 100000 in
set_option maxHeartbeats 0 in
theorem rank_histogram_block_26 :
    fastRankHistogramBlock 13312 512 = (0, 0, 256, 256) := by
  decide +kernel

set_option maxRecDepth 100000 in
set_option maxHeartbeats 0 in
theorem rank_histogram_block_27 :
    fastRankHistogramBlock 13824 512 = (0, 0, 256, 256) := by
  decide +kernel

set_option maxRecDepth 100000 in
set_option maxHeartbeats 0 in
theorem rank_histogram_block_28 :
    fastRankHistogramBlock 14336 512 = (0, 0, 256, 256) := by
  decide +kernel

set_option maxRecDepth 100000 in
set_option maxHeartbeats 0 in
theorem rank_histogram_block_29 :
    fastRankHistogramBlock 14848 512 = (0, 0, 256, 256) := by
  decide +kernel

set_option maxRecDepth 100000 in
set_option maxHeartbeats 0 in
theorem rank_histogram_block_30 :
    fastRankHistogramBlock 15360 512 = (0, 16, 304, 192) := by
  decide +kernel

set_option maxRecDepth 100000 in
set_option maxHeartbeats 0 in
theorem rank_histogram_block_31 :
    fastRankHistogramBlock 15872 512 = (0, 16, 304, 192) := by
  decide +kernel

end Spin.UniformMixer
