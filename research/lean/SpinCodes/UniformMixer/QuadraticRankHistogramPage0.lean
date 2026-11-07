import SpinCodes.UniformMixer.QuadraticRankHistogram

/-! Kernel-checked rank histogram, quadratic coefficients 0 through 4095. -/

namespace Spin.UniformMixer

set_option maxRecDepth 100000 in
set_option maxHeartbeats 0 in
theorem rank_histogram_block_1 :
    fastRankHistogramBlock 512 512 = (0, 16, 304, 192) := by
  decide +kernel

set_option maxRecDepth 100000 in
set_option maxHeartbeats 0 in
theorem rank_histogram_block_2 :
    fastRankHistogramBlock 1024 512 = (0, 16, 304, 192) := by
  decide +kernel

set_option maxRecDepth 100000 in
set_option maxHeartbeats 0 in
theorem rank_histogram_block_3 :
    fastRankHistogramBlock 1536 512 = (0, 16, 304, 192) := by
  decide +kernel

set_option maxRecDepth 100000 in
set_option maxHeartbeats 0 in
theorem rank_histogram_block_4 :
    fastRankHistogramBlock 2048 512 = (0, 16, 304, 192) := by
  decide +kernel

set_option maxRecDepth 100000 in
set_option maxHeartbeats 0 in
theorem rank_histogram_block_5 :
    fastRankHistogramBlock 2560 512 = (0, 16, 304, 192) := by
  decide +kernel

set_option maxRecDepth 100000 in
set_option maxHeartbeats 0 in
theorem rank_histogram_block_6 :
    fastRankHistogramBlock 3072 512 = (0, 16, 304, 192) := by
  decide +kernel

set_option maxRecDepth 100000 in
set_option maxHeartbeats 0 in
theorem rank_histogram_block_7 :
    fastRankHistogramBlock 3584 512 = (0, 16, 304, 192) := by
  decide +kernel

end Spin.UniformMixer
