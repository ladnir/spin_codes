import SpinCodes.UniformMixer.QuadraticRankHistogram

/-! Kernel-checked rank histogram, quadratic coefficients 20480 through 24575. -/

namespace Spin.UniformMixer

set_option maxRecDepth 100000 in
set_option maxHeartbeats 0 in
theorem rank_histogram_block_40 :
    fastRankHistogramBlock 20480 512 = (0, 16, 304, 192) := by
  decide +kernel

set_option maxRecDepth 100000 in
set_option maxHeartbeats 0 in
theorem rank_histogram_block_41 :
    fastRankHistogramBlock 20992 512 = (0, 0, 256, 256) := by
  decide +kernel

set_option maxRecDepth 100000 in
set_option maxHeartbeats 0 in
theorem rank_histogram_block_42 :
    fastRankHistogramBlock 21504 512 = (0, 16, 304, 192) := by
  decide +kernel

set_option maxRecDepth 100000 in
set_option maxHeartbeats 0 in
theorem rank_histogram_block_43 :
    fastRankHistogramBlock 22016 512 = (0, 0, 256, 256) := by
  decide +kernel

set_option maxRecDepth 100000 in
set_option maxHeartbeats 0 in
theorem rank_histogram_block_44 :
    fastRankHistogramBlock 22528 512 = (0, 0, 256, 256) := by
  decide +kernel

set_option maxRecDepth 100000 in
set_option maxHeartbeats 0 in
theorem rank_histogram_block_45 :
    fastRankHistogramBlock 23040 512 = (0, 16, 304, 192) := by
  decide +kernel

set_option maxRecDepth 100000 in
set_option maxHeartbeats 0 in
theorem rank_histogram_block_46 :
    fastRankHistogramBlock 23552 512 = (0, 0, 256, 256) := by
  decide +kernel

set_option maxRecDepth 100000 in
set_option maxHeartbeats 0 in
theorem rank_histogram_block_47 :
    fastRankHistogramBlock 24064 512 = (0, 16, 304, 192) := by
  decide +kernel

end Spin.UniformMixer
