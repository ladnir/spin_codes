# Riffle RM512-RandomStepConv g=4 sigma=20

Status: PAUSED_MODEL_EXPLORATION

This model groups the 16384 data symbols into 4096 groups of four. Each group
forms the 256-bit input to one exact RM(4,9) \([512,256,32]\) constituent.
The two 64-bit field-parity symbols retain one 256-bit rate-half output block.
The outer therefore still produces 524352 four-bit packets.

The RM weight spectrum is exact and authenticated against its Gleason
formula. Packet support after the local permutation is counted exactly
conditional on binary weight. The determined parity block is bounded
worst-case.

At 9% relative output weight, the floating first-moment upper bound is
\(+3003.95\). The comparable Outer256 model gave \(+4710.17\). The exact
RM512 spectrum therefore improves the diagnostic by 1706.22 bits, but it does
not close the bound.

The bulk \(h=32767\) obstruction disappears under a pointwise inner tilt: its
relaxed contribution becomes \(-7787.91\). The remaining sampled obstruction
is near data support 95 and uses about six live episodes.

An explicit multi-episode family has expected-count exponent \(+8.28\) when
the parity block has zero packet support. Four parity packets make this family
negative. The current question is therefore the joint spectrum of one RM data
group and its two field parities.

See REPORT.md, receipts/rm512_randomstepconv_distance09.json, and
receipts/inner_episode_structure_audit.json.

This model was paused after the RM minimum distance was judged too weak for
the structural proof gym. The successor is
`Riffle FrozenRandom512-P2-RandomStepConv g=4 sigma=20`.
