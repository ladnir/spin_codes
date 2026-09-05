# Exact RM512 outer with one-lap RandomStepConv

## Result

The exact RM(4,9) \([512,256,32]\) spectrum improves the 9% first-moment
upper bound, but it does not close it:

\[
\log_2\mathbb E[Z_{0.09}]\le 3003.95.
\]

The comparable Outer256 model gave \(4710.17\). Doubling the local input and
output widths therefore recovers 1706.22 bits.

The dominant data-packet support remains in the interval 16384--32767. Its
occupation mode is 343 active RM groups.

These values are floating diagnostics. A positive upper bound does not refute
the construction.

## Construction model

The model groups four 64-bit data symbols into each 256-bit RM input. The 4096
data groups produce

\[
4096\cdot512=2097152
\]

bits. The two 64-bit field-parity symbols retain one 256-bit rate-half output
block. The total remains 2097408 bits, or 524352 four-bit packets.

The data constituent is RM(4,9). Its complete weight spectrum has integer
mass \(2^{256}\), minimum distance 32, and complement symmetry. The committed
spectrum matches the Type-II Gleason formula exactly.

The parity block is not assigned the RM spectrum. Its value is determined by
the data, and its correlation with the data blocks is not enumerated. The
upper bound therefore discards its support when the support tilt is below one
and charges all 64 packets when the tilt exceeds one.

## Full-spectrum accounting

The evaluator sums every RM block-weight and packet-support configuration.
For RM multiplicity \(A_w\) and exact conditional packet-support probability
\(P_{w,s}\), it expands

\[
\left(1+\sum_{w>0}\sum_s A_wP_{w,s}x^s\right)^{4096}.
\]

Thus the calculation does not replace an active block by its minimum weight.

An unconditioned active RM block has mean binary weight 256 and mean packet
support 120. Its central 98% weight interval is 228--284.

Conditioning toward the dominant low-output saddle changes the mean weight to
203.51 and the mean packet support to 95.40. Its median weight is 204, and its
central 98% weight interval is 172--236. The saddle selects moderately lighter
but still numerous RM words rather than minimum-weight words.

The conditioned mean number of active groups is 343.48, with standard
deviation 17.74. Its 1%--99% interval is 303--385 active groups.

## Upper-bound audit

The representative relaxed data support is \(h=32767\). At this support, the
exact RM data coefficient and the one-lap inner bound give:

| Quantity | Exponent |
|---|---:|
| Relaxed RM data coefficient | 77769.24 bits |
| Inner probability upper bound | -74777.34 bits |
| Pointwise first-moment upper bound | 2991.90 bits |
| Whole support-segment upper bound | 3003.95 bits |

The segment aggregation loses 12.05 bits relative to the relaxed pointwise
calculation. The support-tail Chernoff step accounts for 9.74 bits. The
reciprocal-binomial secant has no loss at \(h=32767\), which is the segment
endpoint.

Consequently, the positive result does not come from coarse outer
configuration accounting. It comes from the balance between an exact large
RM coefficient and the current high-support inner probability upper bound.

The episode audit later corrected this interpretation. The interval envelope
reuses parameters optimized at support 16384 throughout the interval. At
support 32767, pointwise optimization changes the inner upper exponent from
\(-74777.34\) to \(-85557.15\). The apparent bulk contribution changes from
\(+2991.90\) to \(-7787.91\). Thus the bulk obstruction was an
interval-parameter artifact.

## Concrete late-placement test

The FFT coefficient counts data support \(32767\). The determined parity
block may add 64 packets. The lower event therefore charges full support

\[
H=32831.
\]

Fix any counted data configuration. Consider the following event over the
global packet permutation and the independent inner maps:

1. All active packets land in the final 100477 packet positions.
2. The state remains nonzero after the first active packet.
3. The emitted binary weight is at most 188766.

This event forces a zero prefix of 423875 packet positions. Its probability is
at least \(2^{-86569.36}\) for every parity value. Combining the event with
the relaxed RM data coefficient gives expected-count exponent

\[
77769.24-86569.36=-8800.12.
\]

The event is therefore not a refutation witness. It requires zero
post-activation terminations.

The old interval inner upper bound at nearby data support has exponent
\(-74777.34\). A pointwise inner tilt removes almost all of this gap.

## Episode structure audit

After pointwise optimization, the largest sampled contribution moves to data
support \(h=95\). The exact relaxed RM coefficient has exponent 228.92. The
fixed-support inner upper bound has exponent \(-197.23\), giving pointwise
sum \(+31.69\).

One active RM group contributes the full 228.92-bit outer exponent to
displayed precision. Two active groups are 40.89 bits smaller, and higher
occupations are negligible. The residual question is therefore a one-group
joint-spectrum problem.

The termination-count decomposition has mode five, mean 5.66, and central 80%
range four through eight. The obstruction is therefore not one uninterrupted
late suffix. It consists of several short live episodes whose total live
length is about 94161 packet positions.

The audit also constructs a probability lower bound. Fix full support \(H\),
termination count \(e\), and total live length \(L\). Count trajectories with
\(e\) terminated episodes followed by one episode that survives to the final
position. Their number is

\[
\binom{L-1}{e}
\binom{N-L+e}{e}
\binom{L-e-1}{H-e-1}.
\]

The first factor chooses positive episode lengths. The second chooses the
off gaps. The third places non-start active packets inside live episodes.
Multiplying by

\[
2^{-20e}(1-2^{-20})^{L-e}
\cdot
\Pr[\operatorname{Bin}(4L,1/2)\le188766]
\]

and dividing by \(\binom NH\) gives the event probability. An exhaustive
state-path enumeration verifies the trajectory count for all instances
through seven packet positions.

For data support 95 and a zero-support parity block, summing termination
counts zero through 15 gives expected-count exponent \(+8.28\). The dominant
terms have five or six terminations. This is explicit refutation evidence for
the parity-relaxed model.

The same family is sensitive to parity support:

| Parity packet support | Expected-count exponent |
|---:|---:|
| 0 | 8.28 |
| 1 | 5.90 |
| 2 | 3.51 |
| 3 | 1.13 |
| 4 | -1.26 |
| 64 | -145.02 |

Thus four parity packets kill this explicit family. The current analysis does
not determine the parity support associated with the RM words contributing at
data support 95.

## Interpretation

Moving from 256-bit modeled blocks to exact 512-bit RM blocks improves the
wide-envelope upper bound by 1706 bits. Pointwise inner optimization then
removes the apparent bulk obstruction entirely.

Late activation remains the live-length mechanism, but the relevant explicit
family uses about six live episodes at data support 95. Whether this family
applies to the construction depends on the determined parity block.

The next proof/refutation step should compute or bound the joint spectrum of
one RM data group and its two field parities. In particular, it should bound
the RM multiplicity at data packet support near 95 when parity support is at
most three. The parity-zero kernel is the first case.
