# High distance after freezing the high-rate winner

The [83.33%-rate checkpoint](GF8_HIGH_RATE_WINNER.md) is complete for now.
The new objective is substantially higher symbol distance, ideally close to
50%, starting at rate 1/2. Use requested K20, 128-bit application records,
one regional routing stage, and a 40-bit whole-code setup-failure target.
None of the new results below is a complete distance certificate or a
native encoder timing.

## Alphabet choice changes the target

At rate 1/2, the random-linear-code benchmark solves `H_q(delta)=1/2`,
where `H_q` is the q-ary entropy function. Numerical inversion gives:

| Scalar alphabet | Asymptotic random-code distance benchmark |
|---|---:|
| GF(2^8) | 38.0472% |
| GF(2^16) | 43.8191% |
| GF(2^32) | 46.8838% |
| GF(2^64) | 48.4386% |
| GF(2^128) | 49.2189% |

These values are comparison points, not distance upper bounds or finite
guarantees for SPIN. See [Guruswami's GV notes](https://www.cs.cmu.edu/~venkatg/teaching/codingtheory/notes/notes2.pdf).
The function `gv_distance` in the new screen reproduces them.

Sixteen byte lanes sharing a GF256 code do not inherit the last row of this
table. For a fixed code defined over GF256, scalar extension preserves its
minimum symbol distance. Approaching the larger-alphabet benchmark therefore
requires a genuinely different sampled code, not a wider payload type.

The existing [medium-alphabet comparison](MEDIUM_ALPHABETS.md) already provides
correct GF32 arithmetic controls. Uniform GF32 multipliers were substantially
more expensive than byte multipliers in that isolated pass. That result
discourages a blind substitution, but does not measure a redesigned encoder.

## A GF256 starting candidate

Keep the frozen profile's four-input/four-state full-Cauchy map and its
independent nonzero input and state diagonals. Change the outer to
GRS[64,16] and retain only output coordinates 0 and 1 from each packet.
The state remains continuous across the 64 independently permuted regions.
This gives exactly rate 1/2 at K20:

| Quantity | Value |
|---|---:|
| Input dimension | 1,048,576 |
| Outer rows | 65,536 |
| Raw inner inputs | 4,194,304 |
| Retained outputs | 2,097,152 |
| Initial distance target | 36% |
| Required minimum output weight | 754,975 |

Simply halving the rate with GRS[64,32] and retaining all four outputs did
not pass the sparse screen at 30%. GRS[64,24] with three retained outputs
also failed that tested bound. The stronger GRS[64,16] outer with two retained
outputs improves the source spreading while keeping final rate 1/2.
It pays for that improvement with twice as many raw inputs as final outputs.

For total active-row count `A`, the floating sparse screen at 36% gives:

| A | Sampled class margin |
|---:|---:|
| 64 | 149.32 bits |
| 96 | 150.96 bits |
| 128 | 175.42 bits |

The sampled middle-weight classes from raw-input density .03 through .95
also have positive bounds. At 37%, the sampled `A=64` margin falls to
11.8 bits. These numbers concern individual classes, not their union.
Untested classes, coverage joins, and directed rounding remain unresolved.
A failed upper bound does not exhibit a low-distance codeword.

The screen changes the exact forbidden output-weight threshold in the
existing full-Cauchy transfer. It does not edit the frozen 10% proof helpers.
Three regression tests check agreement at 10%, threshold monotonicity, and
the random-code reference. The 36%, `A=64` result was independently rerun.

```text
python research/workstreams/char2_field_spin/theory/high_distance_restart_screen_tests.py
python research/workstreams/char2_field_spin/theory/high_distance_restart_screen.py --outer 64 16 --outputs 2 --distance 9/25 --active 64 96 128 --density .03 .10 .30 .50 .80 .95 --iterations 45
```

## Next decision

An inexpensive intermediate control is GF(2^16) in the quadratic tower
`u^2=u+delta` over AES GF256, with `delta=0x20`. For coefficient `a+b*u`
and input `x0+x1*u`, multiplication gives

```text
low  = a*x0 + delta*b*x1
high = b*x0 + (a+b)*x1
```

Setup can precompute byte pairs `[a,a+b]` and `[delta*b,b]`. Two GFNI
products, a byte swap, and XOR then implement the product across packed
16-bit lanes. The scalar formula matched the existing quadratic-tower
oracle for all 65,536 coefficients on all sixteen binary basis inputs.
This validates the algebra, not a native kernel or its performance. Uniform
sampling must exclude the one zero coefficient. The fixed GF256 outer and
Cauchy maps still act bytewise in this representation.

Use the 36% GF256 candidate as a proof control, not the final near-50%
answer. Before a full closure campaign, compare the same topology using
genuinely larger-field diagonals. Include GF(2^16) as the low-cost intermediate
control and GF(2^32) as the first near-50% comparison, with
GF(2^64) and GF(2^128) retained for the closer-to-50% target. Keep the fixed outer and
Cauchy maps over an embedded GF256 where that makes arithmetic cheaper.
Compare both optimized bounds and matched rank tilts. If the failing
classes choose tilts below 255 already, more coefficient values alone may
not address their loss.

In parallel, examine structured coefficient sets that mix the larger
alphabet. A small set of distinct field elements gives only its actual
pivot-cancellation probability; it must not borrow the probability of
uniform full-field sampling. The rank-potential and dense-membership
arguments must be rederived wherever they use the full coefficient law.
Do not invest in a full encoder until sparse and middle classes both look
credible and the complete arithmetic and memory costs are accounted for.

Keep regional coordinate spreading. The earlier
[uniform-route suffix obstruction](HIGH_DISTANCE_GLOBAL.md) is a genuine
construction failure and does not disappear with a larger coefficient
alphabet or a stronger causal inner.
