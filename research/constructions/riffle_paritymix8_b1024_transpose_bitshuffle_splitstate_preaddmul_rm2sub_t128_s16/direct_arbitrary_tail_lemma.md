# Direct envelope for arbitrary tail supports

This note closes the mixed class in the endpoint/body decomposition. It works
inside the recorded RM2Sub transfer model. It does not establish the model's
constituent bounds or replace the conjectured outer spectrum.

## Setup

There are `B=1024` transposed regions. Each region contains `R=2048`
positions, one for each Riffle outer block. The construction samples an
independent uniform permutation of the `R` positions in every region.

Fix a set of `b>=1` central-body blocks. Replace their outer words temporarily
by independent ambient-uniform words in `{0,1}^B`. Also fix every other outer
word and every bit of those words. These fixed words include all endpoint-tail
words. The proof averages over the ambient body words, the region
permutations, and the randomness represented by the inner transfer model.

For a Chernoff parameter `z` and `0<=w<=R`, let `F_z[w]` be the nonnegative
`2x2` transfer matrix for one region. The two states record whether the inner
state is zero or live. Each entry upper-bounds the output-weight moment and the
next-state event when the region input has a uniform support of weight `w`.
The analyzer obtains `F_z[w]` by composing the fixed-weight epoch matrices with
the exact hypergeometric weights between epochs.

## Region lemma

Fix one region before its position permutation. Suppose exactly `j` fixed
tail positions contain one. Let `X` be the number of ones among the `b` body
positions. The ambient reference gives

```text
X ~ Bin(b,1/2).
```

The region permutation is independent of `X`. Conditioned on `X=x`, it maps
the `j+x` nonzero positions to a uniform support of size `j+x`. Therefore the
region transfer is entrywise at most

```text
G_z[b,j] = sum_x Pr[X=x] F_z[j+x].
```

The matrices satisfy

```text
G_z[0,j]   = F_z[j],
G_z[b+1,j] = (G_z[b,j] + G_z[b,j+1]) / 2.
```

The recurrence follows by exposing one additional fair body bit. Define the
entrywise envelope

```text
H_z[b] = max_{0 <= j <= R-b} G_z[b,j].
```

The maximum is entrywise; its four entries may use different values of `j`.
Consequently, `H_z[b]` upper-bounds one region for every fixed tail pattern.
It also permits the tail count to vary adversarially between regions.

The construction samples the region permutations independently. Thus the
ambient-reference moment for all `B` regions is at most

```text
[1,0] H_z[b]^B [1,1]^T.
```

Multiplying this moment by `z^(-d)` gives the usual Chernoff upper bound for
output weight at most `d`. Optimizing over the recorded Chernoff grid gives the
quantity called `arbitrary_tail_inner_log2_upper` in the receipt.

## Restoring the body spectrum

Let `U` be the ambient-uniform distribution on one outer word. To define
`Q_body`, sample a uniform nonzero local message and apply the independently
sampled local encoder. Retain the output only when its weight is central-body.
Thus `Q_body` is a subdistribution. After the output-coordinate permutation,
define `L=dQ_body/dU`. For `p>1`, Holder's inequality gives

```text
E_U[1_bad product_i L_i]
  <= E_U[L^p]^(b/p) Pr_U[bad]^((p-1)/p).
```

The complete central-body spectrum computes `E_U[L^p]`. The direct region
lemma supplies `Pr_U[bad]` uniformly for every fixed collection of tail words.
No distributional assumption about their supports is needed at this step.

Let `T` be the expected number of endpoint-tail codewords in one sampled local
DBA code. The construction samples the two DBA interleavers independently for
each Riffle outer block. After choosing the `b` body locations, the expected
number of nonempty tail assignments to the remaining `R-b` locations is

```text
(1+T)^(R-b) - 1.
```

Multiplying this count by the Holder bound and summing over `1<=b<R` covers
every message with at least one body block and at least one tail block.

## Numerical result and scope

At relative distance `0.09`, the aggregate mixed-class margin is
`235.523500071` bits. The dominant term has one body block. Its optimized
inner margin before the outer count is `766.265056786` bits.

The recurrence audit checks two identities at the dominant Chernoff value.
First, `G_z[b,0]` agrees with the existing fair-body transfer table for every
`b`; the maximum log-domain difference is `3.83e-12`. Second, selected
`G_z[b,j]` entries agree with direct binomial sums to `1.78e-15`.

This result is a nearest-binary64 certificate under the recorded transfer
model. A final theorem still requires an outward-rounded calculation, a formal
statement of the inner constituent bounds, and a fixed replacement for the
conjectured EBCH64 spectrum.
