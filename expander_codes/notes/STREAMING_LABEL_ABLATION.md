# Column-label ablation for the streaming field encoder

## Question

The streaming field encoder assigns a nonzero label to each expander
coordinate before convolution. The first implementation sampled a full-field
label. The current implementation samples a sign in `{+1,-1}`. Signing uses
conditional field negation instead of field multiplication.

This note asks whether full-field column labels improve the sampled codes. It
compares independent labels, SplitMix-style labels, affine labels, random
signs, and the all-one sequence. Every comparison holds the topology, edge
signs, convolution taps, and public seed fixed.

## Uniform-tap calculation

Fix one message and let `z[t]` denote its expander output at coordinate `t`.
Let `lambda[t]` be a nonzero column label and define

```text
x[t] = lambda[t] z[t].
```

For fresh uniform feedback coefficients `a[t,1],...,a[t,m]`, write the
convolution recurrence as

```text
y[t] = x[t] + sum_j a[t,j] y[t-j].
```

If the preceding state is nonzero, the inner product with the fresh feedback
vector is uniform in the field.  The conditional distribution of `y[t]` then
does not depend on `x[t]`.  If the preceding state is zero, then `y[t]=x[t]`.
Because `lambda[t]` is nonzero, `x[t]` is zero exactly when `z[t]` is zero.

Thus, for a fixed expander output, the distribution of the convolution's zero
pattern does not depend on the nonzero column labels.  This statement assumes
fresh uniform feedback vectors.  It does not apply directly to a deterministic
counter realization.  Independent edge labels remain different: they act
before colliding edges are summed and can change whether `z[t]` is zero.

## Reduced experiments

The meet-in-the-middle search used ternary `[42,21]` codes with degree 6/3,
region size seven, and memory four.  It selected 16 zero coordinates in each
of 32 trials.  Over seeds 1 through 4, the mean weights found were:

| column-label sequence | mean found weight |
|---|---:|
| independent nonzero labels | 9.25 |
| SplitMix-style labels | 9.00 |
| affine labels | 8.25 |
| all-one labels | 9.00 |

Extending the independent and all-one cases through seed 8 gave mean weight
8.75 for both.  Every matrix had full row rank.  This search gives upper bounds
on the sampled distances; it does not certify minimum distance.

Exact enumeration supplied three smaller controls:

| field and dimensions | seeds | independent mean | sign mean | all-one mean | independent/sign/all-one weight-one samples |
|---|---:|---:|---:|---:|---:|
| `F_3`, `[20,10]`, degree 10/5 | 128 | 3.8359 | 3.8359 | 3.8984 | 2 / 4 / 4 |
| `F_3`, `[18,9]`, degree 6/3 | 512 | 3.9453 | 3.8438 | 3.9043 | 3 / 3 / 6 |
| `F_5`, `[12,6]`, degree 6/3 | 512 | 3.2305 | 3.2344 | 3.1895 | 11 / 11 / 9 |

The first and third experiments contained rank-deficient topology samples.
Every label sequence had the same deficient seeds because nonzero column
scaling cannot change rank. Neither the random-sign sequence nor the all-one
sequence has a consistent distance penalty across the three controls.

A sparse-message screen used degree 26/13, region size 29, memory four, and
`F_127`.  For each of 20 seeds, it exhausted every generator row and every
projective combination of two rows.  Independent labels gave mean pair weight
686.95 and minimum 675 at length 754.  All-one labels gave mean 688.10 and
minimum 684.  The screen found no lighter pair caused by all-one labels.

The exact libOTe generator audit supplies a control for the implemented sign
stream. Across eight seeds, its mean pair weight was 667.5 at region size 29.
The previous full-field implementation gave 667.875 on the same seeds. At
region size 79, the corresponding means were 1810.875 and 1812.125. All 32
matrices in these two before-and-after comparisons had full row rank.

## Performance

The C++ benchmark used the Goldilocks 26/13/m4 profile with
`k=1,048,567` and `n=2,097,134`. The current row reports the best of 30 serial
runs; the earlier rows report the previous ten-run measurements.

| mode | scalar encode | degree-two extension encode | receiver pair |
|---|---:|---:|---:|
| full-field labels, earlier implementation | 0.061 s | 0.111 s | 0.153 s |
| all-one ablation | 0.057 s | 0.084 s | 0.115 s |
| random signs, current implementation | 0.057 s | 0.101 s | 0.122 s |

Random signs reduce the receiver-pair time by about 20 percent relative to the
earlier measurement. They retain most of the all-one speedup on that path.
The paired encoder applies each sign after the coordinate completes its
convolution updates. This fusion removes a separate pass over three field
elements. Single-stream encoding uses a separate sign pass because it was
faster in benchmarks.

## Implementation status

The C++ API does not expose the label ablations as separate encoders.
`RegularEcStreamingFieldCode` is the sole heuristic and uses random coordinate
signs. `RegularEcFieldCode` is the reference implementation. It uses balanced
regional assignments and independently generated nonzero edge labels. The
benchmark selects the heuristic by default; `-reference` selects the reference
implementation.

The all-one experiment remains useful evidence about the effect of column
labels. It is not an implementation option or an instance of the proved
edge-labeled ensemble. The ideal calculation assumes fresh uniform taps that
are independent of the labels. A joint attack on the public topology, sign
sequence, and deterministic tap stream remains outside these tests. No
finite-field protocol currently instantiates either encoder.
