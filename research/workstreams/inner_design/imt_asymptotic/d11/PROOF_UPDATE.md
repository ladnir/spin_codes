# IMT at 11%: refined outer bound

The refined outer-spectrum certificate closes the remaining dense inequalities
at relative distance 0.11. This extends the
[10.99% proof draft](../ASYMPTOTIC_PROOF.md) without changing the construction.
The resulting 11% argument remains a proof draft pending analytic review;
the numerical certificates do not machine-check the written reductions.

## Claim and unchanged parameters

Use the shared Golay--BA-3 outer, randomized bit-transpose permutation, and
IMT inner defined in the base draft. For L=128m, choose the least multiple
of 24 satisfying b >= (39/4)log2(Lb), and put N=Lb. The IMT maps remain the
fixed balanced expansion and `weight5_seed0` feedback with t=128,s=19.
There is one independent sampled transvection per epoch, persistent state
across regions, and no flush.

The updated claim is

```text
Pr[d_min <= floor(0.11 N)] = o(1),     m -> infinity.
```

Probability is over the one sampled outer, route permutations, and IMT
mixers, all shared across messages. Rate remains exactly 1/2; ordinary and
transposed encoding retain the base draft's O(N) work bound. This does not
assert 11% for the finite BCH-256 instance or provide a finite margin.

## Why the old outer bound was insufficient

At full occupancy, the exact uniform-input inner moment gives the exponent

```text
a_hat(x)-h(x)+ln(2)/2 + h(0.11)-ln(2)/2.
```

The second term is approximately -5.82533613e-5. The coarse piecewise-linear
outer bound sometimes spent more than this amount between its support points.
For example, at x approximately 0.4812883, the resulting exponent was positive.
Changing the inner is unnecessary if a finer certified outer bound removes
that interpolation slack.

## Refined concave majorant

Retain all affine supports of the old majorant. Add five left-half supports
with slopes 0.09, 0.07, 0.05, 0.03, and 0.01. For each slope s, the proposed
intercept is an outward-rounded rational upper bound on

```text
log(1+exp(-s)) - ln(2)/2 + 1/100000.
```

The entropy-tangent expression proposes useful supports; it does not prove
that they bound the BA spectrum. `outer_majorant.py` independently checks
each new active interval against the BA variational objective

```text
sup_{a,b} [g(a)+p(a,b)+p(b,x)].
```

The functions g and p and their outward box formulas are those in the frozen
outer-majorant verifier. Every feasible a,b is included. The initial domain
for each support is [0,1]^2 times its exact rational active interval.
Recursive binary splits prove the objective lies strictly below the line.
The saved leaf paths describe a complete partition; replay rejects missing,
duplicate, or overlapping paths and recomputes every leaf bound.

The left active intervals are approximately

| Slope | Active interval |
|---|---|
| 0.09 | [0.4754601242, 0.4800108262] |
| 0.07 | [0.4800108262, 0.4850046232] |
| 0.05 | [0.4850046232, 0.4900014164] |
| 0.03 | [0.4900014164, 0.4950002083] |
| 0.01 | [0.4950002083, 0.4977500052] |

Exact endpoints and coefficients are in `OUTER_REFINED.json`. Each inherited
active interval is contained in its previously certified interval. The
central constant remains the code-size bound. Reflection is valid because
p(b,x)=p(b,1-x); it does not assume that each realized BA code contains the
all-ones word.

The lower envelope has 39 segments. It is continuous and concave because it
is the minimum of affine functions, and it is no larger than the old bound.
All active pieces are certified, so it still majorizes the BA exponent.
The old uniform finite-enumerator argument consequently applies to this
refined majorant as well. No additional setup event is needed.

The producer processed 10,721 boxes and accepted 5,363 leaves, with no
unresolved domain. All leaf bounds passed replay at higher working precision.
The maximum accepted upper residual is approximately -1.7073533154e-9.

## Dense certificate and composition

`dense_refined.py` partitions the old 11% diagnostic cover at the new outer
breakpoints and refines five boxes further. It proves exact coverage of
alpha in [10^-4,1] and mean row weight x in [0.104,0.896]. As in the base
draft, each fixed-witness exponent is convex in (alpha,alpha*x) within one
affine outer segment. Four endpoint inequalities therefore cover each box.

The resulting 1,023 boxes comprise 433 occupation-transfer witnesses,
307 Fourier witnesses, and 283 scalar uniform-input witnesses. Every box
passes at 256 bits and replays at 512 bits. The maximum outward exponent
per output bit is approximately -4.1033541376e-7. There are no unresolved
boxes. This exponent is not a finite security margin.

The sparse part of the base draft already used delta=0.11 in its inequalities.
In particular, its uniform coefficient is below -0.01340384, including the
old outer likelihood excess and support cost. Its exact interval certificate,
cutoff Q=4096, and geometric union are unchanged. The Q=1,2 and fixed-Q>=3
outward checks also already used delta=0.11. They may continue using the old,
looser outer majorant; no common majorant across regimes is required.

Thus only the positive-occupancy bound in the base draft changes. Substitute
the refined majorant, delta=0.11, and the new negative dense exponent. The
same finite-cutoff union proves the updated claim, subject to the analytic
review obligations of the base draft.

## Reproduction and status

From the repository root, with the existing local evidence and dependencies:

```text
python -B workstreams/inner_design/imt_asymptotic/d11/outer_majorant.py --verify workstreams/inner_design/imt_asymptotic/d11/OUTER_REFINED.json --output OUTER_REPLAY_NEW.json
python -B workstreams/inner_design/imt_asymptotic/d11/dense_refined.py --mode replay --input workstreams/inner_design/imt_asymptotic/d11/DENSE_OUTWARD.json --output DENSE_REPLAY_NEW.json
python -B workstreams/inner_design/imt_asymptotic/d11/verify_d11.py --output ASSEMBLY_D11_NEW.json
python -B -m unittest discover -s workstreams/inner_design/imt_asymptotic/d11 -p "test_*.py" -v
```

Use unused output paths; all producers refuse overwrites. To rediscover the
cover, run `dense_refined.py --mode cover`; then use `--mode certify` on that
output. The verifier authenticates old sparse evidence and recomputes its
11% numerical coefficients, but the complete sparse polynomial replay remains
available through the base assembly command. No benchmark is implicit.

No implementation, paper theorem, or global default changes here. The next
step is analytic review of the combined proof, followed by paper integration.
