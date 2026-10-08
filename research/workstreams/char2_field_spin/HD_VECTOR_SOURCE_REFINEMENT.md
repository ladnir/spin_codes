# The vector bound is not limited by row minimum weight

At the saved K20 witnesses, enforcing every active outer row's minimum
weight cannot materially improve the failed vector-state rank bound.
A directed calculation bounds the possible gain by less than `2^-42` bits
for 16 active rows and `2^-16` bits for 64 active rows. The respective
deficits are about 582 and 5,256 bits.

This is a diagnosis of the positive proof expression in
[HD_VECTOR_OBSERVATION_PROOF.md](HD_VECTOR_OBSERVATION_PROOF.md), not a
lower bound on the encoder's actual failure probability. It does not rule
out every reoptimized source bound. No new construction or native timing
is introduced here.

## A distinguished row gives an exact diagnostic

Use the earlier half-rate GRS[128,64] geometry over GF65536, with state
dimension eight. Fix an active-row count `A`, an output mark `z`, and a
rank mark `tau`. There are `R=16384` rows and `b=128` regions. Let `K_h`
be the positive rank kernel for a region with `h` live inputs, including
its regional binomial normalization.

The relaxed source bound includes every binary support pattern on the
`A` selected rows. A pattern of total weight `W` receives source weight
`tau^W`, apart from the common factor `tau^(-64A)`. Its other weight is
the nonnegative terminal mass of the chronological regional matrix product.
Normalize these weights to obtain an auxiliary probability measure on
support patterns. This is an analytical measure, not the sampled code's
setup distribution.

Distinguish one of the `A` rows. In each region, its bit is either zero or
one. Summing the other `A-1` row bits gives two matrices:

\[
 R_0=\sum_{h=0}^{A-1}\binom{A-1}{h}\tau^h K_h,\qquad
 R_1=\sum_{h=0}^{A-1}\binom{A-1}{h}\tau^{h+1}K_{h+1}.
\]

The scalar polynomial

\[
 J(v)=(1,0)(R_0+vR_1)^b\binom11
\]

records the distinguished row's degree. Its normalized coefficient of
`v^w` is exactly the probability that this row has weight `w` under the
auxiliary measure. The matrices preserve chronological state carry.
At `v=1`, Pascal's identity gives the original relaxed regional matrix,
so `J(1)` is its complete terminal mass.

This costs a one-variable matrix polynomial of degree `b`, rather than
enumeration over all `2^b` row supports or all regional occupancy profiles.
The helper uses polynomial powering and truncates only above the required
degree.

## Bound the benefit of enforcing every row

An active GRS[128,64] row has weight at least 65. Define

\[
 p_{\rm bad}=\frac{\sum_{w=0}^{64}[v^w]J(v)}{J(1)}.
\]

The auxiliary measure is invariant under permuting the `A` row labels.
Each row therefore has the same failure probability `p_bad`. Let `G` be
the event that all row degrees are at least 65. A union bound gives

\[
 \max(0,1-Ap_{\rm bad})\le\Pr[G]\le1-p_{\rm bad}.
\]

No independence or negative association between rows is assumed. In
particular, the inner's matrix weights can correlate their degrees.

Filtering the positive source sum to `G` implements the exact per-row
minimum-weight restriction. If the original bound is `B`, the filtered
expression is exactly `B Pr[G]`. Thus, when `A p_bad<1`, its improvement
in margin is at most

\[
 -\log_2(1-Ap_{\rm bad}).
\]

This bounds the value of exact row truncation at fixed marks, including
any method that computes that same restriction more efficiently.

## Results at the saved witnesses

The row-degree summaries below use the previous floating witnesses without
optimization. The minimum feasible degree is 65.

| Active rows | Mean degree | Degree variance | Probability a specified row is too light | Maximum all-row-filter gain |
|---:|---:|---:|---:|---:|
| 16 | 104.1247 | 19.4523 | about `1.2761e-15` | about `2.9457e-14` bits |
| 64 | 93.2529 | 25.3146 | about `3.3730e-8` | about `3.1144e-6` bits |

The helper also replays nearby fixed rational marks with 4,096-bit outward
Arb arithmetic:

| Active rows | Output mark `z` | Rank mark `tau` |
|---:|---|---|
| 16 | `595180216727/596607537200` | `11070105190285/993243241034` |
| 64 | `982888335359/989914484955` | `2850663908977/519167271057` |

The directed replay proves gains below `2^-42` and `2^-16` bits,
respectively. After exact row filtering, the corresponding proof expressions
remain at least `2^582` and `2^5256`. These are lower bounds on an upper-bound
formula, not evidence that the actual code fails with those probabilities.

The total-weight mark previously optimized to one because the relaxed
measure already places almost all its mass above the required row weights.
The distinguished-row calculation shows that this was not hiding a large
population of individually invalid rows compensated by heavier neighbors.

## Exact field-valued multiplicities are also a small correction

For the separate full-field first moment, consider a feasible support of
size `w` in an MDS[b,k] row. Its shortened space has dimension `w-b+k`.
Each supported coordinate is a nonzero linear functional on that space.
A uniform vector therefore has a zero in any specified coordinate with
probability `1/q`. The union bound gives

\[
 q^{w-b+k}(1-w/q)\le a_w\le q^{w-b+k},
\]

where `a_w` is the exact number of words with that support.
After restricting every row to feasible support, replacing the upper
multiplicities by their exact values can improve the first moment by at
most

\[
 -A\log_2(1-b/q).
\]

At `b=128,q=65536`, this is below 0.046 bits for `A=16` and 0.181 bits
for `A=64`. This statement concerns field-valued multiplicities on feasible
supports; it is separate from the rank-envelope row-degree measure above.

In the rank bound, `tau` is a real potential parameter, not a field size.
The initial weight `tau^(w-b+k)` already uses the exact shortened dimension.
Substituting an MDS enumerator with formal field size `tau` would require a
new argument and is not justified by the full-field counting formula.

## A more relevant refinement: retain state-image rank

The two-state envelope records only whether the whole state is zero or
live. It forgets the dimension of the span of its coordinate forms after
earlier zero equations. That information changes both observation and return
costs.

For example, a live input leaving a zero state produces `S'=dX`. On this
branch all state forms are multiples of the required-live scalar `X`, so
the state-image rank is one. During subsequent zero inputs, a zero output
that preserves liveness occurs with probability exactly `1/q`. The current
live-state envelope instead uses its worst case over all positive image
ranks, roughly `1/tau`. It discards this rank-one information at birth.

A second exact calculation quantifies the return cost. On the current
surviving source space, let `X` be a required-live scalar form and define

\[
 h=\dim\bigl(\operatorname{span}(S_1,\ldots,S_r,X)
                    /\operatorname{span}(X)\bigr).
\]

Impose the return equations `S_j+d_jX=0`, with uniform fresh `d` in F^r.
Modulo `X`, these equations have rank `h` for every `d`. Exactly `q^h`
values of `d` give a consistent solution with `X` nonzero. Those values
impose rank `h`; every other value makes `X` vanish on the solution space.
Consequently, retaining only input liveness gives the exact expected rank
factor

\[
 q^{h-r}\tau^{-h}.
\]

Additional earlier live-history requirements can only lower this value.
The old uniform charge `tau^-r` loses a factor `(q/tau)^(r-h)` whenever
`h<r`. The helper checks this formula by enumerating all update vectors
over GF4 for several dependent and independent state/input configurations.

A complete refinement would propagate the state-image rank and its overlap
with each new input form. The formulas above do not supply that transition
system, nor do they close a distance certificate. They identify a concrete
source of slack that is much larger than missing row truncation.

## Reproduction and next step

```text
python -B research/workstreams/char2_field_spin/theory/hd_vector_source_tests.py
python -B research/workstreams/char2_field_spin/theory/hd_vector_source.py
python -B research/workstreams/char2_field_spin/theory/hd_vector_source.py --directed
```

Five tests pass. They compare the distinguished-row polynomial with complete
small support enumeration, compare logarithmic and exact arithmetic, replay
the directed ceilings, check MDS multiplicity inequalities, and enumerate
the conditional return-rank moment.

The earlier vector helper and witnesses are unchanged. No broad parameter
search, native build, or benchmark was run. If this route is revisited,
state-image rank should be the next proof variable. Exact MDS row-degree
profiling alone is not a promising repair at the current witnesses.

## Interface bridge to fixed-Cauchy packets

The cheaper 2K controls in [HD_FAST_PROOF.md](HD_FAST_PROOF.md) deserve a
rank-aware proof before the full-vector inner receives native work.
Their fixed maps use fresh nonzero diagonal scalars, rather than fresh
uniform coefficient vectors. The vector transition probabilities therefore
cannot be copied directly.

Write a fixed-Cauchy packet as

\[
 X=D_xx,\quad R=D_sS,\quad Y=FX+A_eR,\quad S'=R+CX.
\]

Two local interfaces nevertheless carry over.

First, in a quiet packet (`x=0`), the state update is `S'=D_sS`.
If the old state-image rank is `H` and the selected output equations have
rank `rho` on that image, the successor state-image rank is exactly
`H-rho`. A branch with `rho=H` violates old-state liveness and is discarded.
Repeated quiet restrictions cannot keep charging against an
unchanged unrestricted live state. For births from zero, the successor rank
is at most the number of active input forms, not automatically the full
state dimension.

Second, a singleton input admits a sharp return interface. Condition on
its fresh input diagonal and call the resulting required-live form `X`.
The relevant column `c` of `C` has all entries nonzero. A return imposes
`D_s S=-cX`, so every state coordinate must be live on a surviving witness.
Set `h=dim(span(S_1,...,S_s,X)/span(X))`. The affine image of `S` on
`X=1` has dimension `h`. An injective projection onto `h` coordinate
positions shows that it contains at most `(q-1)^h` full-support points.
The vector `-D_s^-1 c` is uniform in `(F*)^s`. Each feasible point imposes
rank `h`, giving the return-moment upper bound

\[
 (q-1)^{h-s}\tau^{-h},\qquad 1\le\tau\le q-1.
\]

Expose the return constraints before output restrictions when using this
bound; the state diagonal must still be fresh. After a singleton return,
`Y=(F-A_eC)X`. The full-Cauchy return map has a nonzero entry in every
retained row of that column, so all retained outputs are live. If `m`
outputs are retained, multiply the return bound by `z^m`.

The recommended first proof control is width/state 4/4 with all four
outputs retained. Prove the quiet and singleton rank transitions first;
leave multi-input packets under the existing safe envelope. Width/state
8/8 is a secondary control, not another simultaneous parameter search.
Before claiming composition, establish a simultaneous weighted upper-tail
bound for successor ranks. Independent worst-case bounds for each rank
entry need not preserve the continuation cone.

This bridge is a local proof assessment, not a new transfer implementation
or certificate. Packet-aligned cut obstructions still need to be checked
for the chosen geometry. No additional screen or timing was run.
