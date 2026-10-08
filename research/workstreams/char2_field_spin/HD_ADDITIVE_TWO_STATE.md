# Two additive states with a varying observation

The proposed regional construction is rejected for the 40% target.
An interval cut that allows unsafe holes gives setup-failure probability
at least 0.9985378145964314; see [the finite proof](HD_TWO_STATE_CUTS.md).
The earlier cut using only wholly safe runs missed this obstruction.
The exact local support law and quotient-rank helpers below remain useful
research components. No native implementation or timing was produced.

## Map and implementation opportunity

Work over F=GF(2^16). Independently at each coordinate, sample
(c,d1,d2,f) uniformly from the nonzero tuples satisfying

    c+d1+f*d2 != 0.

Starting with s1=s2=0, compute

    y = s1+f*s2+c*x,
    s1' = s1+d1*x,
    s2' = s2+d2*x.

Retain all y values and discard the final states. All setup remains fixed
for every message. Each 128-bit record contains eight parallel F lanes.
The overall input/output map is triangular with nonzero diagonal c. The
additional displayed condition makes the three-by-three input/state map
invertible as well.

Both state streams are XOR prefixes of independent products. They contain
no serial field multiplication. Compute d1*x, d2*x, and c'*x in batches,
then form the two prefixes and the output, where

    c' = c+d1+f*d2,
    y = s1'+f*s2'+c'*x.

The change c to c' is an involution of the admissible tuple set. It therefore
preserves the setup distribution. This schedule requires four field products
per input record: three input products and one product by the second prefix.
That is an operation count, not a latency prediction.

The observation coefficient f must vary. With fixed f, the single state
r=s1+f*s2 satisfies r'=r+(d1+f*d2)x and y=r+c*x. That special case collapses
to a scalar inner. A checked GF4 fixture has an inverse off-diagonal minor
of rank two when f varies; the new state is not merely redundant storage.

The candidate full encoder uses GRS[64,32] rows, coordinate-major placement,
two independently weighted prefixes, one regional routing stage, and this
inner. Global carries cross region boundaries. A region can finish both
prefixes and the post-inner before the next region starts, since routing
does not cross regions. This scheduling observation does not establish
the full encoder's distance or measured performance.

## Exact support law

Let q=|F|, C=q^2-3q+3, and T=(q-1)^2*C. There are exactly T admissible
coefficient tuples. A support state records which of s1 and s2 are nonzero.
For each fixed (d1,d2), there are C choices of (c,f); hence d1 and d2 are
independent uniform nonzero scalars under the conditioned law.

For zero input, the state is unchanged. State00 outputs zero. States10
and01 output a nonzero value. State11 outputs zero with probability1/(q-1).
For nonzero input, the following counts are divided by T. All values within
the same incoming support class have the same law by coordinate scaling.

| Old state | New state | Total tuple count | Tuples with zero output |
|---|---|---:|---:|
|00|11|T|0|
|10|01|(q-1)C|(q-1)^2|
|10|11|(q-2)(q-1)C|(q-2)^2(q-1)|
|01|10|(q-1)C|(q-1)^2|
|01|11|(q-2)(q-1)C|(q-2)^2(q-1)|
|11|00|C|0|
|11|10 or01, each|(q-2)C|(q-2)^2|
|11|11|(q-2)^2 C|(q-2)(q^2-5q+7)|

`theory/hd_additive_two_state.py` implements this law. Independent enumeration
checks it over GF4, GF8, and GF16. The retained tests exhaust all GF4 local
values, all two-step GF4 setups and messages, and an inverse-tail fixture.
This is a local-law result, not an independence claim about complete words.

## Why the first run cut was insufficient

The inverse has at most two tail dimensions across any interval boundary.
If safe check columns lie in a space U of dimension h, a union of r safe
intervals gives a final check rank of at most h+2r. The scalar proof charged
only h+r. The difference prevents that proof from rejecting40% here.

For orientation, at safe density alpha, selecting all runs of length at
least three has asymptotic expected useful dimension alpha^3*N after
charging two tails per run. Beating h approximately alpha*N/2 requires
alpha>1/sqrt(2). At that threshold, h+2r is approximately56.07% of N,
not below40%. A bounded mean-based longest-run screen also did not reject
40%. These necessary checks are not a distance lower bound; more general
cuts were therefore still necessary. The later interval-with-holes cut
supplies one and rejects the complete proposed regional construction.

## Rank envelope and initial failed screens

A first moment over field-valued messages can heavily overcount a rare
bad subspace. The separate helper `theory/hd_additive_rank_screen.py` instead
counts exact supports and their shortened source dimensions. It uses a
rank mark tau>=1. Every newly imposed zero equation drops the surviving
source dimension by zero or one, unless the branch dies.

The conceptual history retains every previously declared live linear form.
After each zero equation, discard a history if any required live form is
identically zero on the remaining source space. For a fresh nonzero scalar a,
an equation M+a*L with L live has a rank-zero increment for at most one a.
Its expected rank factor is therefore at most

    g = 1/(q-1) + (1-1/(q-1))/tau.

Retaining the required-live scalar L sharpens this factor to

    g_sharp = max(1/(q-1), 1/tau).

If M and L are independent on the current source space, every coefficient
imposes one equation and L remains live. If M is a multiple of L,
every noncancelling coefficient kills L and that history is discarded.
At most one coefficient leaves rank unchanged. These cases give the two
terms in the maximum. This argument requires a specified live scalar;
a state vector known only to be nonzero does not supply that condition.
The optional `--sharp` screen preserves the previous weaker default.

At each position expose the output equation before the two new-state zero
equations. A live input lets its fresh c supply the output pivot. With
zero input and both states live, fresh f supplies that pivot. The two state
zeros use fresh d1 and d2. Only four current state flags index the envelope;
they do not replace the full live-form history in the proof.

For this bound first sample all four coefficients independently from F*.
Keep only support paths allowed by local invertibility, and divide the final
probability bound by Z^N, where

    Z = (q^2-3q+3)/(q-1)^2.

This accounts for conditioning the complete coefficient experiment. AtK20
the cost is approximately46.167 bits. It is paid once for the complete
encoder, not independently at each zero equation.

The first prefix's MDS row transfer uses tau^(w-b+k) for each feasible
nonzero support of size w, and weight1 for the zero row. For two prefixes,
a product-vector envelope majorizes the row transfer by

    T0^(tensor b) + tau^(-(b-k))*(T0+tau*T1)^(tensor b).

Here T0 and T1 are the four-state support operators of two successive
weighted prefixes. A product-vector row ratio is a sum of two products.
It is convex in the entering-state histogram, so its maximum is bounded by
the largest homogeneous histogram. This avoids an exponentially large
state matrix, at the cost of additional proof slack.

All numerical screens currently use a single **uniform global** permutation
between precode and post-inner. Each middle weight h therefore contributes
one factor1/binom(N,h). The faster regional candidate instead needs the
joint regional weights and their separate binomial normalizations. These
global results cannot be used as a certificate for regional routing.

The bounded global screens fail: with one prefix, GRS[128,64], tau16, and
middle density .4, the upper-bound margin is about-958677 bits. With two
prefixes, GRS[64,32], tau64, and density .03, it is about-168587 bits.
Other tested marks also fail. These values are evaluations of loose upper
bounds, not codeword witnesses or impossibility results. No complete class
cover, interval verification, or native implementation was attempted.

The sharper charge does not rescue those global-route bounds. For the
two-prefix [64,32] control, tau256 gives margins about -205737 bits at
middle density .03 and -1611001 bits at density .4. Both fits converged.
The exact regional contraction and its separate per-region normalization
are recorded in [HD_TWO_STATE_REGIONAL.md](HD_TWO_STATE_REGIONAL.md).

## Disposition

Do not implement this candidate for the 40% target or continue tightening
its lower-distance bound. The finite cut is an actual obstruction,
independent of the coefficient values. A replacement must change the
shared row geometry or give the post-inner sufficiently stronger coupling.
The local law, conditioning factor, and live-history proof contract remain
available for a changed construction.

```text
python -B research/workstreams/char2_field_spin/theory/hd_additive_two_state_tests.py
python -B research/workstreams/char2_field_spin/theory/hd_additive_rank_tests.py
python -B research/workstreams/char2_field_spin/theory/hd_additive_rank_screen.py --pre 2 --outer 64 --tau 64 --density .03 .4
```
