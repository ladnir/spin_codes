# A late-group obstruction for fixed RS8 Turbo SPIN

## Claim and scope

Fix the current RS8 local outer, 64 packet regions, independent setup per
group, and any causal inner with zero initial state and zero output before
the first nonzero routed input. Let L=K/128 and N=256L. Under the ideal
setup law, as L grows,

\[
 \Pr_\theta[d_{\min}(E_\theta)\le (5/64)N]\longrightarrow 1.
\]

In particular this fixed construction cannot have a greater-than-10% or
greater-than-11% asymptotic relative-distance theorem with probability
1-o(1). This says nothing adverse about its proved K=2^16, 2^18, and
2^20 finite certificates. It also does not rule out a Turbo-like family
whose local packet-support floor grows or is enforced by construction.

## One fixed message per group

In one RS[8,4] row, use the degree-three polynomial
`p(X)=X(X+1)(X+2)` over GF(16), and set the other seven rows to zero.
The row is specified by its values at 0,1,2,3, so this is one fixed
nonzero 128-bit group message before any setup is sampled. Its eight
outer columns are zero at 0,1,2 and nonzero at 3,4,5,6,7. Direct
evaluation in the implementation's `0x13` basis gives the eight symbols
`(0,0,0,6,1,6,4,5)`, confirming exactly five active 32-bit columns.

For any fixed nonzero column, the independent uniform nonzero GF(2^32)
scalar sends that column to a uniform nonzero 32-bit vector; the binary
adjoint has the same fixed-input law. Exactly `8*15=120` such vectors
have one active four-bit packet. The five column events are independent,
so they all yield one active packet with probability

\[
 q_5=\left(\frac{120}{2^{32}-1}\right)^5>0.
\]

Conditional on this event, the group's uniform 64-packet shuffle places
those five active packets in a uniform five-subset of the 64 regions.
The probability that they occupy exactly the final five regions is

\[
 p_5=\frac{q_5}{\binom{64}{5}}
     \approx 2.2330316172\cdot10^{-45}
     \approx 2^{-148.32776}.
\]

For the less restrictive 10% witness, it suffices that all five land in
the final six regions. Its probability is
`p_6=6p_5` (about `2^-145.74280`).

## Why the inner cannot repair this event

Put the fixed group message above in one group and zero in every other
group. Regional permutations change only the positions *within* each
region. On the `p_5` event, all routed input before the final five regions
is zero. The inner begins with zero state and applies a linear update, so
its state and output remain zero throughout that prefix, regardless of
its expansion, feedback, or sampled state maps. Each region contains
`4L` binary coordinates. Hence the complete output is confined to the
last `5*4L` coordinates and has weight at most

\[
 20L=(5/64)N=0.078125N.
\]

The per-group events use independent column maps and group shuffles.
They do not depend on the shared regional permutations or inner maps.
Thus for L candidate groups,

\[
 \Pr_\theta[d_{\min}(E_\theta)\le(5/64)N]
 \ge 1-(1-p_5)^L\longrightarrow1.
\]

The crossover is astronomical: `-log2(p_5)=148.32776`, and the
rough half-probability point is `log2(K)=154.799`. The finite K20
certificate is therefore not in tension with this asymptotic result.
The obstruction is a true low-weight codeword event, rather than a loose
first-moment bound. It is also independent of the choice of s15 or s22
inner and survives arbitrarily strong causal state randomization.

## Design consequence

With 64 fixed regions, a local packet-support floor of at least seven
removes this particular causal counterexample at a 10% target; at least
eight are needed for an 11% target because `7/64<0.11`. Such floors are
necessary to defeat this mechanism, not sufficient for a full distance
proof. An alternative is a growing outer whose probability of these
late, low-support words vanishes fast enough as L increases. Merely
strengthening the current inner cannot fix an all-zero prefix.

The next proof/design step is to compare cheap local ways to enforce the
eight-region floor, using the existing RS[16,8] outer with minimum nine
active packets as a control. Retain RS8 for the finite performance ladder;
there is no reason to discard its existing certificates.
