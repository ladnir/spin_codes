# Active-packet collapse

## Exact local reduction

Fix one transposed region and a shared-group profile

\[
n=(n_1,n_2,n_3,n_4).
\]

Let the Bernoulli envelope parameter be `p`, and define

\[
q_r=1-(1-p)^r.
\]

A width-`r` group produces a nonzero packet with probability `q_r`.
PacketMul makes that packet uniform over the 15 nonzero elements of `GF(16)`
without changing whether it is zero. Therefore the inner transfer conditioned
on `K=k` surviving packets is the universal matrix `U_k`. It contains no
reference to `n`.

This removes the packet widths from the inner interface. The widths remain
only in the survivor law

\[
K_n=\sum_{r=1}^4 \operatorname{Binomial}(n_r,q_r).
\]

The profile can be summed out exactly for one region. Mark the number of
active outer blocks by `x` and the number of surviving packets by `y`. One
four-block group contributes

\[
F(x,y)=\sum_{r=0}^4\binom4r x^r
\bigl((1-q_r)+q_r y\bigr).
\]

Using the definition of `q_r`, this simplifies to

\[
F(x,y)=\bigl(1+(1-p)x\bigr)^4
+y\left((1+x)^4-\bigl(1+(1-p)x\bigr)^4\right).
\]

Therefore the coefficient

\[
[x^h y^k]F(x,y)^{2048}
\]

is the complete multiplicity-weighted law of `k` surviving packets at
occupation `h`. No `n_r` variable appears.

## Why the profile does not disappear automatically

The four-block groups are fixed across all 256 regions. Hence one outer block
support selects one profile `n`, and the same profile determines the
survivor law in every region. The required outer sum is

\[
\sum_n N(n)\,\mathbf e_Z^{\mathsf T}R_n^{256}\mathbf 1,
\qquad
R_n=\mathbb E[U_{K_n}].
\]

Averaging `R_n` over profiles before taking the 256th power would resample
the four-block grouping independently in every region. That is a different
ensemble. Thus the one-region active-packet collapse is exact, but the
256-region calculation must preserve the correlation caused by the fixed
outer support.

## Profile-free exponential envelope

For `0<theta<=1`, define the entrywise envelope

\[
M(\theta)_{ab}=\max_k \frac{(U_k)_{ab}}{\theta^k}.
\]

Then

\[
R_n\preceq M(\theta)
\prod_{r=1}^4
\left(1-q_r+q_r\theta\right)^{n_r}.
\]

The scalar factor commutes with the matrix product. For occupation

\[
h=\sum_{r=1}^4 r n_r,
\]

the complete profile sum becomes the coefficient

\[
[x^h]\left(
1+\sum_{r=1}^4
\binom4r
\left(1-q_r+q_r\theta\right)^{256}x^r
\right)^{2048}.
\]

This is a rigorous reduction from four profile variables to one occupation
variable. A stronger version assigns a separate `theta_ab` to each of
the four state transitions and sums binary state paths by their transition
counts. Every resulting outer sum is still one univariate coefficient.

## Numerical audit at occupation 128

The exact profile sum has 7,105.481 bits of margin at the recorded 9% proof
point. The best common-rate active-count envelope has margin

\[
-14{,}897.282\text{ bits}.
\]

The four-rate state-path envelope improves this to

\[
-12{,}034.886\text{ bits},
\]

but remains unusable. Its largest term stays in the zero state for 255
regions and activates in the last region. The exponential envelope loses the
sharp tail information needed to suppress this late-start path.

This failure does not refute the construction. It shows that a single
exponential summary of the active-packet count is too coarse. The exact
profile computation at the same point remains positive.

## Next proof target

Keep `U_k` as the inner interface, but preserve more of the distribution of
`K` near the zero-state boundary. The next useful test is a small piecewise or
coefficient-exact treatment of zero-state runs, followed by the coarser
active-count envelope after the first activation. This targets the observed
loss without restoring the complete `n1,n2,n3,n4` enumeration.

The reproducer is
`scripts/analyze_riffle_ldpcsplitstate_packetmul_active_count_envelope.py`.
Its receipt is `receipts/active_count_envelope_occupation128.json`.
