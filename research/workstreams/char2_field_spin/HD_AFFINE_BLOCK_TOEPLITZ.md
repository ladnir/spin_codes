# Affine block Toeplitz outer proposal

Four independent Toeplitz4 blocks can replace the bit-syndrome mixer and
both outer diagonals at comparable arithmetic cost. The blocks use a
different fixed partition of the sixteen row wires in successive columns.
Their pointwise atom bound yields a complete-group counting interface
without a union over constituent-row sets. This is a proposed ensemble,
not a distance certificate or a measured implementation.

## Map and fixed wiring

Work over F=GF256. Each of the 8192 groups still contains sixteen
RS[16,8] rows, hence 128 input and 256 output symbols. Keep the current
256-region paired route, fair pair swaps, and single persistent inner
state. There are no independent inner chains.

Identify row labels with pairs (u,v) in GF4 squared. This GF4 is used
only to specify fixed wire partitions. Its four elements may be represented
by 0,1,omega,omega+1, with omega squared=omega+1. Define five partitions:

\[
 \mathcal P_\infty=\{\{(b,v):v\in\mathrm{GF4}\}:b\in\mathrm{GF4}\},
 \qquad
 \mathcal P_s=\{\{(u,su+b):u\in\mathrm{GF4}\}:b\in\mathrm{GF4}\}.
\]

Here s ranges over GF4. Use the ordered list
(infinity,0,1,omega,omega+1) three times, followed by infinity, for the
sixteen RS columns. Each partition contains four disjoint four-wire blocks.
Every pair of distinct wires belongs to a common block in exactly one
of the five partitions.

Independently for every source group, column, and block, sample a
Toeplitz4 matrix uniformly conditioned on invertibility. Apply it to the
four RS symbols in that block. Write its four outputs back to those
same ordered wires. This replaces D2 H* D1 entirely; it does not add
Toeplitz blocks around H*. Setup is sampled once and reused by the code.

The partitions are deterministic. They do not resample the active row
set in each column. A native implementation can name the four inputs
and outputs explicitly, so the partition requires no record shuffle.
Whether the compiler preserves that schedule remains an implementation
question.

## Pointwise whole-group bound

Put q=256 and

\[
 \rho=\frac{q^{-4}}{1-q^{-1}}.
\]

The established [Toeplitz atom argument](HD_GROUPED_OUTER.md#invertible-toeplitz-column-mixers)
gives, for every fixed nonzero x in F^4 and every y in F^4,
Pr[T x=y]<=rho. The four block samplers are independent. Invertibility
also preserves whether an entire block is zero.

Fix an arbitrary target word y in F^256 before exposing its mixers.
For column j, let b_j in {0,1,2,3,4} count its nonzero output blocks
under the actual partition for that column. Let S8(y) be the sum of
the eight largest values among b_0,...,b_15. Break ties by column index.
Then the expected number of nonzero source messages encoding to y is

\[
 \mu(y)\le
 \mathbf 1\{\#\{j:b_j>0\}\ge9\}\,\rho^{S_8(y)}.
\]

To prove the claim, condition on all mixers in the eight columns with
smallest b_j. Inverting those columns determines eight symbols of each
RS row, hence at most one complete source message. For that candidate,
the unexposed block inputs are fixed. Each nonzero target block requires
an independent hit of probability at most rho. A zero target block
requires a zero block input; inconsistency only reduces the probability.
Multiplying the remaining hit bounds gives rho raised to S8(y).
Finally, a nonzero source has a nonzero RS row, which occupies at least
nine columns. Every column mixer is invertible, so y occupies at least
nine columns too.

This conditions on columns chosen from the fixed target, not from the
random setup. No message-count factor, row-set factor, or additional RS
density is missing: the eight conditioned columns already determine the
entire 128-symbol source. Independence across source groups gives the
usual product comparison. Existing receipts cannot be relabelled as
receipts for this new sampler.

## Cheap scalar proof kernel

Let 0<=t<=1 be a common scalar occupancy atom from an inner comparison.
The marked mass of a nonzero four-symbol block is

\[
 u(t)=(1+255t)^4-1.
\]

For a column-activity histogram n=(n_0,...,n_4), put
L(n)=sum_b b n_b and h(n)=16-n_0. Let S8(n) be the sum of the largest
eight entries of the corresponding sixteen-entry activity sequence.
The complete-group scalar majorant is the positive polynomial

\[
 G(t)=\sum_{\substack{n_0+\cdots+n_4=16\\n_0\le7}}
 \frac{16!}{\prod_b n_b!}
 \prod_{b=0}^4\binom4b^{n_b}
 \rho^{S_8(n)}u(t)^{L(n)}.
\]

There are only binom(20,4)=4845 histograms before the nine-column
restriction. Precompute their contributions by total L to obtain a
degree-64 positive polynomial in u. This counts the labelled column and
block choices once. The scalar mark is symmetric, so each fixed
partition has the same block enumerator; no symmetry of the actual code
or resampling of rows is asserted.

If the unchanged regional operators satisfy
R[m,f]v<=eta t^m v for all required occupancies and source profiles,
with v_0=1, the fixed-r source union is bounded by

\[
 N_r\,\frac{z^{-734003}}{255\min_i v_i}\,
 \eta^{256}G(t)^r.
\]

Here N_r is the number of source-group sets in the declared scope.
For all sets it is binom(8192,r). The pointwise outer comparison is
constant across amplitudes of a fixed support, so its comparison measure
has independent uniform nonzero amplitudes conditional on that support.
This permits the existing forced-nonzero regional interface. The formula
contains no extra factor 255^(-8a).

The simpler inequality S8>=L/2 gives

\[
 G(t)\le\sum_{h=9}^{16}\binom{16}{h}
 \left[(1+\sqrt\rho\,u(t))^4-1\right]^h.
\]

That relaxation can be costly exactly where the column floor matters.
For nine full columns, S8=32 but L/2=18: it loses rho^(-14), about 448
bits per group. Use the exact 4845-histogram polynomial in the first
control. It is too small to justify paying that relaxation.

## Cost and decision boundary

The existing Toeplitz4 circuit uses nine variable GFNI products and
fifteen binary XORs. Four blocks therefore cost 36 variable products and
60 XORs per sixteen-symbol column. The normalized bit-syndrome path
uses 31 variable diagonal products, eight fixed products, and 66 XORs.
The proposed arithmetic saves three product sites and six XORs per column.
Both retain the sixteen RS[16,8] rows and unchanged payload traffic.

Prepared keys increase from 31 to 36 bytes per source group and column:
3.875 MiB becomes 4.5 MiB at K20. More variable-key traffic and a different
dependency schedule may erase the arithmetic saving. These counts are
not a latency prediction.

The inspected history contains larger dense/Toeplitz mixers and the
[B128 smaller-group screen](HD_SHIFTED_B128.md), but not this fixed
256-region schedule. B128 changed the number and length of regions;
this proposal does not. Its earlier medium-weight failure remains a
warning about the entropy cost of shorter constituent codes.

The first useful test is the exact scalar G with a cached, explicitly
scoped inner bank, compared with the half-score relaxation. A failure
would first reject that scalar comparison, not the sampled code. The
fixed partition schedule may support stronger row-floor constraints,
but none are charged in G. The zero-GF deterministic input-order proposal
and its separate rare-outer-event bound remain cheaper priorities before
any native implementation of this fallback.
