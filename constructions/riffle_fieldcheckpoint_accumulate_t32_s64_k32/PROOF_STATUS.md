# Distance proof status

The checkpoint randomization gives an exact two-state inner model. It does
not by itself complete the outer-spectrum sum.

## Exact epoch transfer

Fix one epoch input and write \(m=16\) for the number of visits to every
lane. For lane \(i\), let \(b_{i,j}\) be the input-prefix parity at its
\(j\)-th visit and let

\[
 a_i:=\lvert\{j:b_{i,j}=1\}\rvert.
\]

If the state at the start of the epoch is uniform over the nonzero states,
then its output-weight moment is

\[
 \frac{\prod_{i=1}^{64}
       \left(z^{a_i}+z^{m-a_i}\right)
       -z^{\sum_i a_i}}
      {2^{64}-1}.
\]

Let \(v\) be the vector of final input parities in the 64 lanes. Starting
from a live state, the epoch ends in zero only when its initial state equals
\(v\). The transition is impossible when \(v=0\), and otherwise has
probability exactly \(1/(2^{64}-1)\) before the output-weight tilt is
applied.

At the following checkpoint, zero remains zero. For every fixed nonzero
state \(q\), a uniform nonzero multiplier makes \(L_a^Tq\) uniform over the
nonzero states. A fresh checkpoint therefore erases all information about a
nonzero endpoint except that it is nonzero. This establishes exact closure
on the two classes zero and live, including their output-weight moments.

For an empty epoch, put \(D=2^{64}-1\). Its transfer is diagonal, with live
entry

\[
 G_0(z)=\frac{(1+z^{16})^{64}-1}{D}.
\]

For an epoch containing one impulse at a uniform position, the zero-to-live
entry is

\[
 \frac1{16}\sum_{j=1}^{16}z^j,
\]

and the live-to-zero entry is

\[
 \frac1{16D}\sum_{j=0}^{15}z^j.
\]

The remaining entry follows by summing over the other nonzero states. The
implemented formulas agree with exhaustive enumeration at state width three
to within \(6.7\times10^{-16}\).

## Exact transfer for arbitrary epoch occupancy

Condition on exactly \(r\) impulses at a uniform \(r\)-subset of an epoch.
For a subset \(S\) of the 16 visits to one lane, let \(a(S)\) denote the
number of nonzero prefix parities. Define

\[
 A_k(z):=\sum_{\substack{S\subseteq[16]\\|S|=k}}z^{a(S)},
 \qquad
 B_k(z):=\sum_{\substack{S\subseteq[16]\\|S|=k}}z^{16-a(S)}.
\]

Write \(A(u)=\sum_k A_k(z)u^k\), \(B(u)=\sum_k B_k(z)u^k\), and let
\(E_A\) retain the even coefficients of \(A\), while \(O_B\) retains the odd
coefficients of \(B\). If \(M_r(z)\) is the zero/live epoch transfer and
\(C_r={1024\choose r}\), then

\[
 (M_r)_{00}=\frac{[u^r]E_A(u)^{64}}{C_r},
\]

\[
 (M_r)_{01}=\frac{[u^r]\left(A(u)^{64}-E_A(u)^{64}\right)}{C_r},
\]

\[
 (M_r)_{10}=
 \frac{[u^r]\left((E_A(u)+O_B(u))^{64}-E_A(u)^{64}\right)}{D C_r},
\]

and

\[
 (M_r)_{10}+(M_r)_{11}=
 \frac{[u^r]\left((A(u)+B(u))^{64}-A(u)^{64}\right)}{D C_r}.
\]

Thus four powers of degree-16 lane polynomials produce all 1025 occupancy
transfers. Exhaustive checks at state width three and epoch length six agree
to within floating-point roundoff. At the target dimensions and \(z=1\),
all entries are finite and nonnegative, and every row sum differs from one
by at most \(1.2\times10^{-15}\).

## One-active result

An inactive region has transfer \(R_0(z)=M_0(z)^8\). A region with one
uniform impulse has transfer

\[
 R_1(z)=\frac18\sum_{j=0}^{7}M_0(z)^jM_1(z)M_0(z)^{7-j}.
\]

Conditioned on outer weight \(w\), the complete inner moment is

\[
 F_w(z)=
 e_0^T\frac{[u^w](R_0(z)+uR_1(z))^{256}}{{256\choose w}}\mathbf1.
\]

The Chernoff calculation uses

\[
 \Pr[W\le d\mid w]\le \inf_{0<z<1}z^{-d}F_w(z).
\]

At output length \(2^{21}\), the modeled spectrum gives 81.2667 bits of
one-active margin at relative distance 0.09. Weight 38 is dominant. The
40-bit frontier lies between 0.14049 and 0.14050: the respective margins are
40.0007 and 39.9944 bits.

The 9% one-active margin varies little across the tested checkpoint
geometries:

| State bits | Checkpoint steps | Epoch bits | Margin bits |
|---:|---:|---:|---:|
| 32 | 8 | 256 | 81.4799 |
| 64 | 16 | 512 | 81.4407 |
| 64 | 32 | 1024 | 81.2667 |
| 64 | 64 | 2048 | 80.9236 |
| 128 | 128 | 4096 | 80.6411 |

This stability indicates that one-active failure is governed mainly by the
outer word's region placement and late start. Checkpoint extinction is not
the leading term at these parameters.

The weight-two cancellation audit reaches the same conclusion for two
minimum-weight outer words. A weight-two region fails to activate with
probability \(2^{-9.09}\). Nevertheless, nonactivation through region 210
has exponent 240.23 bits. Charging the modeled multiplicity of all
weight-38 two-block messages leaves 166.68 bits before the output-tail bound.
See `CANCELLATION_ANALYSIS.md` for the probability space and calculation.

## Remaining proof problem

If a region contains \(h\) active bits, its exact transfer is

\[
 R_h(z)=\frac{[u^h]
   \left(\sum_{r=0}^{1024}{1024\choose r}M_r(z)u^r\right)^8}
   {{8192\choose h}}.
\]

The remaining difficulty is the outer average. Every active outer word of
weight \(w\) selects a uniform \(w\)-subset of the 256 regions. Several
active outer blocks therefore produce an exchangeable, correlated vector of
region occupancies. The matrices \(R_h(z)\) need not commute, so an unordered
occupancy histogram alone does not immediately determine their ordered
product.

This is now the only conceptual compression problem. The inner contribution
for every fixed occupancy is exact and has dimension two.

A promising next step is to find one positive two-component potential that
simultaneously dominates the live portions of all relevant \(R_h(z)\). The
bound should condition on the first activating epoch, so it retains the late
start cost, and then scalarize the remaining ordered product. It should be
compared against exact one- and two-active calculations before being used in
the full spectrum sum.

## Claim scope

The formulas above are exact for the stated random permutations and
checkpoint multipliers. The reported distance values are ordinary
floating-point diagnostics for one active outer block and a modeled outer
spectrum. No complete distance certificate is claimed.

## Progress toward a full certificate

The complete two-active modeled-spectrum diagnostic now sums all 4278
unordered outer-weight pairs. It includes every support intersection. Its
margin is 159.8280 bits, and the pair \((38,38)\) dominates.

A pointwise spectrum-density envelope gives a stronger compression. Every
regular support in the modeled spectrum has density at most twice the
aggregate density of a uniform random 256-bit outer word. The factor-two
loss applies independently to each regular active block. This reduces all
regular weight tuples to one active-block count.

The envelope calculation for occupations one through 64 has 55.9008 bits of
margin at 9% relative distance. Occupation one dominates. Occupation 64 has
pointwise exponent about \(-3393\) bits.

## Middle-density result

The exponent-safe checker computes every regular region matrix in the log
semiring. It uses the exact hypergeometric recurrence across the eight epochs
of a region. On occupations one through 128, it agrees with the earlier
ordinary-arithmetic calculation; at occupation 64 the two reported exponents
agree to the displayed precision. A small reference instance agrees to
within \(10^{-15}\).

This computation identifies a middle-density obstruction at 9% relative
distance. For 1024 regular active outer blocks, the optimized pointwise
first-moment exponent is approximately \(+16039\) bits. The optimizer is
strictly inside the tested Chernoff grid. Replacing the density envelope by
the exact modeled spectrum can recover at most approximately one bit per
active block, which is much smaller than this deficit. This result does not
prove that the construction has a low-weight word. It proves that the current
first-moment route cannot establish 9% distance.

At 6% relative distance, the same route covers every regular occupation.
The ranges are as follows.

- Occupations 1 through 128 inherit at least the 55.9008-bit bound already
  proved at the larger 9% threshold.
- Occupations 129 through 1024 have 1947.16 bits of aggregate margin.
- Occupations 1025 through 5466 have 834.4399 bits of aggregate margin. The
  dominant count is 1237.
- Starting at occupation 5467, an invertibility and Hamming-ball-volume bound
  has 107.1058 bits of aggregate margin.

The unique all-one outer word requires separate treatment. If a region has
\(a\) fair regular candidates and one forced impulse, its transfer satisfies

\[
 F_{a,1}(z)=2F_{a+1,0}(z)-F_{a,0}(z).
\]

The log-domain finite difference agrees with the independent mixed enumerator
to within \(2\times10^{-15}\). For one all-one block and 65 through 5487
regular blocks, the aggregate margin at 6% is 1161.5591 bits. The earlier
mixed calculation covers every total occupation through 64 at the larger 9%
threshold. Starting at 5488 regular blocks, the dense-subspace bound sums
every possible number of all-one blocks and has 91.3725 bits of aggregate
margin.

The conditional 6% certificate therefore remains open only for configurations
with at least two all-one outer words and fewer than 5488 regular outer words.
The retained calculations must then be repeated with outward-rounded
arithmetic. See `proof/FULL_CERTIFICATE_PLAN.md`.
