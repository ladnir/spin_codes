# Finite (k=2^{20}) EBCH32--ParityFanout--BA setup

## Status and purpose

This document freezes the probability space and the outer-spectrum identity
for the power-of-two finite candidate.  It does not claim the final distance
certificate.  That claim additionally requires outward certificates for every
occupation (Q=1,\ldots,8192), the RM2Sub interface audit, and an efficient
test for the setup event defined below.

## Parameters

Set

\[
 k=2^{20},\qquad B=256,\qquad L=8192,\qquad
 N=BL=2^{21},\qquad d=\lfloor 0.11N\rfloor=230686.
\]

One outer row contains eight independent coordinate blocks of length 32.
Each block uses a fixed binary linear \([32,16,8]\) extended BCH code with
weight enumerator

\[
 A(z)=1+620z^8+13888z^{12}+36518z^{16}
      +13888z^{20}+620z^{24}+z^{32}.
\]

Thus the unrandomized direct-sum row code has dimension 128 and exact weight
enumerator (A(z)^8).  No shortening is needed at this geometry.

## Row setup experiment

The following random choices are mutually independent between rows and
between stages unless stated otherwise.

1. Sample (S\) uniformly from the 31-subsets of \([256]\).
2. Conditional on (S\), sample (T\) uniformly from the 33-subsets of
   \([256]\setminus S\).
3. Apply the linear parity-fanout map
   \[
   F_{S,T}(x)=x+\left(\sum_{i\in S}x_i\right)1_T.
   \]
4. Sample a uniform permutation of the 256 coordinates and apply a terminated
   rate-one accumulator.
5. Independently repeat the preceding permutation-and-accumulator stage.

Because (S\cap T=\varnothing), the parity on (S\) is unchanged by the
addition on (T\).  Hence (F_{S,T}^2\) is the identity.  Every sampled
fanout map is therefore an invertible linear map.  Each subsequent
permutation and terminated accumulator is also invertible.  Every row code
has dimension 128 for every setup outcome.

The phrase "independent row draw" below means an independent draw of all four
row-random objects: ((S,T)) and the two accumulator interleavers.  A
per-constituent permutation before fanout is omitted.  It has no effect on the
law because the constituent is fixed and (S,T) are sampled uniformly over
the full row.

## Exact expected spectrum identity

Let (C_w=[z^w]A(z)^8\).  Fix an input word of weight (w).  Write

\[
 a=|S\cap\operatorname{supp}(x)|,
 \qquad b=|T\cap\operatorname{supp}(x)|.
\]

The source-overlap probability is

\[
 p_w(a)=
 \frac{\binom wa\binom{256-w}{31-a}}{\binom{256}{31}}.
\]

If (a) is even, fanout leaves the word unchanged.  If (a) is odd, then,
conditional on (a),

\[
 p_w(b\mid a)=
 \frac{\binom{w-a}{b}
       \binom{256-31-(w-a)}{33-b}}
      {\binom{256-31}{33}},
\]

and the output weight is (w+33-2b\).  With the convention that an invalid
binomial coefficient is zero, the exact expected multiplicity immediately
after fanout is

\[
 D_h=\sum_{w=0}^{256} C_w\left(
  1_{h=w}\sum_{a\text{ even}}p_w(a)
  +\sum_{a\text{ odd}}p_w(a)
       \sum_b p_w(b\mid a)1_{h=w+33-2b}
 \right).
\]

This identity includes the unique all-ones word.  For (w=256), one has
(a=31) and (b=33) deterministically, so its fanout image has weight
(256+33-66=223).  The word does not remain at weight 256 and is not removed
by assumption.

For a terminated accumulator preceded by a uniform interleaver, the exact
weight transition is

\[
 P_{u,h}=\frac{
   \binom{h-1}{\lceil u/2\rceil-1}
   \binom{256-h}{\lfloor u/2\rfloor}}
  {\binom{256}{u}}
\]

for (u>0), with (P_{0,0}=1); invalid binomial coefficients again denote
zero.  Consequently, the expected final row spectrum is exactly

\[
 \bar A_h=\sum_{u=0}^{256}\sum_{v=0}^{256}
 D_u P_{u,v}P_{v,h}.
\]

All identities above are finite sums of rational numbers.  The companion
verifier evaluates nonnegative upper bounds with directed binary64 rounding
and checks the exact integer Vandermonde normalizations used by the fanout
law.

## Conditioning law

Let (Z\) be the number of nonzero final row-code words whose weights lie
outside

\[
 W=[24,232].
\]

The accepted-row event is (G_{256}=\{Z=0\}\).  Markov's inequality gives

\[
 \Pr[G_{256}]\ge
 1-\sum_{h\in\{1,\ldots,23,233,\ldots,256\}}\bar A_h.
\]

The construction-level setup law is independent rejection sampling: draw a
row as above until (G_{256}) holds, independently for each of the 8192
rows.  The mathematical conditional law is therefore exact once a decision
procedure for (G_{256}) is available.  The expected number of trials per
row is at most the reciprocal of the displayed lower bound.

An efficient exact or one-sided-safe decision procedure for (G_{256}) has
not yet been supplied.  The probability calculation alone does not establish
a linear-time setup algorithm.

## Whole-code probability space

After the 8192 accepted row draws, sample the existing region permutations,
bit-transpose/bit-shuffle permutations, split-state choices, and nonzero
RM2Sub-S19 multipliers according to their specified independent laws.  The
finite distance probability is over this entire setup experiment.  Message
words are enumerated; they are not sampled.

The target statement is

\[
 \Pr_{\text{setup}}[d_{\min}\le 230686]\le 2^{-40}.
\]

The first-moment proof partitions nonzero parent messages by the number
(Q\) of active outer rows.  A complete certificate must cover every
(Q=1,\ldots,8192) and every admitted row-weight mixture, without sampled
gaps.  Only then may Markov's inequality convert the expected number of bad
words into the target probability.

## Length rule

This frozen candidate is an admissible power-of-two point, not yet an
arbitrary-length family.  Its native message length is (LB/2=2^{20}), and
its native output length is (LB=2^{21}).  General lengths require the
separate adjacent-size, shortening, padding, or puncturing wrapper in
`ARBITRARY_LENGTH_WRAPPER.md`; no such wrapper is implicit in this finite
claim.

## Proof ledger

- Proved algebraically here: the fanout map is invertible; the exact fanout
  overlap identity; the all-ones transition (256\mapsto223); the two-stage
  accumulator spectrum identity; and preservation of row dimension.
- Machine-checkable once the receipt is present: an outward upper bound on
  the rejected expected spectrum and a lower bound on
  \(\Pr[G_{256}]\).
- Diagnostic only at present: the claimed 11% margins obtained from
  nearest-binary64 optimization.
- Open: outward certificates for all occupations and mixtures; the
  RM2Sub/parity-fanout implementation-interface audit; an efficient
  (G_{256}) test; and an implementation benchmark for ordinary and
  transposed encoding.
