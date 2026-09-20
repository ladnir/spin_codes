# Transfer-weighted EA route

## Objective

Reduce the XOR cost of the certified one-stage sparse expand--accumulate
outer without demanding simultaneous upper bounds on every realized weight
shell.

The target instance remains

\[
k=2^{20},\qquad N=2^{21},\qquad
D=228{,}590>0.109N.
\]

The approved outer family is unchanged.  One setup attempt samples a
right-degree-\(r\) sparse map

\[
E:\mathbb F_2^{256}\longrightarrow\mathbb F_2^{512},
\]

applies the zero-initialized accumulator (A), and sets (C=AE).  Bounded
setup accepts the first full-rank constituent.  The encoder repeats that one
accepted constituent in all 4096 outer rows.  This route does not resample a
constituent per row.

## Direct object

Let (A_w(C)) be the nonzero-message spectrum of (C).  For occupation one,
let (t_1(w)) be the probability upper bound supplied by the uniform routing
and RandomStepConv transfer for an outer word of weight (w).  Define

\[
F_1(C):=4096\sum_{w=1}^{512}A_w(C)t_1(w).
\]

Conditional on (C), (F_1(C)) upper-bounds the occupation-one bad-word
probability.  Therefore

\[
\mathbb E_C F_1(C)
=4096\sum_{w=1}^{512}\mathbb E_C[A_w(C)]t_1(w).
\tag{1}
\]

Equation (1) is valid for one sampled constituent.  It uses linearity of
expectation, not independent resampling across outer rows.

If (Z(C)) is the number of nonzero kernel messages, then

\[
\Pr[\operatorname{rank}C<256]\le \mathbb E Z(C).
\]

For an accepted attempt, the right side of (1) may be divided by
(1-\mathbb E Z(C)).  Sixteen independent attempts add an abort probability
at most ((\mathbb E Z(C))^{16}).

## Why this can improve the degree

The existing degree-33 certificate first forces every realized shell into a
simultaneous cap.  It also forces all weights outside 42 through 470 to have
zero multiplicity.  These requirements are stronger than the final SPIN
union bound.  A rare low-weight outer word is harmless when its multiplicity
times its actual inner transfer is small enough.

The direct functional retains that compensation.  The first experiment is
therefore the exact-in-model occupation-one sweep in
`evaluate_ea_combined_q1.py`.

## Reuse obstruction

Occupation one is not the whole distance proof.  When (Q\ge2) outer rows
are active, the same sampled constituent is evaluated on several local
messages.  Averaging only (F_1(C)) is then insufficient.  The next proof
object is a finite family of positive transfer functionals (F_j(C)) whose
products dominate every sparse occupation.  The required setup statement is
a joint moment or high-probability bound for that family.

For (Q=2), the exact ordered-pair EA law is already known: a message-pair
type has three Krawtchouk biases, and the accumulated pair has four states.
This gives an exact route to second moments of transfer-weighted functionals.
For larger (Q), the proof must avoid a (2^Q)-state tuple kernel.  The
intended route is to prove a positive majorant that reduces every occupation
to powers of a small certified functional family.

## Status

- Proved identity: the one-word and ordered-pair laws of the independent
  right-degree sparse map followed by an accumulator.
- Proved reduction: equation (1), including bounded rank-test conditioning.
- Diagnostic: the degree sweep uses nearest long-double arithmetic and
  transfer factors recovered from an existing binary64 receipt.
- Open: outward rounding of the best occupation-one candidate.
- Open: exact transfer-weighted second moment for (Q=2).
- Open: a positive functional majorant covering (3\le Q\le4096).
- Open: replacement of RandomStepConv by the practical RM2Sub inner.

## First finite results

The full occupation-one diagnostic at block size 512 gives:

| degree | raw XORs per constituent | Q1 margin |
|---:|---:|---:|
| 21 | 10,751 | 39.174 bits |
| 23 | 11,775 | 44.826 bits |
| 25 | 12,799 | 50.505 bits |
| 33 | 16,895 | 73.408 bits |

Thus direct weighting lowers the Q1 threshold from degree 33 to degree 23.
This is a real reduction in the proof requirement, but it does not survive
the first reuse test at the same block size.

For two active rows carrying the same local message, the complete
block-size-512 rank-one Chernoff sum gives:

| degree | Q2 rank-one margin |
|---:|---:|
| 23 | 19.174 bits |
| 25 | 23.746 bits |
| 29 | 33.088 bits |
| 33 | 42.707 bits |

The dominant outer weight is 16.  Therefore degree 23 cannot close the
40-bit target through this witness.  This is not a lower bound on the true
bad-code probability; a different transfer inequality could improve it.

Larger rate-half constituents improve this rank-one sum while leaving the
leading raw work approximately (Nr).  The following entries retain only
outer weights at most 128, or at most 64 for the two largest blocks:

| constituent | degree | retained rank-one margin | interpretation |
|---|---:|---:|---|
| [1024,512] | 31 | 43.611 bits | candidate |
| [2048,1024] | 27 | 39.102 bits | this witness misses |
| [2048,1024] | 29 | 44.124 bits | candidate |
| [4096,2048] | 25 | 38.045 bits | this witness misses |
| [4096,2048] | 27 | 43.240 bits | candidate |
| [8192,4096] | 23 | 35.919 bits | this witness misses |
| [8192,4096] | 25 | 41.248 bits | candidate |

The retained sums are partial values of a positive Chernoff union bound.
A value below 40 bits rules out closure with this witness.  A value above 40
bits does not prove closure because the omitted weights and rank-two message
pairs remain positive.  The data show a stable tradeoff: each block-size
doubling has so far reduced the candidate odd degree by about two.
