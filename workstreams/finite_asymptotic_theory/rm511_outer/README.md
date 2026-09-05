# RM(5,11) Outer Route

Status: active proof investigation.

## Construction under test

The outer code is the direct sum of 1024 copies of one fixed binary
Reed--Muller constituent

\[
  \operatorname{RM}(5,11)=[2048,1024,64].
\]

The 1024 input coordinates index monomials in eleven variables of degree at
most five.  The 2048 output coordinates are evaluations at the points of
\(\mathbb F_2^{11}\), in a fixed order.  The outer has no setup randomness.

The first transfer keeps the existing target

\[
  k=2^{20},\qquad N=2^{21},\qquad D=228{,}590,
\]

and the uniform-routing, RandomStepConv-M22 probability space.  Each outer
row receives an independent uniform permutation of its 2048 coordinates.
Each transposed region receives an independent uniform permutation of the
1024 row positions.  The inner independently samples one uniform binary
\(23\)-by-\(23\) linear map at each output position.

## Exact low-weight information

Kasami and Tokura give the complete enumerator below twice the minimum
distance.  The resulting nonzero RM(5,11) coefficients below 128 are

| weight | exact multiplicity |
|---:|---:|
| 64 | 113,562,778,208 |
| 96 | 36,668,966,832,250,368 |
| 112 | 54,200,225,260,621,496,320 |
| 120 | 1,213,348,790,503,877,902,336 |
| 124 | 3,841,018,795,264,309,198,848 |

`verify_rm511_low_weight_formulas.py` evaluates the integer formulas.  As an
independent regression test, it reproduces all four corresponding nonzero
coefficients of the authenticated RM(4,9) spectrum exactly.

The 1976 Kasami--Tokura--Azumi result extends exact enumeration through
weights below \(2.5d=160\).  Its longer formulas have not yet been transcribed
and checked here.  No sampled estimate will be used as a certificate value.

## First transfer result

A binary64 witness search gives the following shell margins at M22 and the
10.9% target:

| weight | margin |
|---:|---:|
| 64 | 83.3743 bits |
| 96 | 138.7030 bits |
| 112 | 166.4248 bits |
| 120 | 181.3770 bits |
| 124 | 189.5777 bits |

The outward replay proves that the aggregate over these five shells is at
most

\[
  2^{-83.3743314000}.
\]

It is dominated by weight 64.  The receipt is
`rm511_q1_low_randomstepconv_M22_d109_outward.json`.  This is a certificate
for the stated partial sum, not for all of occupation one.

A total-mass bound for every word of weight at least 128 is far too weak: it
assigns all \(2^{1024}\) words to weight 128 and gives a positive log bound.
Therefore, the unknown bulk cannot be discarded without a spectrum envelope.

## Proof route

RM(5,11) is a Type-II self-dual code.  Its weight enumerator therefore has a
Gleason expansion

\[
  W(x,y)=\sum_{j=0}^{85}c_j\Phi(x,y)^{256-3j}\Psi(x,y)^j,
\]

where

\[
  \Phi=x^8+14x^4y^4+y^8,
  \qquad
  \Psi=x^4y^4(x^4-y^4)^4.
\]

The exact coefficients below weight 128 determine \(c_0,\ldots,c_{31}\).
A numerical linear-programming probe over the remaining nonnegative Gleason
enumerators did not give a useful envelope: the relaxation can place nearly
all code mass at weight 128.  This rejects the generic Type-II relaxation as
the next certificate mechanism.  It does not assert that the true RM
enumerator has such a coefficient.

A distance-only Johnson-space packing envelope also fails.  At weight 128,
the packing cap exceeds the cap needed for a 40-bit pointwise transfer by
268.01 bits.  At weights 160, 192, 224, 256, and 320, the corresponding
excesses are 301.07, 320.66, 327.67, 322.29, and 227.52 bits.  At weight 384,
the total-mass bound is already sufficient.  The diagnostic receipt is
`rm511_q1_packing_gap_diagnostic.json`.

The viable proof target is therefore an RM-specific spectrum or cumulative
weight bound for weights 128 through 380.  The exact
Kasami--Tokura--Azumi coefficients below 160 would reduce the first part of
this interval, but they do not by themselves close it.

## Open obligations

1. Obtain a finite, explicit RM-specific upper envelope for weights 128
   through 380.
2. Replay the resulting complete Q1 sum with outward arithmetic.
3. Certify Q2 using the same envelope.
4. Build an RM-compatible high-occupation transfer.
5. Only after proof viability, benchmark the length-2048 Möbius transform.
