# Exact outer definition for the finite \(k=2^{20}\) certificate

## Object being defined

This document defines the current power-of-two outer map. It separates the
deterministic repeated code from the setup-time routing. RandomStepConv is
the inner map and is not part of this definition. The earlier certificate
with 176 zero rows is a historical, different-length instance.

Set

\[
 K=64,\qquad B=128,\qquad L=16384.
\]

The message dimension and routed outer length are

\[
 k=KL=1048576,
 \qquad
 N=BL=2097152=2^{21}.
\]

## The one fixed constituent

All polynomial arithmetic in this section is over \(\mathbb F_2\).  Let

\[
 g(X)=\sum_{e\in E}X^e,
\]

where

\[
\begin{split}
 E=\{&0,1,2,3,4,9,11,13,19,20,22,24,27,28,29,31,\\
     &35,36,40,42,44,46,50,55,58,60,61,62,63\}.
\end{split}
\]

Equivalently, when bit \(e\) is the coefficient of \(X^e\),

\[
 g=\mathtt{0xf4845518b9582a1f}.
\]

For an input block \(u=(u_0,\ldots,u_{63})\), put

\[
 u(X)=\sum_{i=0}^{63}u_iX^i,
 \qquad
 p(X)=u(X)g(X).
\]

Since \(\deg p\le126\), this multiplication requires no reduction.  Define
the constituent encoder \(E_{\rm BCH}:\mathbb F_2^{64}\to\mathbb F_2^{128}\)
by

\[
 E_{\rm BCH}(u)_j=[X^j]p(X)\quad(0\le j<127),
 \qquad
 E_{\rm BCH}(u)_{127}=\sum_{j=0}^{126}[X^j]p(X).
 \tag{1}
\]

Thus coordinate 127 is the extension-parity coordinate.  An equivalent
integer implementation is

\[
 E_{\rm BCH}(u)
 =\bigoplus_{i:u_i=1}\left((g\mathbin{\mathtt{<<}}i)
       \mathbin{\mathtt{|}}2^{127}\right).
 \tag{2}
\]

The 64 generator rows in (2), serialized in increasing row order as
16-byte little-endian integers, have SHA-256 digest

```text
972cfc1c6de12e4ddc0c67680fd8e3cddedbd84d0a409c7061e09a86b818ca56
```

This is the extended primitive narrow-sense BCH code with parameters
\([128,64,22]\).  One conventional field presentation uses
\(\mathbb F_2[T]/(T^7+T+1)\), primitive polynomial `0x83`, and primitive
designed distance 21 before adding the extension-parity coordinate.  The
polynomial and coordinate rules in (1) are the normative encoder definition;
the field presentation is descriptive.

The same map \(E_{\rm BCH}\) is used in every row.  No BCH generator, BCH
code, or constituent permutation is sampled as a new code.

## Repetition

Write the message as consecutive blocks

\[
 m=(u_0,\ldots,u_{L-1}),
 \qquad u_i\in\mathbb F_2^{64}.
\]

The bit \(m_{64i+a}\) is bit \(a\) of \(u_i\).

For every \(0\le i<L\), compute

\[
 c_i=E_{\rm BCH}(u_i)\in\mathbb F_2^{128}.
\]

Before routing, the outer map is therefore the direct sum of 16,384 copies
of one fixed constituent. There are no zero padding rows.

## Setup-time routing

The proof ensemble samples the following mutually independent objects once.

1. For each row \(i\in\{0,\ldots,L-1\}\), sample
   \(\pi_i\) uniformly from the symmetric group \(S_{128}\).
2. For each coordinate region \(j\in\{0,\ldots,127\}\), sample
   \(\sigma_j\) uniformly from \(S_L\).

Use the output-to-source convention

\[
 a_{i,j}=c_{i,\pi_i(j)},
 \qquad
 x_{jL+t}=a_{\sigma_j(t),j}
 \quad(0\le j<128,\ 0\le t<L).
 \tag{3}
\]

Equation (3) first permutes the 128 coordinates independently within each
row, then transposes the rows into 128 consecutive regions, and finally
permutes the \(L\) positions independently within each region.  The routed
outer output is

\[
 x=(x_0,\ldots,x_{N-1})\in\mathbb F_2^N.
\]

This \(x\) is the input sequence consumed by the inner encoder.  Replacing
every permutation in (3) by its inverse gives the same setup distribution,
but (3) fixes one convention for implementations and proof statements.

The permutations are sampled once and reused for all messages.  Conditional
on them, the outer map is deterministic and linear.  They collectively form
a permutation of the \(BL\) coordinates, so the routed outer code has
dimension \(2^{20}\) and outer minimum distance 22.

## What is random and what is not

The constituent polynomial and all 16,384 uses of its encoder are
deterministic. Only the permutations \((\pi_i)_i\) and \((\sigma_j)_j\) are
outer-side setup randomness. The RandomStepConv matrices are independent
inner-side setup randomness.

The parity pivots used in the distance proof are auxiliary random variables
introduced to represent a uniform even reference row.  They are not sampled
or stored by the encoder and do not change (1)--(3).

## Spectrum premise used by the certificate

The distance proof uses the complete ordinary weight enumerator \((A_w)\)
stored in `scripts/EBCH128_64.wd`.  That file has SHA-256 digest

```text
f633eb9a2f76c64c1d74313b066c2613adba642162e04b67466714f3c30c4b8a
```

It states \(A_0=A_{128}=1\), \(A_w=0\) for unlisted weights, minimum
nonzero weight 22, and total mass \(\sum_w A_w=2^{64}\).  The finite-distance
verifier authenticates and consumes this table.

There is one remaining provenance obligation.  The current local certificate
does not derive the complete table from the concrete generator (1), nor does
it contain a separate formal proof that the table is the weight enumerator of
that row space.  Therefore the fully concrete interpretation of the distance
theorem uses the following explicit premise:

> The row space generated by (2) has the weight enumerator in
> `scripts/EBCH128_64.wd`.

The repository independently reconstructs the BCH generator and checks its
row digest, and the table is the imported published enumerator for the stated
extended BCH parameters.  Those facts are strong identification evidence,
but they are not yet a local derivation of the spectrum-generator binding.

## Compact interface

For setup

\[
 \rho=((\pi_i)_{i=0}^{L-1},(\sigma_j)_{j=0}^{127}),
\]

the exact outer interface is

\[
 \operatorname{Outer}_\rho:\mathbb F_2^{1048576}
 \longrightarrow\mathbb F_2^{2097152},
 \qquad
 m\longmapsto x,
\]

where (1) defines the sole constituent and (3) defines every output
coordinate.
