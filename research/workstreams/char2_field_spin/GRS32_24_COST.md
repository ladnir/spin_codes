# A fixed parity circuit for systematic GRS[32,24]

The verified circuit below computes the parity contribution with 34 fixed GF256 products and 114 Boolean instructions per four source rows.
It needs another 24 XOR instructions to add the systematic contribution.
This is an arithmetic cost result for the outer map; it supplies no distance bound for the complete routed encoder and no native timing.

## Outer map and additive cosets

Use the AES field

\[
F=\mathbb F_2[X]/(X^8+X^4+X^3+X+1).
\]

A byte denotes its polynomial-basis field element, so addition is bytewise XOR.
Let \(U=\{0,\ldots,7\}\), an eight-element additive subspace.
The four cosets are \(U\), \(8+U\), \(16+U\), and \(24+U\).
They are not cosets of a GF8 subfield: GF256 has no such subfield.

For source values \(m_0,\ldots,m_{23}\), let \(f\in F[Z]\) be the unique polynomial of degree below 24 with \(f(i)=m_i\).
The systematic outer returns

\[
(m_0,\ldots,m_{23},f(24),\ldots,f(31)).
\]

Write its parity matrix as \(P\in F^{8\times24}\), where

\[
P_{r,i}=\prod_{\substack{0\le j<24\\j\ne i}}
             \frac{24+r+j}{i+j}.
\]

All additions in this formula are field additions, including those between byte labels.
For the eight parity inputs \(x\in F^8\), the transpose computes \(P^Tx\in F^{24}\).
The 24 systematic inputs are added afterward.

Partition the 24 outputs into three consecutive blocks of eight.
Each block is symmetric XOR-circulant: its entry at row \(i\), column \(j\) is \(c_{i\mathbin{\mathrm{xor}}j}\).
Translation by an element of \(U\) preserves the source and parity point sets, which gives this symmetry.
The independently interpolated first rows are:

| Source block | First row, as AES-field bytes |
|---|---|
| 0 through 7 | 208, 230, 192, 134, 42, 68, 78, 81 |
| 8 through 15 | 162, 127, 147, 40, 38, 203, 227, 105 |
| 16 through 23 | 52, 172, 158, 17, 209, 32, 250, 29 |

The XOR of each row is one.

## Transform and shared products

Index the eight coordinates by subsets of \(\{0,1,2\}\).
For any vector \(a\), define the superset transform

\[
(Ma)_S=\sum_{T\supseteq S}a_T.
\]

The transform costs 12 XORs and satisfies \(M^2=I\) in characteristic two.
For one circulant block, put \(\gamma=Mc\) and \(z=Mx\).
Its output is \(Mr\), where

\[
r_S=\sum_{I\subseteq S}\gamma_I z_{S\setminus I}.
\]

This is multiplication in
\(F[u_0,u_1,u_2]/(u_0^2,u_1^2,u_2^2)\), after substituting \(u_i=g_i+1\) in the group algebra.
Here \(\gamma_\varnothing=1\), so the eight products by that coefficient are free.
One block has 19 remaining scalar products; three independent blocks would have 57.

For each fixed \(z_S\), the circuit shares products between the three blocks.
Multiplication is linear over GF2, so \((a+b)z_S=az_S+bz_S\).
The reference contains the complete fixed coefficient bases and XOR instructions.
Their exact costs are:

| Input subset mask | Fixed products | XORs constructing other scalar products |
|---:|---:|---:|
| 0 | 7 | 18 |
| 1 | 7 | 2 |
| 2 | 7 | 2 |
| 3 | 2 | 1 |
| 4 | 7 | 2 |
| 5 | 2 | 1 |
| 6 | 2 | 1 |
| 7 | 0 | 0 |
| Total | 34 | 27 |

The input transform costs 12 binary XORs.
The three output transforms cost 36 binary XORs.
The three products in the nilpotent basis require 21 binary XORs and 18 three-input XORs to accumulate their terms.
Together with scalar sharing, this gives 96 binary XORs and 18 three-input XORs: 114 Boolean instructions.
Without a three-input XOR instruction, the same circuit uses 132 binary XORs.

## SIMD and memory costs

For a batch of four source rows with 128-bit records, each scalar product becomes one GFNI512 instruction.
Its coefficient is identical across all bytes; the four rows do not require separate coefficient preparation.
The circuit uses 18 distinct nontrivial multiplication constants, occupying 1152 bytes if each is expanded to a full ZMM operand.

The complete outer transpose loads eight parity vectors and 24 systematic vectors, then stores 24 output vectors.
That is 3584 bytes per batch, or \(112/3\) bytes per source symbol, excluding cache write allocation and constant-operand reads.
A standard four-record output transpose adds 48 lane-shuffle instructions per batch.
There are no outer coefficients sampled or streamed separately for each source row.

For comparison, the rectangular Toeplitz 8-by-24 circuit needs 81 variable products and 133 XORs when it shares its input transform.
It has the same record traffic, systematic additions, and output transpose dimensions.
The fixed GRS circuit reduces these arithmetic counts to 34 products and 114 Boolean instructions.
These counts do not establish a runtime ratio: register allocation, constant loads, and intermediate storage remain unmeasured.
In particular, a native implementation must schedule shared products rather than retain every intermediate simultaneously.

All operations use GF256 scalars.
They act componentwise in compatible GF256-module representations of GF(2^m), for \(m\in\{8,32,64,128\}\).
No general multiplication in those larger fields is required.
Input scalars, routing, inner arithmetic, and regional state orientation are outside these outer counts.

## Verification and next gate

The scalar reference constructs \(P\) independently by Lagrange interpolation.
It checks all circulant entries, all transform constants, and every instruction in the scalar-sharing program.
It then checks zero, all 255 nonzero scalar multiples of each of eight coordinate basis vectors, and 256 combined inputs.
The 2297 exact comparisons pass, and the operation counter verifies the stated count on every combined input.
Both maps are GF256-linear, so the coordinate-basis checks also determine equality on every input.

```sh
python -B research/workstreams/char2_field_spin/theory/grs32_24_parity_reference.py --check
```

The next gate is a complete higher-rate or higher-distance analysis of the proposed inner and route.
That analysis must cover multiple active source rows and the dependence created by sampling state orientation once per region.
This cost result does not justify applying a transfer that assumes fresh uniform nonzero state at every step.
