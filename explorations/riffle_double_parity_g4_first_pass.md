# First proof pass on Riffle DP g=4

## Result

**Riffle DP g=4 survives the obstruction that ruled out Riffle DP g=8.**

This result is a proof gate, not an end-to-end certificate. The exact MDS
shell calculation has the correct entropy slope at packet width four. Exact
enumeration also controls the first three layers of the rare weight-three
tail. The weakest of those layers has about 27 bits of margin.

The next proof step should replace further subset enumeration with a
weighted bound for the remaining weight-three layers. A complete
packet-profile certificate would be premature.

## Candidate under audit

The candidate name is **Riffle DP g=4**. Fix

\[
K=2^{20},\qquad B=K/64=16{,}384,\qquad
\mathbb F=\mathbb F_{2^{64}}.
\]

For data symbols \(m_0,\ldots,m_{B-1}\in\mathbb F\), define

\[
p_0:=\sum_{i=0}^{B-1}m_i,
\qquad
p_1:=\sum_{i=0}^{B-1}\alpha_i m_i,
\]

where \(\alpha_0,\ldots,\alpha_{B-1}\) are distinct and nonzero. Encode each
of the \(B+2\) symbols with the binary extended BCH \([128,64,22]\) code.

This pass keeps the sequential outer layout from the Double-Parity
benchmark. It splits each BCH word into 32 contiguous four-bit packets. A
uniform permutation acts on all packets before the recursive inner map.
No striped layout or additional mixer is part of this candidate.

The unpunctured binary length is

\[
N'=128(B+2)=2{,}097{,}408.
\]

Set \(d:=\lfloor0.09N'\rfloor=188{,}766\). The packet permutation has

\[
M=N'/4=524{,}352
\]

positions. The final \(T:=\lfloor d/64\rfloor=2{,}949\) inner nodes contain

\[
L=16T=47{,}184
\]

packet positions.

## Exact terminal calculation

Fix an outer codeword with \(r\) nonzero four-bit packets. Under the uniform
packet permutation, all active packets enter the final \(T\) nodes with
probability

\[
P_4(r):=\frac{\binom{L}{r}}{\binom{M}{r}}.
\]

On this event, the inner output has weight at most \(64T\le d\). Each outer
codeword therefore contributes at least \(P_4(r)\) to the first moment.

The Double-Parity outer is an \([B+2,B,3]\) MDS code over \(\mathbb F\). Let
\(A_w\) denote its exact number of words with block weight \(w\). Every
active BCH word occupies at most 32 packets. Hence

\[
\mathbb E[Z_d]\ge A_wP_4(32w)
\]

for each \(w\ge3\). Exact arithmetic gives:

| Outer block weight \(w\) | \(\log_2(A_wP_4(32w))\) |
|---:|---:|
| 3 | -230.231606973875 |
| 4 | -265.504413728495 |
| 5 | -301.127813132726 |
| 6 | -337.042932694766 |
| 7 | -373.209151656445 |
| 8 | -409.596743929833 |
| 9 | -446.183010723088 |
| 10 | -482.950051442122 |

These lower bounds decrease with \(w\). At \(g=8\), the corresponding
values increased by about 20 bits per block and crossed the target at
\(w=5\). The prior positive entropy slope is therefore absent at \(g=4\).

The change has a simple scale explanation. The terminal fraction is

\[
\rho:=L/M\approx0.089985.
\]

Placing all 32 packets from one BCH word in the terminal region costs about

\[
32\log_2(1/\rho)\approx111.17
\]

bits. One nonzero field value supplies about 64 bits of value entropy. At
\(g=8\), the analogous placement cost is only about 55.59 bits.

The table contains lower bounds. Its decrease does not upper-bound the full
first moment. BCH words that occupy much fewer than 32 packets have larger
terminal probabilities and can dominate the sum.

## Exact local tail

The existing exact byte-support enumeration contains every nonzero BCH word
on at most eight bytes. It contains 12,587 words. Their four-bit packet
supports are:

| Four-bit packets | Count |
|---:|---:|
| 11 | 5 |
| 12 | 89 |
| 13 | 607 |
| 14 | 2,236 |
| 15 | 5,018 |
| 16 | 4,632 |

Any word on at most eight nibbles would also occupy at most eight bytes.
This calculation originally proved a lower bound of nine on nibble support.

A new enumerator works directly with the 32 four-bit packets. For each packet
subset, it computes the kernel of the restricted parity-check matrix. It then
authenticates every emitted word by re-encoding its recovered message.

The exact scan through support 13 visited 809,785,132 packet subsets. It gave:

| Four-bit packets | Exact count |
|---:|---:|
| 1--10 | 0 |
| 11 | 20 |
| 12 | 1,526 |
| 13 | 37,014 |

The maximum shortened dimension was two. The local BCH code therefore has
exact nibble distance 11. The earlier byte enumeration contains five of the
20 support-11 words, which provides an independent cross-check.

## Exact minimum outer tail

Every block-weight-three outer word contains three nonzero BCH words. Its
total nibble support is therefore at least 33. Equality requires three local
words of support 11.

The 20 support-11 messages permit exact enumeration of all four outer
categories. The result is:

| Category | Words of total support 33 |
|---|---:|
| one data position | 26 |
| two data positions, only \(p_0\) active | 0 |
| two data positions, only \(p_1\) active | 0 |
| three data positions | 0 |

The exact terminal contribution of these 26 words is

\[
26P_4(33)=2^{-109.961733072748}.
\]

The next two layers admit rigorous upper bounds from the exact local counts.
For a fixed ordered value triple, at most

\[
\binom B2+2B+1
\]

outer supports can realize that triple across the four categories. Applying
this bound to every ordered local triple gives:

| Total nibble support | Bound type | Terminal contribution |
|---:|---|---:|
| 33 | exact | \(2^{-109.961733072748}\) |
| 34 | upper bound | \(2^{-70.332635254404}\) |
| 35 | upper bound | \(2^{-67.155898463807}\) |

Thus each authenticated low layer lies below \(2^{-40}\). These rows do not
bound total support 36 or greater.

The obvious completion is too loose. For support \(s\ge14\), upper-bound the
nibble count by every BCH word whose byte support lies between
\(\lceil s/2\rceil\) and \(\min(s,16)\). Contracting this envelope across
three local words gives an aggregate upper exponent of about \(+55.19\).
This value rejects the envelope, not the construction. The next bound must
retain more nibble information.

## Weight-three sampling

Block weight three is the minimum MDS shell. Its supports fall into four
categories: one data position, two data positions with \(p_0=0\), two data
positions with \(p_1=0\), and three data positions with both parities zero.

The diagnostic sampled 300,000 words independently in each category. It
used the deterministic coefficients \(\alpha_i=x^i\). For each sample, it
computed the exact BCH codewords and their total nibble support.

| Category | Smallest sampled support | Mean support | Estimated terminal contribution |
|---|---:|---:|---:|
| one data position | 70 | 90.0033 | \(2^{-183.45}\) |
| two data positions, only \(p_0\) active | 76 | 89.9951 | \(2^{-190.99}\) |
| two data positions, only \(p_1\) active | 72 | 90.0031 | \(2^{-177.09}\) |
| three data positions | 75 | 89.9995 | \(2^{-175.38}\) |

The estimates have large apparent margin. They do not control rare events.
The g=8 investigation already showed that a small algebraic tail can exceed
a Monte Carlo estimate.

## Proof status and next gate

The current evidence supports four statements.

1. The MDS shell obstruction that rejected Riffle DP g=8 does not transfer
   to packet width four.
2. Typical block-weight-three words have ample terminal-placement margin.
3. The exact minimum layer and the next two bounded layers have ample margin.
4. No end-to-end \(2^{-40}\) bound has been proved for Riffle DP g=4.

The next gate is a nibble-aware weighted upper bound for all remaining
block-weight-three layers. Brute-force packet-subset scans grow too quickly
after support 13, and the byte-support envelope is insufficient. The bound
should combine local packet-support information with the four exact MDS
support categories. The candidate should advance to a larger certificate
only if the complete weight-three shell remains below \(2^{-40}\).

The reproducible tools are:

- `scripts/enumerate_bch_g4_low_support.cpp`;
- `scripts/analyze_riffle_double_parity_g4.py`.

On Windows, the exact low-tail receipt can be reproduced with:

```powershell
cl /nologo /std:c++20 /O2 /EHsc /W4 `
  /Fe:enumerate_bch_g4_low_support.exe `
  scripts\enumerate_bch_g4_low_support.cpp
.\enumerate_bch_g4_low_support.exe `
  --max-support 13 --output bch_g4_le13.bin
python scripts\analyze_riffle_double_parity_g4.py `
  --skip-local-tail `
  --exact-low-message-file bch_g4_le13.bin `
  --exact-low-max-support 13
```
