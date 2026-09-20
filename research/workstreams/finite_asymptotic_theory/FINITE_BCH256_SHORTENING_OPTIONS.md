# BCH constituents near length 256

## Decision

The best current constituent is the six-coordinate shortening

\[
  C_{250}:=\operatorname{Shorten}_{\{0,\ldots,5\}}(C_{256}),
\]

where (C_{256}) is the parity extension of the primitive narrow-sense
binary BCH code with parameters ([255,131,\ge37]). The resulting code has
parameters

\[
  C_{250}:[250,125,\ge38].
\]

This option gives up a power-of-two block length but keeps exact rate one
half. It is a genuine shortened BCH code. It does not use a BA layer or an
arbitrary codimension-three subcode.

The exact weight distribution of (C_{250}) remains unknown. The finite
11% proof is therefore not complete. However, the low-occupation calculation
now closes without that distribution. The only material distance obligation
is a high-occupation transfer that keeps the row-weight mixture inside the
250-region recurrence.

## Exact constituent

Let (C_{255}subseteq\mathbb F_2^{255}) be the primitive narrow-sense BCH
code with designed distance 37. Let (C_{256}subseteq\mathbb F_2^{256})
be its parity extension. The standard BCH bound and parity extension give

\[
  C_{256}:[256,131,\ge38].
\]

The local construction script builds an exact parity-check matrix of rank
125. It then imposes (c_i=0) for (i=0,\ldots,5). The six added checks
increase the rank from 125 to 131. Deleting the six zero coordinates preserves
dimension. Thus (C_{250}) has dimension 125 exactly.

Every word of (C_{250}) has even weight. The parent (C_{256}) contains
the all-one word. If a shortened word has weight (w), its lift to
(C_{256}) has six leading zeros. Complementing that lift gives a parent
word of weight (256-w). Consequently every nonzero shortened word satisfies

\[
  38\le w\le218.
  \tag{1}
\]

The receipt `shortened_bch256_parameter_audit.json` records the exact rank
ladders and generator hashes. The script
`audit_shortened_bch250_125.py` regenerates that receipt.

The published parameter evidence is consistent with this construction.
CodeTables records a ([256,128,38]) subcode of the extended
([256,131,38]) BCH code. Wambach proves that the unextended
([255,131,37]) code has true minimum distance 37. The Okayama weight-
distribution index lists the middle-dimension code, but it does not publish
a linked distribution for dimension 131.

- [CodeTables entry for length 256 and dimension 128](https://www.codetables.de/BKLC/BKLC.php?k=128&n=256&q=2)
- [Wambach's minimum-distance result](https://kups.ub.uni-koeln.de/54679/)
- [Okayama weight-distribution index](https://isec.ec.okayama-u.ac.jp/home/kusaka/wd/)

These sources authenticate parameters and minimum distance. They do not
authenticate a weight enumerator for (C_{256}) or (C_{250}).

## Finite construction and probability space

Fix one copy of (C_{250}) and reuse it for every outer block. No outer code
is sampled. For an admissible row count (L), the encoder processes
(125L) message bits as follows.

1. Encode each 125-bit row with the fixed map for (C_{250}).
2. Sample one independent coordinate permutation in (S_{250}) for each row.
3. Transpose the resulting array into 250 regions of length (L).
4. Sample one independent permutation of the (L) positions in each region.
5. Apply the frozen RM2Sub-S19 inner map with independent nonzero multipliers.

The distance probability is over the permutations and RM2Sub multipliers.
The constituent code and message set are fixed. For

\[
  N:=250L,
  \qquad
  d:=\lfloor0.11N\rfloor,
\]

the target statement is

\[
  \Pr[\exists m\ne0:\operatorname{wt}(\operatorname{Enc}(m))\le d]
  \le2^{-\lambda}.
  \tag{2}
\]

The present receipts use (t=128). Therefore (L) must be divisible by
128. The nearest full-rate admissible point to (k=2^{20}) is

\[
  L=8448,
  \quad k_0=125L=1{,}056{,}000,
  \quad N=2{,}112{,}000,
  \quad d=232{,}320.
  \tag{3}

For exactly (k=2^{20}) input bits, pad 7,424 zero bits before the outer
encoder. The padded messages form a subspace of the message space in (3).
Thus any distance bound for the full-rate instance also bounds the padded
instance. The effective rate of the padded instance is approximately
0.496485.

For arbitrary (k), define

\[
  L(k):=128\left\lceil\frac{\lceil k/125\rceil}{128}\right\rceil,
  \qquad N(k):=250L(k).
  \tag{4}

Zero padding to (125L(k)) message coordinates gives the arbitrary-length
wrapper. Both the number of XORs and the permutation work remain (O(N(k)))
for this fixed constituent.

## Comparison of nearby options

| constituent | exact status | rate | consequence |
|---|---:|---:|---|
| deterministic Eq3 subcode ([256,128,\ge38]) | implemented; not a pure BCH code | (1/2) | power-of-two geometry, but unknown spectrum and an all-one exceptional word |
| three-coordinate shortening ([253,128,\ge38]) | exact genuine shortening | (128/253\) | rate exceeds (1/2); reject for an 11% finite proof |
| six-coordinate shortening ([250,125,\ge38]) | exact genuine shortening | (1/2) | recommended finite proof target |
| primitive BCH ([255,123,39]) | true distance 39 is published | (123/255) | more rate slack, but dimension is farther from 128 and the all-one word remains |

At (k=2^{20}), the ([253,128]) option has (L=8192) and
(N=2{,}072{,}576). Even an ideal random-code first-moment benchmark misses
relative distance 11% by approximately 12,113 bits. This failure is caused
by rate, not by rounding or by the BCH spectrum. Puncturing or shortening only
three coordinates is therefore not a useful 11% route.

The power-of-two Eq3 code is an exact ([256,128,\ge38]) construction already
present in the repository. It intersects (C_{256}) with three equality
checks on the first four coordinates. Its exact spectrum is not known. The
three checks also make it a BCH-derived subcode rather than a genuine BCH
constituent.

## Spectrum-free low-occupation bounds

For every binary code of length (n), distance at least 38, and fixed weight
(w), its weight-(w) words form a constant-weight code of Johnson distance
at least 19. A radius-nine packing argument gives

\[
 A_w
 \le
 \left\lfloor
 \frac{\binom{n}{w}}
 {\sum_{i=0}^{9}\binom wi\binom{n-w}{i}}
 \right\rfloor.
 \tag{5}
\]

For (C_{250}), combine (5), the dimension bound (A_w\le2^{125}-1),
even weight, and (1). The receipt
`shortened_bch250_125_constant_weight_envelope.json` contains the resulting
exact integer shell bounds. It is an upper envelope, not a claimed weight
distribution.

Nearest-binary64 RM2Sub calculations give the following results at (3).

- A random-like modeled spectrum gives 59.587 bits of margin for occupation
  one. This value is only a shape diagnostic.
- The exact integer packing envelope gives 13.777 bits of margin for
  occupation one.
- The same envelope gives 19.044 aggregate bits of margin for occupations
  2 through 100.
- Combining those two partial sums gives approximately 13.740 bits.

The first two envelope results use no BCH spectrum assumption. They are not
yet formal certificates because the RM2Sub transfer and logarithmic
optimization use nearest-binary64 arithmetic.

## Dense obstruction

Applying one Bernoulli support envelope independently to every active row is
too loose at high occupation. At (Q=L=8448), the best tested dense tilt
still gives

\[
  \log_2 \text{bound}=446{,}723.720.
  \tag{6}

The failure in (6) is not an outer all-one obstruction. The shortening
removed that word, and every row satisfies (1). The main loss comes from
maximizing the shell likelihood before multiplying the 250 region transfers.

For the packing envelope, the half-Bernoulli comparison costs
174.032 bits per active row. The exact row mass is only

\[
  \log_2|C_{250}\setminus\{0\}|<125.

Replacing the pointwise charge by exact total mass would recover about
414,224 bits at (Q=L). The current transfer would still miss by about
32,500 bits. Thus a total-mass correction alone does not close the proof.

The required lemma must keep the row-weight mixture inside the serialized
250-region transfer. One suitable formulation partitions the permitted
weights into finitely many bands. Let (M_b) denote the unknown number of
nonzero codewords in band (b). The admissible band masses satisfy

\[
  M_b\le U_b,
  \qquad
  \sum_bM_b=2^{125}-1,
  \tag{7}

where (U_b) is the sum of the exact shell bounds in (5). A complete proof
may maximize the coefficient transfer over the polytope (7). This route does
not require the exact BCH spectrum.

## Proof status and next step

The following facts are proved or exactly audited.

- The selected constituent is ([250,125,\ge38]).
- It is even and every nonzero weight lies in ([38,218]).
- The admissible schedule (3) and wrapper (4) are exact.
- The shell multiplicity integers in (5) are valid upper bounds.

The following statements remain numerical diagnostics.

- The 13.777-bit occupation-one margin.
- The 19.044-bit occupation-2-through-100 margin.
- The dense failures and their decomposition.

The remaining distance obligations are:

1. Prove a band-coupled or exact coefficient transfer for every
   (101\le Q\le8448).
2. Convert the low-occupation calculations to outward arithmetic.
3. Sum all occupation bounds and state the resulting value of (lambda) in
   (2).
4. Bind the RM2Sub receipts and permutation law to the implementation.

The next proof task should attack the all-active coefficient problem first.
If the band-mass polytope in (7) closes (Q=L), extend the same witnesses to
occupation intervals. If it does not close, computing the exact BCH spectrum
is unlikely to be the cheapest next step; a lower-rate shortening should be
tested against the same coefficient bound first.
