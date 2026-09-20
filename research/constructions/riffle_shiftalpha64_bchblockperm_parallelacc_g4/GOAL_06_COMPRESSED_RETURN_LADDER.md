# Goal 06: compressed return ladder

## Target

Replace the exact zero-return coordinate in the Goal 05 kernel by a scalar
tilt, reproduce the three- and four-block reference calculations, and then
increase the number of active minimum-weight blocks until a new obstruction
or a new compression opportunity appears.

The output threshold remains

\[
D=188{,}766=\lfloor 0.09L\rfloor,
\qquad L=2{,}097{,}408.
\]

This goal studies the inner kernel for fixed block profiles. It does not yet
sum the outer spectrum.

## First compression: zero returns

Let (H) be the number of active packets and let (m) be the number of
active prefixes at accumulator state zero. The exact gap sum contains

\[
\binom{N-H+m}{m}.
\]

For every (0<\theta<1), coefficient extraction from
((1-z)^{-(N-H+1)}) gives

\[
\binom{N-H+m}{m}
\le
\theta^{-m}(1-\theta)^{-(N-H+1)}.
\tag{1}
\]

Thus every return to zero receives the multiplicative weight
\(\theta^{-1}\). The dynamic program no longer retains (m). It retains the
total input weight, packet support (H), and the current state in
\(\mathbb F_2^4\).

At the optimized tilt values, (1) gives:

| profile | exact Goal 05 bound | return-compressed bound | loss |
|:---|---:|---:|---:|
| ((24,24,24)) | -40.2950 bits | -36.7901 bits | 3.5049 bits |
| ((22,22,22,22)) | -44.6721 bits | -43.7832 bits | 0.8889 bits |

The return coordinate is therefore removable at modest proof cost.

## Second compression: total binary weight

For a polynomial (F(x)) with nonnegative coefficients and every (x>0),

\[
[x^W]F(x)\le x^{-W}F(x).
\tag{2}
\]

Applying (2) to the packet-value enumerator removes the exact total binary
weight. The remaining dynamic-program state is only

\[
(H,q),\qquad q\in\mathbb F_2^4.
\]

The optimized double-compressed reference bounds are -32.4042 bits for
((24,24,24)) and -39.1982 bits for ((22,22,22,22)). Relative to the
return-compressed values, the second envelope costs 4.3859 and 4.5850 bits,
respectively. The two compressions together cost 7.8908 and 5.4740 bits
relative to the exact Goal 05 sums.

There is no additional equal-block balance loss after (2). If there are
(s) blocks of weight (h), then

\[
\frac{1}{\binom{128s}{sh}}
\frac{\binom{128s}{sh}}{\binom{128}{h}^{s}}
=
\frac{1}{\binom{128}{h}^{s}}.
\tag{3}
\]

The receipt continues to display the two factors separately for comparison
with Goal 05, but (3) shows that the reported `conditioning_penalty_bits` is
bookkeeping, not an independent proof penalty.

## Minimum-weight ladder

The following fixed tilts were used for the climb:

\[
x=0.334,
\]

with positive-cost and zero-return parameters extrapolated from the optimized
three- and four-block calculations. Each row is a rigorous bound even though
the fixed parameters need not be optimal.

| active weight-22 blocks (s) | bound on (K_D(22^s)) | two-array state |
|---:|---:|---:|
| 4 | (2^{-39.1982}) | 32.3 KiB |
| 5 | (2^{-49.8750}) | 40.3 KiB |
| 6 | (2^{-60.4130}) | 48.3 KiB |
| 8 | (2^{-82.1280}) | 64.3 KiB |
| 12 | (2^{-131.7594}) | 96.3 KiB |
| 16 | (2^{-188.9353}) | 128.3 KiB |
| 24 | (2^{-317.6716}) | 192.3 KiB |
| 32 | (2^{-458.2976}) | 256.3 KiB |
| 64 | (2^{-1073.4895}) | 512.3 KiB |

The kernel exponent strengthens throughout the tested range. From 32 to 64
blocks it gains 615.2 bits. The current evidence therefore does not show a
large-(s) inner-kernel obstruction for the equal minimum-weight profile.

## What was compressed, and what remains

The original exact calculation indexed four quantities:

\[
(W,H,m,q).
\]

Goal 06 reduces this to

\[
(H,q).
\]

A direct four-coordinate table has storage proportional to (s^3). After the
first compression, the ((W,H,q)) table has quadratic storage and roughly
cubic packet-by-packet running time. The second compression reduces these to
linear storage and roughly quadratic running time. The 64-block calculation
completes in a few seconds on the reference machine.

The next computational compression target is (H). It should not be removed
yet by a crude worst-case maximization because the global placement factor

\[
\frac{(1-\theta)^{-(N-H+1)}}{\binom NH}
\tag{4}
\]

depends materially on (H). A promising exact representation is

\[
\frac{1}{\binom NH}
=(N+1)\int_0^1 t^H(1-t)^{N-H}\,dt.
\tag{5}
\]

Equation (5) converts the support sum into a one-dimensional integral of a
16-state transfer operator with a scalar support fugacity. It could permit
matrix powering while preserving the global placement law. This becomes
worth implementing when the ladder or the outer-spectrum sum requires much
larger (s).

## Reproduction

Reference optimized calculations:

```powershell
python scripts/analyze_riffle_shiftalpha64_compressed_return.py --part-count 3 --part-weight 24 --compress-weight --optimize --optimizer-maxiter 50
python scripts/analyze_riffle_shiftalpha64_compressed_return.py --part-count 4 --part-weight 22 --compress-weight --optimize --optimizer-maxiter 50
```

The 64-block fixed-tilt endpoint:

```powershell
python scripts/analyze_riffle_shiftalpha64_compressed_return.py --part-count 64 --part-weight 22 --distance 188766 --scaled-cost 585 --zero-scale 335.7 --weight-tilt 0.334 --compress-weight
```

The summary receipt is
`receipts/goal06_compressed_return_ladder.json`.

The analysis script has SHA-256
`D1EB63EEBA9B91378D464FACFD436300983C6FA089E3DC98FB19F2571CBB4B41`.

## Next goal

Connect this inner-kernel ladder to a coarse outer block-occupation spectrum.
Begin with dyadic occupation bins and keep equal-weight 22 as the pessimistic
inner profile inside each bin. At every bin enlargement, compare the added
outer multiplicity with the measured gain in the inner exponent. Implement
the support-integral compression in (5) only when the retained (H)
coordinate becomes the limiting cost.
