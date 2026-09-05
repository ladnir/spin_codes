# Goal 02 report: compressed all-support inner envelope

## Result

Goal 02 is complete.  Eighteen intervals cover every integer support

\[
18\le h\le 524352.
\]

For each distance threshold and interval, the receipt stores one output tilt
and one coefficient radius.  These two values define an analytic upper bound
for every support in that interval.  The receipt also gives an
outward-rounded cap for the complete interval.

The calculation concerns the RandomStepConv inner encoder.  It does not use
the outer spectrum.

## Why endpoint checks suffice

Fix a threshold \(D\) and tilts \(0<z,r<1\).  Goal 02 defines

\[
B_{z,r,D}(h):=
\frac{z^{-D}}{\binom Nh}
\frac{r^{-(N-h)}G_z(r)^h}{1-r}.
\]

Goal 01 implies \(p_h(D)\le B_{z,r,D}(h)\).  For fixed tilts,

\[
\log B_{z,r,D}(h)
=C+h\log\!\bigl(rG_z(r)\bigr)-\log\binom Nh.
\]

The binomial coefficients are log-concave.  Hence the displayed function is
convex in \(h\).  Its maximum on an integer interval occurs at an endpoint.

The checker verifies the underlying cross-product identity exactly.  For
\(0\le h\le N-2\), the difference between the two cross-products is

\[
(N-h)(h+2)-(h+1)(N-h-1)=N+1>0.
\]

Thus a floating-point convexity test is not part of the proof.

## Target-size envelope

The table shows selected interval caps.  Each entry bounds the probability
for every support in the stated interval.

| Support interval | 5% threshold | 9% threshold | 12% threshold |
|---:|---:|---:|---:|
| 18--23 | \(2^{-31.2486}\) | \(2^{-16.4648}\) | \(2^{-9.3722}\) |
| 24--31 | \(2^{-49.9085}\) | \(2^{-30.0918}\) | \(2^{-20.5525}\) |
| 32--47 | \(2^{-75.0183}\) | \(2^{-48.4809}\) | \(2^{-35.6718}\) |
| 64--91 | \(2^{-176.8079}\) | \(2^{-123.3241}\) | \(2^{-97.3886}\) |
| 128--255 | \(2^{-383.1053}\) | \(2^{-275.5832}\) | \(2^{-223.2945}\) |
| 1024--2047 | \(2^{-3328.7825}\) | \(2^{-2458.5245}\) | \(2^{-2034.3030}\) |
| 131072--262143 | \(2^{-541513.5671}\) | \(2^{-393315.2559}\) | \(2^{-318775.6501}\) |
| 262144--524352 | \(2^{-1043360.3230}\) | \(2^{-806185.4623}\) | \(2^{-668736.4429}\) |

The weakest interval is 18--23 at all three thresholds.  Every later interval
has a strictly smaller recorded cap.

## Compression loss

The shared-radius formula is weaker than the Goal 01 formula, which assigns a
separate radius to each termination count.  At the 33 Goal 01 anchor points,
the loss ranges from 0.318 to 0.484 bits.

The checker also optimizes an independent certificate at each interval
midpoint.  The largest midpoint loss is 6.45% of the pointwise exponent.
This loss occurs in a wide interval at large support, where the absolute
failure exponent is already hundreds of thousands of bits.  A later outer
calculation can split any interval that requires a tighter local bound.

## Numerical status

Floating-point arithmetic selects the tilts.  Any admissible tilts give a
valid bound.  The verifier converts each selected tilt to an exact decimal
rational.  It evaluates both interval endpoints with 90-digit decimal
arithmetic and directed rounding.

The receipt contains 54 interval certificates: 18 intervals for each of the
three distance thresholds.  It covers 524335 integer supports without a
support-by-support table.

The result is not an end-to-end distance proof.  The next phase must combine
this inner envelope with the exact outer spectrum and the two parity
constraints.
