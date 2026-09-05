# Goal 03 report: BCH moments for the two parity blocks

## Result

The ordinary BCH weight spectrum supports a useful outer interface.  A joint
BCH-weight profile is not yet necessary for a positive-distance attempt near
the 5% output threshold.

The strongest interface uses Hölder for one active data block.  It uses a
restricted Cauchy inequality for two or more active data blocks.  The latter
keeps the selected data values nonzero on one side of the inequality.  This
removes an artificial zero-value term that dominated the first implementation.

The exact inequalities are combined with the Goal 02 inner envelope.  The
current numerical shell values use floating-point optimization.  They are
diagnostic, not outward-rounded end-to-end certificates.

## Restricted Cauchy interface

For $0<t\le1$, let $Q_x(t)$ be the packet-support moment of the permuted BCH
encoding of $x$.  Define

\[
M_1(t):=\sum_{x\ne0}Q_x(t),
\quad
M_{2,*}(t):=\sum_{x\ne0}Q_x(t)^2,
\quad
M_{2,\mathrm{all}}(t):=\sum_x Q_x(t)^2.
\]

For an outer support containing $r\ge2$ data positions, fix $r-2$ data
values.  The remaining two values determine the parity pair through an
invertible affine map.  Cauchy--Schwarz gives

\[
Z_r(t)
\le
M_1(t)^{r-2}M_{2,*}(t)M_{2,\mathrm{all}}(t).
\]

The earlier bound used $M_{2,\mathrm{all}}(t)^2$.  At the low-support saddle,
that moment was dominated by the zero field element.  Retaining the nonzero
domain on the data side improves the $r=2$, 5% shell by about 20 bits before
the pointwise inner refinement.

For $r=1$, Hölder gives

\[
Z_1(t)
=\sum_{x\ne0}Q_x(t)^2Q_{\alpha x}(t)
\le\sum_{x\ne0}Q_x(t)^3.
\]

Both bounds use only the ordinary BCH spectrum.

## Low-shell results

The table reports the logarithm base two of each complete shell contribution.
The count includes the choice of active data positions.

| Active data blocks $r$ | 5% | 9% | 12% |
|---:|---:|---:|---:|
| 1 | -78.3113 | -43.9241 | -26.0125 |
| 2 | -18.8176 | 5.4298 | 20.6850 |
| 3 | -21.7846 | 18.5900 | 42.5694 |
| 4 | -30.1411 | 31.9049 | 67.3546 |
| 5 | -34.1617 | 41.5337 | 87.5258 |
| 6 | -44.8813 | 55.3708 | 106.7450 |
| 7 | -49.5288 | 62.2216 | 132.5283 |
| 8 | -61.9866 | 83.3712 | 170.4584 |
| 16 | -89.6168 | 172.2629 | 348.2406 |

Every tested shell is below one at the 5% threshold.  The weakest tested
shell is $r=2$, with contribution $2^{-18.8176}$.  The sum of the bounds for
$1\le r\le8$ is $2^{-18.6434}$.

The sampled 5% exponents improve after $r=2$.  The $r=16$ sample continues
that trend.  These samples do not cover all occupations and therefore do not
prove an end-to-end distance.

At 9%, the one-data shell is below $2^{-40}$.  The two-data shell is above
one, and later tested shells grow.  The present spectrum-only compression
does not support a 9% proof.

## Comparison with worst-case parity blocks

Using the data spectrum and only deterministic parity activity is much
weaker.  At 5%, its shell exponents for $r=1,2,3,4$ are

\[
-13.1878,\quad 3.0153,\quad 14.5430,\quad 9.2763.
\]

Restricted Cauchy changes them to

\[
-78.3113,\quad -18.8176,\quad -21.7846,\quad -30.1411.
\]

Thus worst-casing the parity blocks discards too much information.  The
second-moment interface recovers enough information without a joint profile.

## Remaining work

The next proof target should remain at the 5% threshold.  It must bound every
occupation from $r=1$ through 16384.  A continuous or dyadic occupation
envelope should show that the shell exponent keeps decreasing after the
observed low-shell bottleneck.

If that closure succeeds, the construction has a credible route to a small
linear-distance proof using only the ordinary BCH spectrum.  Improving the
threshold toward 9% would require a stronger interface.  Candidate
improvements include support bands, the excluded affine lines in the Cauchy
image, or a targeted joint profile for the first few occupation shells.
