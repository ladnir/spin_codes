# Goal 04 report: all-occupation closure

## Result

The all-occupation diagnostic closes at the 5% output threshold.  The
floating-point result is

\[
\log_2\mathbb E[\#\{\text{bad nonzero outer words}\}]
\le -19.7044.
\]

The algebra includes every occupation $1\le r\le16384$.  It also includes
every packet support $18\le h\le524352$.  The dominant support interval is
24--31, and its occupation saddle is $r=2$.

This result is not yet a completed distance proof.  The optimization and the
final evaluation use ordinary floating-point arithmetic.  A separate
outward-rounded pass must reproduce a negative exponent.

## Occupation collapse

Goal 03 gives

\[
Z_r(u)
\le
M_1(u)^{r-2}M_{2,*}(u)M_{2,\mathrm{all}}(u)
\]

for every $r\ge2$.  Summing over the choices of active data positions gives

\[
\sum_{r=2}^{K}\binom Kr Z_r(u)
\le
\frac{M_{2,*}(u)M_{2,\mathrm{all}}(u)}{M_1(u)^2}
\left((1+M_1(u))^K-1-KM_1(u)\right).
\]

The one-data contribution is at most $KM_3(u)$.  Thus one scalar expression
covers all nonzero messages.

The implementation checks the stable evaluation of the binomial tail against
direct sums for twelve small instances.  Every check matches within the
recorded floating-point tolerance.

## Reciprocal-binomial secant

The first implementation used the smallest value of $\binom Nh$ throughout
each support interval.  That substitution lost tens of thousands of bits on
wide intervals and produced a useless result.

The successful calculation uses the convex function

\[
f(h):=-\log\binom Nh.
\]

On an interval $[a,b]$, the secant $A+Bh$ satisfies $f(h)\le A+Bh$.
Consequently,

\[
\binom Nh^{-1}\le e^A(e^B)^h.
\]

The extra factor $e^B$ becomes part of the packet-support moment.  This keeps
the exact variation of the inner certificate up to the secant error.

## Target results

| Relative output threshold | Complete first-moment exponent | Dominant support | Occupation saddle |
|---:|---:|---:|---:|
| 5% | -19.7044 | 24--31 | 2 |
| 9% | 8097.1781 | 16384--32767 | 1401 |
| 12% | 213169.0914 | 524335--524352 | 16384 |

At 5%, the leading interval contributions are:

| Packet support | Logarithm base two of contribution |
|---:|---:|
| 24--31 | -20.7484 |
| 32--47 | -21.4252 |
| 18--23 | -22.0175 |
| 48--63 | -26.3306 |
| 64--91 | -31.4294 |
| 92--127 | -41.6746 |
| 128--255 | -43.0269 |

The sum of all refined intervals is $2^{-19.7044}$.

## Interpretation

The 5% result supports the intended proof route.  Neither a joint BCH-weight
profile nor a case analysis over outer occupations is required at this
threshold.  If outward rounding preserves even a small part of the current
margin, the first-moment argument proves the existence of a fixed setup with
binary distance greater than 104870.

The 9% obstruction is not confined to one or two active data blocks.  Its
dominant occupation is approximately 1401.  A targeted joint profile for the
first few shells therefore cannot close the present 9% bound.  Reaching that
threshold requires a stronger bulk interface or a sharper inner analysis.

## Next proof step

The next goal should certify the 5% result.  It must use directed rounding for
the local moments, the occupation tail, the secant parameters, both support
tail optimizations, and the final interval sum.  Floating point may select
the tilts because every admissible tilt remains valid.  The verifier must
evaluate the selected decimal tilts outward.
