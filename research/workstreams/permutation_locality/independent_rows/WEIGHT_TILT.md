# Coupling the inner bound to total input weight

The transfer bound maximizes over packet-weight patterns. Without a constraint
on total input weight, different steps can select expensive patterns that are
not representative of the outer words being counted. A weight tilt couples
these two parts of the argument. It changes the bound, not the code or its
setup distribution.

## Cancellation of the extra factor

Fix an active group of four BCH rows. Let `W` be their total Hamming weight,
`U` the size of their union after independent row shuffles, and `J` the size
of their common intersection. Thus `J` counts the regions receiving a full
four-bit packet. Packet routing and lane permutations preserve `W` and `J`.

Choose `a > 0` and `0 < rho <= 1`. For an inner step with nonzero packet
weights `b_1,...,b_j`, multiply its transfer coefficients by

\[
 a^{b_1+\cdots+b_j}\rho^{\#\{i:b_i=4\}}.
\]

The factor is applied before maximizing over packet-weight patterns. Across
all steps, the factors multiply to `a^(sum_g W_g) rho^(sum_g J_g)`. The outer
count therefore uses the reciprocal weight

\[
 a^{-W}\gamma^J,\qquad \gamma:=\rho^{-1}.
\]

There is no additional normalization factor. The usual Chernoff factor for
the output-weight threshold and the coefficient-extraction denominators are
unchanged. The all-zero four-row tuple has `W=J=0`, so its weighted mass is
still one; it is subtracted once when counting active groups.

## Exact averaged outer counts

Let `A_w` count the words of weight `w` in one BCH row. Independent uniform
coordinate permutations assign each support of size `w` mass
`A_w / binom(256,w)`. Replacing the row spectrum by

\[
 \widetilde A_w:=A_w a^{-w}
\]

incorporates the total-weight factor into the existing intersection-weighted
enumeration. Its support shell at `u` is exactly

\[
 \sum_{w_1,\ldots,w_4}
   \left(\prod_{i=1}^4 A_{w_i}\right)a^{-\sum_i w_i}
   \mathbb E\left[\gamma^J\mathbf 1\{U=u\}\mid w_1,\ldots,w_4\right].
\]

The expectation is over the four independent uniform subsets of the stated
sizes. This is a weighted count, not a probability distribution.

The implementation clears a common denominator `D` from the rational row
spectrum, runs the existing exact integer inversion, and divides by `D^4`.
Componentwise upper bounds on `A_w` yield shell upper bounds because the
displayed sum has nonnegative coefficients. This justification applies to
the original positive sum, not to arbitrary substitutions into an alternating
inversion formula.

## The known-total CDF bound

The additional identity `sum_w A_w = 2^128` can improve cumulative counts.
Fill the shell caps from the smallest weight upward until this total is
reached. The resulting row-weight distribution is stochastically no larger
than the true one.

First suppose `a >= 1`.
Couple these row weights monotonically and use prefixes of a common random
coordinate ordering for their supports. Reducing the row weights can only
reduce `W` and `U`. Consequently

\[
 a^{-W}\mathbf 1\{U\le u\}
\]

can only increase under this coupling. The greedy spectrum therefore bounds
its cumulative weighted count. Since `J <= U <= u`, multiplying this bound
by `gamma^u` bounds the joint weighted count. The verifier takes the minimum
of this bound and the cumulative sum of the componentwise shell caps.

The greedy spectrum must not be substituted directly into the joint factor
`a^(-W) gamma^J`: reducing row weights can also reduce `J`. An exact small
counterexample is included in the tests. With length two, two rows,
`a=2`, `gamma=8`, actual spectrum `[1,0,1]`, and caps `[1,1,1]`, the nonzero
joint count is `9/2`, whereas direct substitution of the greedy spectrum
`[1,1,0]` gives only `17/8`.

For `0 < a < 1`, even `a^(-W)` increases with `W`, so the preceding weighted
coupling does not apply. Use the unweighted known-total union CDF `C_u`
instead. On `U <= u`, the inequalities `W <= 4u` and `J <= u` imply

\[
 \sum_{\text{nonzero tuples}}\mathbb E\left[
 a^{-W}\gamma^J\mathbf 1\{U\le u\}\right]
 \le C_u(\gamma a^{-4})^u.
\]

Take the minimum with the cumulative sum of the exact joint shell caps.
This bound is intentionally conservative; exact shell caps still apply
without modification. The tests include a separate failure example for
naive greedy substitution below one: with the same spectra, two rows,
`a=1/2`, and `gamma=1`, the true nonzero weighted count is `24`, but the
greedy count is only `8`.

## Implementation review

The independent-row driver passes the same `a` into the base transfer bound,
fresh and window moments, multi-window averaging, collision and zero-return
bounds, mature-tail lift, histogram refinement, full-feedback refinement, and
joint-cancellation refinement. Its direct tail checks also include `a`.
Unweighted census tables remain reusable: their refinements apply the factor
for each packet-weight pattern.

The two-update transformation operates on coefficients that already include
the factor. Its only explicit subtraction concerns an empty step, whose
input-weight factor is one. No further input-weight argument is required.

The operator-cache key includes `a`, the output tilt, `rho`, occupancy degree,
precision, refinement choices, and the fixed two-update choice. Its source
digest covers the parent Python sources, the driver, the cache implementation,
and both map files read by the inner census. The outer counts are recomputed,
not loaded from this cache. Entries store exact dyadic upper endpoints with
a checksum; they are local memoization, not independent proof certificates.
Generator sources must stay unchanged from process startup through cache
generation: the digest describes files, not Python's already imported code.
Final independent reproduction should omit the cache.

Validation currently includes 17 exact support tests, including rational
tilts `a=1/2` and `a=49/50`, a cache round-trip and
corruption test, and ten distinct cache keys after changing mathematical or
source parameters. At `a=1`, the new outer functions reproduce the previous
bounds. None of these checks fills an uncovered support class or occupancy;
only a complete outward cover establishes its stated occupancy bound.
