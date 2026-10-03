# Outward local operators for the 24-bit state

The local verifier encloses the ten-state comparison family without an integer Walsh transform over `2^24` entries. It computes profile polynomials exactly, bounds the native transform's rounding error, and finishes the profile arithmetic with rational numbers. The resulting endpoints are conditional on the explicit binary64 arithmetic contract below. They are not a whole-code distance certificate.

## Exact quantities and comparison basis

Fix a dyadic `z=u/v` strictly between zero and one. Let `X_j` occupy exactly `j` of the eight input bytes, chosen uniformly. Each occupied byte is independent and uniform over the 255 nonzero values. For the literal maps `A,C`, define

\[
 W_j(s)=\mathbb E[z^{\operatorname{wt}(X_j)}\mathbf1_{CX_j=s}],
 \qquad M_j(s)=\mathbb E[z^{\operatorname{wt}(X_j+As)}].
\]

The comparison basis consists of the zero state, the uniform distribution on nonzero states, and eight normalized birth distributions. Birth distribution `i` assigns mass `W_i(s)/mu_i` to each nonzero state, where

\[
 \mu_i=\left(\frac{(1+z)^8-1}{255}\right)^i-W_i(0).
\]

The nonzero comparison rows place the entire emission moment in the uniform coordinate. For a nonempty input they also place its quotient by `2^24-1` in the zero coordinate. These are the existing filled upper-comparison operators, not the exact physical transition probabilities. Lower endpoints in this module enclose that comparison family; they are not physical lower bounds.

For a character restricted to one byte with Hamming weight `r`, set

\[
 Q_r=(v+u)^{8-r}(v-u)^r-v^8.
\]

The occupancy-`j` character moment has numerator `[x^j] product_h(1+x Q_{r_h})` and denominator `C(8,j)(255v^8)^j`. The code constructs these integers for each character profile.

For an expansion profile with byte weights `w_h`, set `I_h=u^{w_h}v^{8-w_h}` and `B_h=(v+u)^8-I_h`. The exact emission moment has numerator `[x^j] product_h(I_h+x B_h)` and denominator `C(8,j)255^j v^64`. These formulas also support the exhaustive smaller-field tests, with eight replaced by the packet bit width where appropriate.

## Walsh-transform error

The census supplies each state's character profile and expansion profile. Integer comparisons and the literal field maps determine these indices; the census does not use floating arithmetic.

Convert each exact character-profile coefficient `f(t)` to binary64 `fhat(t)`. Exact rational comparisons compute

\[
 \delta=\max_t|\widehat f(t)-f(t)|,
 \qquad R=\max_t|\widehat f(t)|.
\]

Let `n=24`, `S=2^n`, and `gamma_n=n*2^-53/(1-n*2^-53)`. Each Walsh output is a signed sum with exactly `n` rounded additions along every input path. Negation is exact. After exact division by `S`, the uniform absolute error is at most

\[
 \varepsilon=\delta+\gamma_n R.
\]

The verifier checks the exceptional-value conditions before applying this bound. If every stored input is an integer multiple of `2^-q`, each subsequent rounded addition remains on that lattice. The test `q+n<=1022` therefore excludes subnormal nonzero intermediates, including after division by `S`. A separate exact magnitude check bounds every intermediate below overflow. Cancellation does not invalidate the lattice argument.

The transform's result can be slightly negative. Projection onto the nonnegative half-line cannot increase its error from the genuine nonnegative mass. Thus clipping is accompanied by the proven `epsilon`; it is not used to discard an unexplained error.

## Compression and normalization

For each expansion profile `p`, let `N_p` count its nonzero states. The verifier sums the clipped transformed values in chunks. If `bhat_p` is the computed sum and `g=gamma_(S+number_of_chunks)`, exact rational endpoints enclose the genuine profile mass:

\[
 \max\left(0,\frac{\widehat b_p}{1+g}-N_p\varepsilon\right)
 \;\le\;\sum_{s\ne0:\,\operatorname{profile}(As)=p}W_j(s)
 \;\le\;\frac{\widehat b_p}{1-g}+N_p\varepsilon.
\]

The operation count is conservative: a within-chunk bin receives at most one addition per input atom, followed by one merge per chunk. The positive sum stays on the checked dyadic lattice. The transform's magnitude bound also covers these sums.

The zero-state mass receives its own absolute interval. It is exactly zero for occupancies one through three, as established by the checked packet-restriction ranks. Single-packet birth masses are computed directly from their 2040 possible nonzero byte inputs.

To enclose a birth-conditioned emission moment, multiply each profile-mass endpoint by that profile's exact emission numerator. Sum these rational products, then divide by the exact emission denominator and the appropriate endpoint for `mu_i`. This keeps the profile dependence; it does not replace emissions by a global maximum. Uniform-state moments are exact sums against the integer expansion census.

Finally, convert each rational endpoint to binary64 and check the inequality by integer cross-multiplication through `Fraction`. Stored floats are exact dyadic endpoints. Exact zeros remain zero. The receipt includes their hexadecimal representations and the exact numerator, denominator, and hexadecimal representation of `z` when `z` itself is representable.

## Arithmetic contract and use

The native circuit assumes IEEE binary64 round-to-nearest, correctly rounded NumPy addition, subtraction, negation, and power-of-two scaling, and no fast-math reassociation. Each weighted `bincount` must use at most one binary64 addition per input atom. The checked range argument excludes subnormal intermediates and overflow. Runtime exception checks supplement, but do not prove, that contract. Vendor-specific arithmetic outside this contract requires a different verifier.

The exact toy tests enumerate all 256 inputs and all 64 states over `GF(4)^3`. They check every emission moment, transformed syndrome mass, compressed interval, and local matrix entry against independent rational enumeration. Tests and source hashes authenticate an implementation; they do not substitute for the arithmetic contract.

Downstream code may thin the dyadic upper physical operators using exact rational coefficients and then round upward. It may also define an artificial birth basis from those chosen upper coefficients. That rebasing requires an exactly normalized rational change of basis and exact intertwining before final upward rounding. Independently rounding normalization weights and treating their sum as one is not justified.

The next step is to verify one G4 fractional bound using these local endpoints. Only after that check should the global verifier cover all outer occupancies and sum their bounds.
