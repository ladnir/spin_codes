# Full certificate plan

## Target statements

Fix message length \(n=2^{20}\), output length \(N=2^{21}\), and distance
threshold \(d=188743\). Fix a binary linear \([256,128]\) outer code with
weight spectrum \((A_w)_{w=0}^{256}\). Assume that \(A_0=A_{256}=1\).

Setup samples the coordinate permutations, region permutations, and nonzero
checkpoint multipliers specified in `CONSTRUCTION.md`. Let \(X\) denote the
number of nonzero messages whose encoded output has weight at most \(d\).
The original 9% certificate target is

\[
 \mathbb E[X] < 2^{-40},
\]

where the expectation is over setup. Markov's inequality then gives

\[
 \Pr[X>0] < 2^{-40}.
\]

The current modeled spectrum is not the spectrum of an identified code. A
certificate using that model remains conditional. An unconditional
construction claim requires an exact spectrum or a proved upper envelope for
an explicit outer code.

The middle-density calculation shows that this first-moment route does not
close at 9%. The active fallback target sets

\[
 d=\lfloor 0.06N\rfloor=125829.
\]

All statements below that claim a 6% margin use this threshold.

## Completed components

The exact epoch transfer closes on the state classes zero and uniform
nonzero. The coefficient formula covers every epoch occupancy from zero
through 1024.

The one-active modeled-spectrum calculation gives 81.2667 bits of margin.
The complete two-active modeled-spectrum calculation sums all 4278 unordered
weight pairs and every support intersection. It gives 159.8280 bits of
margin. The pair \((38,38)\) dominates.

The regular-spectrum envelope removes the outer-weight tuple. For every
support \(S\subsetneq[256]\) other than the empty support, define

\[
 \mu(S):=\frac{A_{|S|}}{{256\choose |S|}}.
\]

For the modeled spectrum and every regular support,

\[
 \mu(S)
 \le
 2\frac{2^{128}-1}{2^{256}}.
\]

Thus each regular active block is dominated by a uniform random 256-bit
outer word at a factor-two cost. Configurations with \(a\) regular active
blocks reduce to one occupation parameter \(a\).

The floating-point envelope sum for \(1\le a\le64\) has 55.9008 bits of
margin at 9%, so the same bound applies at 6%. The log-domain region
recurrence covers every remaining regular occupation at 6%. The exact middle
range through occupation 5466 has 834.4399 bits of aggregate margin. An
analytic Hamming-ball bound covers occupations 5467 through 8192 with
107.1058 bits of aggregate margin.

The exceptional layer with exactly one all-one outer word also closes. The
identity

\[
 F_{a,1}=2F_{a+1,0}-F_{a,0}
\]

gives 1161.5591 bits of aggregate margin from 65 through 5487 regular
blocks. The earlier exact mixed calculation covers total occupation at most
64. Starting at 5488 regular blocks, the analytic dense bound sums every
possible all-one count and retains 91.3725 bits.

## Remaining components

The conditional 6% certificate needs two additional components.

1. Bound configurations with at least two all-one outer words and fewer than
   5488 regular outer words. The density envelope excludes the unique support
   \([256]\).
2. Recompute every retained term with outward-rounded arithmetic and fixed
   Chernoff tilts. The checker must sum all certified upper bounds without
   optimization.

The implemented bulk method computes the regular region matrix \(\overline
R_a(z)\) for every \(a\). A region first selects a uniform \(a\)-subset of
candidate positions, then assigns independent fair bits to those positions.
If \(C_c(z)\) denotes the epoch transfer with \(c\) candidate positions,
then

\[
 C_c(z)=2^{-c}\sum_{r=0}^{c}{c\choose r}M_r(z).
\]

The exact regular region matrix is

\[
 \overline R_a(z)=
 \frac{[u^a]
 \left(\sum_{c=0}^{1024}{1024\choose c}C_c(z)u^c\right)^8}
 {{8192\choose a}}.
\]

The inner moment for occupation \(a\) is

\[
 e_0^T\overline R_a(z)^{256}\mathbf1.
\]

The checker implements this recurrence in the log semiring. It does not use
the exploratory FFT calculation, whose roundoff floor was too large for the
required coefficients.

## Next implementation goal

Extend the stable finite-difference relation

\[
 F_{a,b+1}=2F_{a+1,b}-F_{a,b}
\]

to the first several all-one layers. Compare each layer with the direct mixed
enumerator on the existing low-occupation overlap. If the layer margins grow
after charging their placement multiplicity, prove a monotone tail bound in
\(b\). Otherwise retain exact layers until a separate forced-impulse bound
takes over.
