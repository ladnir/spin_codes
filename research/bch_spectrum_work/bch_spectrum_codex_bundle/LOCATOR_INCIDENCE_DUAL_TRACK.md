# A dual track for the weight-38 shell

Application correction, 2026-09-04: the weight-38 cap below is an allocation
from a joint low-shell budget. It does not alone close the 40-bit Q1 screen
with the current other-shell LP caps. The directed combined upper bound then
has 36.8688 bits. See `generated/endpoint_application_budget_audit.json`.
All locator and incidence statements below concern the weight-38 obligation;
later language about closing the application must be read with this correction.

## The narrow question

Let \(D:=\mathbb F_{256}^{*}\). For each \(18\)-set \(B\subset D\), let
\(g_B\in\mathbb F_{256}[u]\) be the unique polynomial such that

\[
 \deg g_B\leq18,\qquad g_B(0)=1,\qquad
 g_B(u)=u^{146}\quad\text{for every }u\in B.
\]

Define

\[
 N_B:=\#\{u\in D:g_B(u)=u^{146}\},
 \qquad E_B:=N_B-18.
\]

The polynomial \(z^{37}+g_B(z^2)\) is nonzero and has degree 37. The map
\(z\mapsto z^2\) permutes \(D\). Hence \(18\leq N_B\leq37\), so
\(0\leq E_B\leq19\).

Let \(B\) be uniform over the \(18\)-subsets of \(D\). Double counting gives
the exact identity

\[
 R_6:=
 \frac{\mathbb E_B\binom{E_B}{6}}
      {\binom{237}{6}/256^6}
 =
 \frac{I_{24}}{\binom{255}{24}/256^6},
\]

where \(I_{24}\) is the number of \(24\)-sets on which some admissible
polynomial agrees with \(u^{146}\).

Every normalized weight-37 locator contributes \(\binom{37}{24}\) such
sets. Therefore

\[
 L_{37}\leq
 R_6\frac{\binom{255}{24}}
          {256^6\binom{37}{24}}.
\]

An outward-rounded transfer calculation, followed by exact rational
aggregation, gives the allocated integer shell cap

\[
 A_{38}(C)\leq3{,}827{,}351{,}840{,}403.
\]

A sufficient incidence condition for this shell allocation is

\[
 R_6<6.098404530135757.
\]

The displayed decimal is only a readable lower endpoint; the receipt stores
the exact rational endpoint. In particular, the following standalone lemma
suffices for the weight-38 allocation.

**Target lemma.** For the experiment above, \(R_6\leq6\).

The target lemma gives

\[
 \log_2 L_{37}\leq34.0698<34.0933.
\]

This margin is only \(0.02347\) bits. A proof with constant \(5\) would leave
\(0.28650\) bits. These margins already include directed rounding in the
shell allocation. These are not margins for the full Q1 bound. The receipt is
generated/random_inner_threshold_outward.json.

## Heuristic track

The heuristic claim will concern \(R_6\), not the full weight distribution and
not the probability that a candidate splits completely. This distinction is
necessary: the exact \(q=128\) endpoint is strongly enriched even though its
sixth incidence moment is not.

The current size ladder is:

| Field | Locator family | Statistic | Ratio to independent-binomial reference |
|---|---:|---:|---:|
| \(q=8\) | \(t=1,d=3\) | exact order 2 | \(1.829\) |
| \(q=32\) | \(t=3,d=7\) | exact order 4 | \(1.994\) |
| \(q=128\) | \(t=10,d=21\) | sampled order 6 | \(0.784\), normal 95% upper \(0.840\) |
| \(q=256\) | \(t=18,d=37\) | sampled order 6 | \(1.032\), normal 95% upper \(1.136\) |

The first two rows enumerate every base interpolant. The last two rows use
five million independent samples each. The \(q=128\) normalized locator list
has exact size 330, which is \(96.94\) times its random full-splitting
reference. Thus the endpoint itself is not random-like. The lower incidence
moments remain close to, or below, their binomial references. The exact
endpoint contributes only 0.362 to the \(q=128\) sixth-order ratio of 0.784.
Endpoint enrichment therefore does not force inflation past the available
constant-factor budget.

A conservative heuristic assumption is

\[
 \mathsf{H}_{\mathrm{inc}}:\qquad R_6\leq2.
\]

This assumption implies

\[
 L_{37}\leq 6.011\times10^9,
 \qquad
 A_{38}(C)\leq1.256\times10^{12}.
\]

It leaves \(1.608\) bits below the application threshold. The assumption is
falsifiable and refers to the exact functional that enters the application.

The \(q=256\) sampler also stratifies \(B\) by \(|B\cap\mathbb F_{16}^{*}|\).
The observed contribution is not concentrated on bases rich in the largest
proper subfield. This test addresses the most obvious hidden structured
component, but it does not exclude every exceptional component.

Two new strata target the character-sum concerns more directly. In the same
five-million-sample run, interpolants of degree 17 contribute only 0.00463 to
the total sixth-moment ratio; degree-18 interpolants contribute 1.02708.
Conditioning on whether the linear coefficient belongs to the 51-element
image \(\{a^{130}:a\neq0\}\) changes the estimated conditional ratio from
1.024 outside the image to 1.063 inside it. The corresponding normal
95-percent upper endpoint inside the image is 1.295. These diagnostics show
no concentration near the affine Walsh ridge, but they remain sampling
evidence.

## Proof track

The proof track aims first for \(R_6\leq5\), with \(R_6\leq6\) as the minimum
successful result.

Two facts are already exact.

1. If the interpolant satisfies \(g_B'\neq0\), the six incidence equations
   have Jacobian rank six. The nondegenerate incidence locus is smooth of the
   expected codimension.
2. If \(g_B'=0\), then \(g_B=h^2\) with \(\deg h\leq9\). Nine agreement
   points determine \(h\), which gives

   \[
   L_{37}^{\mathrm{deg}}
   \leq
   \left\lfloor
   \frac{\binom{255}{9}}{\binom{37}{9}}
   \right\rfloor
   =87{,}550{,}900.
   \]

The remaining work is a point bound on the smooth locus. Three approaches can
be pursued without changing the target lemma.

1. Expand the six agreement indicators with additive characters. Use the
   binary expansion \(146=128+16+2\) to simplify the resulting sums.
2. Project the smooth incidence variety onto a smaller set of symmetric
   coordinates. Bound the fibers after isolating subfield and
   derivative-zero strata.
3. Construct a computer-assisted orbit count or upper bound. The certificate
   should expose exact field operations and orbit multiplicities, so an
   independent checker can verify the result without rerunning the search.

Generic smooth-variety estimates are unlikely to suffice at \(q=256\), because
their degree-dependent error constants are too large. The proof must exploit
the monomial exponent, the symmetric-set formulation, or both.

The first approach has an exact target. Let
\(\Phi(x_1,\ldots,x_{24})\in\mathbb F_{256}^{6}\) contain the constant and
five high-coefficient defects of the interpolating remainder. Fix a nontrivial
additive character \(\psi\) of \(\mathbb F_{256}\). For
\(a\in\mathbb F_{256}^{6}\), define

\[
 S_a:=
 \sum_{\substack{x_1,\ldots,x_{24}\in D\\
                  x_i\neq x_j\text{ for }i\neq j}}
 \psi\bigl(a\mathbin{\cdot}\Phi(x_1,\ldots,x_{24})\bigr).
\]

Character orthogonality gives

\[
 24!I_{24}
 =256^{-6}\sum_{a\in\mathbb F_{256}^{6}}S_a,
 \qquad
 S_0=(255)_{24}.
\]

Consequently, the target lemma follows from the signed aggregate bound

\[
 \sum_{a\neq0}S_a\leq5(255)_{24}.
\]

The exponent \(146=128+16+2\) has binary weight three. Thus
\(u\mapsto u^{146}\) is cubic when \(\mathbb F_{256}\) is represented as an
eight-dimensional vector space over \(\mathbb F_2\). This low algebraic degree
is invisible to a generic degree-146 estimate and is the main reason to try a
character-sum argument. A termwise absolute bound on all \(S_a\) may still be
too expensive; the required inequality concerns their signed aggregate.

There is a second exact form of the character expansion. Write an admissible
polynomial as

\[
 g(u)=1+\sum_{j=1}^{18}c_j u^j.
\]

For an ordered tuple \(X=(x_1,\ldots,x_{24})\) of distinct nonzero field
elements, define

\[
 K_X:=\left\{\lambda\in\mathbb F_{256}^{24}:
       \sum_i\lambda_i x_i^j=0\text{ for }1\leq j\leq18\right\}.
\]

Summing first over the coefficients \(c_j\) gives

\[
 24!I_{24}
 =256^{-6}\sum_X\sum_{\lambda\in K_X}
   \psi\!\left(\sum_i\lambda_i(x_i^{146}+1)\right).
\]

The Vandermonde matrix has rank 18, so \(K_X\) has dimension six. If
\(P_X(z)=\prod_i(z-x_i)\), every element of \(K_X\) has the unique form

\[
 \lambda_i=\frac{h(x_i)}{x_iP_X'(x_i)},
 \qquad \deg h\leq5.
\]

This kernel form is a useful proof reduction: the six dual variables are the
six coefficients of \(h\). It also exposes a limitation of a direct Walsh
argument. The denominator \(P_X'(x_i)\) couples all 24 points, so the inner
sum does not factor into ordinary one-variable Walsh coefficients.

Two further exact reductions delimit the useful proof routes.

First, let \(B_n\) be the number of admissible polynomials with exactly \(n\)
agreements. Any set of at most 18 nonzero evaluation points imposes
independent conditions on the 18 free coefficients. Therefore

\[
 \sum_{n=0}^{37}B_n\binom ns
 =\binom{255}{s}256^{18-s},
 \qquad 0\leq s\leq18.
\]

An exact rational linear program maximizes
\(\sum_nB_n\binom n{24}\) subject to these identities, \(B_n\geq0\), and
\(n\leq37\). Its primal and dual optima agree exactly. The resulting bound is

\[
 R_6\leq33{,}072{,}431.234.
\]

The extremizer is supported on \(0,1,\ldots,17,37\); it places all permitted
post-18 mass at the endpoint. Thus the universal Reed--Solomon moments alone
are short by about 22.39 bits relative to the factor-6 target. A successful
proof must use a monomial-specific constraint that excludes this endpoint
extremizer. The receipt is generated/rs_coset_moment_lp.json.

Second, the monomial-specific equations have a compact symmetric form. For a
24-set \(T\), let \(e_j(T)\) and \(h_j(T)\) denote its elementary and complete
homogeneous symmetric functions. Lagrange interpolation gives

\[
 [u^j](u^{146}\bmod P_T)
 =\sum_{r=0}^{23-j}e_r(T)h_{146-j-r}(T),
 \qquad 19\leq j\leq23.
\]

The five high-coefficient equations are triangular. Together with the
constant-coefficient equation, admissibility is exactly

\[
 e_{24}(T)h_{122}(T)=1,\qquad
 h_{123}(T)=h_{124}(T)=h_{125}(T)=h_{126}(T)=h_{127}(T)=0.
\]

In characteristic two, the generating functions satisfy
\(H(t)=E(t)H(t)^2\). Hence

\[
 h_n=\sum_{\substack{0\leq i\leq24\\i\equiv n\pmod2}}
 e_i h_{(n-i)/2}^{\,2}.
\]

This identity halves the indices in all six equations and retains the
Frobenius structure of exponent 146. It is a more promising basis for a
computer-assisted count than the degree-146 remainder equations. The
field-convention checks are in generated/symmetric_incidence_identity.json.

An exhaustive Fourier calculation gives a more precise building block. For
\(a\neq0\) and \(b\in\mathbb F_{256}\), define

\[
 W(a,b):=
 \sum_{x\in\mathbb F_{256}}
 (-1)^{\operatorname{Tr}(a x^{146}+b x)}.
\]

Every coefficient belongs to

\[
 \{-32,-16,0,16,32,64\}.
\]

Of the 65,280 coefficients, 26,775 vanish. The value 64 occurs exactly for
the 255 pairs \((a,b)=(a,a^{131})\). Every other coefficient has absolute
value at most 32. The vectorial monomial \(x\mapsto x^{146}\) also has
differential uniformity 6.

These facts are exact, but they do not yet bound \(I_{24}\). The ordinary
Walsh spectrum applies directly only to affine perturbations

\[
 u^{146}+cu+1.
\]

For this slice, the exceptional condition becomes \(c=a^{130}\). Its signed
numerator contribution is

\[
 64\sum_{a:a^{130}=c}(-1)^{\operatorname{Tr}(a)}.
\]

Summing over all \(c\) counts each nonzero \(a\) once and gives \(-64\).
Thus the exceptional ridge cancels in aggregate on the affine slice; it does
not accumulate with one sign. Exact enumeration gives a stronger locator
consequence:

| perturbation degree | largest number of roots |
|---:|---:|
| at most 1 | 4 |
| at most 2 | 7 |
| at most 3 | 10 |

Every one of these values is below 18. Hence all perturbations of degree at
most three contribute zero to both the 18-base experiment and \(I_{24}\).
The receipts are generated/monomial_146_spectrum.json and
generated/low_degree_146_slices.json.

The known locator sample lies at the opposite end of the degree range. The
138 certified affine orbits generate 5,244 normalized locators: 5,206 have
degree 18, 31 have degree 17, and seven have degree 16. This is finite evidence
about discovered locators, not an exclusion theorem for unknown components.

The Frobenius-fixed coefficient family can be exhausted completely. Restrict
all 18 free coefficients of \(g\) to \(\mathbb F_2\). Among the resulting
\(2^{18}\) polynomials, exactly 82 have at least 24 agreements and exactly 10
have 37 agreements. No polynomial has between 29 and 36 agreements. The
complete contribution of this family is

\[
 R_6^{(\mathbb F_2\text{-coeff})}
 =3.3276472005\cdot10^{-9}.
\]

The set-side view is deliberately non-random. Frobenius on
\(\mathbb F_{256}^{*}\) has orbit-size distribution
\(1^1 2^1 4^3 8^{30}\). Exactly 5,365 unions of these orbits have size 24,
and 112 are admissible. The independent-rank expectation is only
\(5365/256^6\), so this stratum is enormously enriched. Nevertheless, those
112 sets contribute only \(1.0462\cdot10^{-17}\) to the global \(R_6\)
normalization. This is direct evidence that exceptional algebraic components
exist but need not threaten the aggregate bound.

The ten binary endpoint locators occupy five affine orbits. Four were absent
from the radius-five search and have now been certified directly. Each has
trivial affine stabilizer and orbit size 65,280. Consequently, the rigorous
search floors improve to

\[
 A_{38}(P)\geq7{,}768{,}320,\qquad
 A_{38}(C)\geq944{,}384.
\]

The receipts are generated/binary_coefficient_146_family.json and
generated/frobenius_invariant_incidence.json. These are lower-bound and
exceptional-stratum results; they do not upper-bound the remaining
coefficient families.

The next subfield stratum is also exact. Under \(x\mapsto x^4\), the nonzero
field elements form 3 singleton, 6 pair, and 60 four-element orbits. Exhaustive
enumeration of all 267,516,561 invariant 24-sets finds 67,900 admissible sets
and 40,712 distinct interpolants with coefficients in \(\mathbb F_4\). Thus the
full \(\mathbb F_4\)-coefficient contribution is

\[
 R_6^{(\mathbb F_4\text{-coeff})}
 =1.1338277402\cdot10^{-8}.
\]

The six-equation hit count is only 1.03963 times its independent
\(\mathbb F_4\) reference. Exactly 34 interpolants have 37 agreements, and
they supply 99.786 percent of the subfield contribution. These endpoints
occupy 15 affine orbits. Ten are new, all with trivial stabilizer, so the
rigorous floors become

\[
 A_{38}(P)\geq8{,}421{,}120,\qquad
 A_{38}(C)\geq1{,}023{,}744.
\]

The receipt is generated/frobenius4_invariant_incidence.json. This exact
subfield result still does not upper-bound coefficients outside
\(\mathbb F_4\).

The \(\mathbb F_{16}\) stratum is too large for the same exhaustive method:
there are 348,770,148,711,057,905 invariant 24-sets. Instead, sample an
invariant 18-set \(B\) and interpolate the unique \(g_B\). Invariance and
uniqueness force every coefficient of \(g_B\) into \(\mathbb F_{16}\).
For each invariant agreeing extension \(T\supset B\), weight the sample by
\(1/m(T)\), where \(m(T)\) counts the invariant 18-subsets of \(T\). This
weight makes the estimator unbiased for the number of admissible invariant
24-sets.

Ten million samples give ratio 1.00624 with standard error 0.00412. The normal
95-percent upper endpoint is 1.01432. These values are statistical evidence,
not a point-count certificate. The same run yields ten exact endpoint
polynomials. They occupy nine previously unknown affine orbits, all with
trivial stabilizer. Their affine closures contain 342 normalized locators,
including 64 with coefficients in \(\mathbb F_{16}\). The exact witness
closures improve the rigorous floors to

\[
 A_{38}(P)\geq9{,}008{,}640,\qquad
 A_{38}(C)\geq1{,}095{,}168.
\]

The receipt is generated/frobenius16_locator_extension_monte_carlo.json. Its
incidence estimate is heuristic; its nine orbit records and shell
consequences are exact.

This calculation also closes off an overly optimistic proof route. One cannot
separate the 255 ordinary Walsh ridge modes inside the full incidence sum
without first reducing the higher-degree phase to an affine one. The next
character-sum step must instead exploit the six-dimensional kernel form,
derive generalized spectra for higher-degree perturbations, or average over
the polynomial coefficients before applying absolute values.

## Evidence discipline

The exact reduction, the endpoint degree bound, the Jacobian theorem, and the
derivative-zero bound are proof claims. The \(q=8\) and \(q=32\) ladder rows
are exhaustive computations. The \(q=128\) and \(q=256\) rows are Monte
Carlo evidence. Neither their normal confidence endpoints nor
\(\mathsf H_{\mathrm{inc}}\) is a proof.

The machine-readable ladder is generated/locator_size_ladder.json. The
sampling source is code/sample_locator_extensions.cpp.
