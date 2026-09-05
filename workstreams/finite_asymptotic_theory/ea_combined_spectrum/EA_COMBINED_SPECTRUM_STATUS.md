# EA combined-spectrum status

**Status: paused by user.**  The route does not currently offer enough
improvement over the certified degree-33 construction to justify completing
the rank-two and high-occupation proofs.

## Result

Applying the Reed--Muller lesson to expand--accumulate improves the
occupation-one analysis but does not by itself reduce the required degree at
the existing 512-bit constituent size.

For one sampled right-degree sparse map followed by an accumulator, repeated
in every outer row, direct transfer weighting gives a 44.826-bit
occupation-one diagnostic at degree 23.  The older shell-cap interface needed
degree 33.

The first repeated-row test reverses most of that gain.  When two active
outer rows carry the same local message, the complete degree-23 rank-one
Chernoff sum has only 19.174 bits.  Degree 33 has 42.707 bits.  The dominant
outer image weight is 16.

Larger constituents help without changing the leading raw XOR count.  At
rate one half, there are (L=N/B) constituents, and one degree-(r) stage
costs (L(Br-1)\approx Nr) XORs.  Low-weight partial sums identify the
following candidates:

- ([1024,512]), degree 31;
- ([2048,1024]), degree 29;
- ([4096,2048]), degree 27; and
- ([8192,4096]), degree 25.

Each candidate exceeds 40 bits in the retained equal-message rank-one sum.
None is a certificate: higher outer weights, distinct-message rank two, and
occupations at least three remain open.

## Exact probability space

Fix a rate-half pair ((B,K)) and (L=N/B).  In each setup attempt, sample
independently

\[
S_j\gets\binom{[K]}r\qquad(0\le j<B)
\]

and define

\[
(Ex)_j=\bigoplus_{i\in S_j}x_i,\qquad C=AE,
\]

where (A) is the zero-initialized prefix accumulator.  Setup makes at most
16 independent attempts and accepts the first full-rank (C).  Every one of
the (L) outer rows uses that same accepted map.

Each outer row receives an independent uniform permutation of its (B)
coordinates.  After transposition, each of the (B) regions receives an
independent uniform permutation of its (L) positions.  RandomStepConv-M22
uses an independent uniform binary (23\)-by-(23) linear map at every one of
the (N) positions.  All setup choices are fixed for the sampled code.

## Combined enumerators

For a nonnegative transfer function (f), define

\[
F_f(C)=\sum_{x\ne0}f(\operatorname{wt}(Cx)).
\]

The exact one-word EA law computes (mathbb E F_f(C)).  For nonzero
distinct (x,y), their type is determined by

\[
h_1=\operatorname{wt}(x),\quad
h_2=\operatorname{wt}(y),\quad
h_3=\operatorname{wt}(x+y).
\]

The sparse-map pair symbol has probabilities

\[
p_{a,b}=\frac14\left(1+(-1)^a\beta_r(h_1)
 +(-1)^b\beta_r(h_2)+(-1)^{a+b}\beta_r(h_3)\right),
\]

where (eta_r(h)=K_r^{(K)}(h)/\binom Kr).  A four-state pair accumulator
then computes

\[
\mathbb E[F_f(C)^2]
=\sum_{x,y\ne0}
  \mathbb E\!left[f(\operatorname{wt}(Cx))
                    f(\operatorname{wt}(Cy))\right].
\tag{1}
\]

`verify_combined_functional_small.py` checks (1) exactly against all 1,024
sparse maps at ((K,B,r)=(4,5,3)), using a non-exponential positive weight
vector.  Both rational moments agree.

For the actual Q2 transfer, equal local messages need only the one-word
spectrum law because the two row permutations are independent conditional on
the common outer weight.  `evaluate_ea_combined_q2_rank_one.py` evaluates
this sector.  Distinct messages require the full pair law above.

## Interpretation

The degree-33 requirement was not solely an artifact of demanding standalone
outer distance.  The simultaneous zero caps below weight 42 were stronger
than needed for Q1, but rare low-weight images become important when the same
local message appears in two active rows.  The combined spectrum sees this
effect directly.

The larger-block trend is nevertheless useful.  It trades constituent size
for degree while retaining sparse linear-time encoding.  It also satisfies
the user's preference for one fixed repeated constituent.  At the finite
(k=2^{20}) target, block sizes through 8192 remain sublinear, although their
cache behavior and optimized XOR schedules have not been measured.

## Proof obligations

1. Optimize and outward-round the rank-one Q2 transfer.  The current common
   tilt is (u=-7.5); a weight-dependent finite witness set may improve the
   margins.
2. Evaluate the distinct-message rank-two combined functional for the best
   larger-block candidate.  A direct enumeration of all message-pair types is
   not acceptable at large (K); the proof needs a positive low-rank or
   character majorant.
3. Construct a finite family of nonnegative one-row functionals whose
   products dominate every occupation (Q\ge3).  Certify that family by its
   first and second moments, then reuse the existing dense-occupation inner
   transfer.
4. Add the omitted outer-weight tail to each truncated candidate receipt.
5. After one candidate closes mathematically, sample one full-rank map,
   optimize its transposed XOR circuit, and benchmark it alone against the
   frozen degree-33 and BCH kernels.
6. Transfer the final outer interface from RandomStepConv-M22 to RM2Sub-S19.

## Recommendation

Keep this route paused.  If it is resumed, do not optimize or benchmark a
larger constituent first.  Improve the rank-one Q2 witness and compute the
rank-two combined bound for the ([2048,1024]), degree-29 candidate.  It is
the smallest candidate with about 4.1 retained bits of slack and is less
likely than the larger blocks to lose its theoretical XOR advantage to code
size and live-set pressure.
