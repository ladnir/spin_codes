# Retain the third age through a full-input pair

The encoder is unchanged. This note improves the geometric transfer used
by the distance proof. A full input count means that all512 input positions
are allowed by the paid set P; it does not assert that an actual encoded
message has full support at that update.

The [quota-potential diagnostic](QUARTET_QUOTA_POTENTIAL.md) exposed a large
loss even with exact selected input counts. Its transfer remembers n1,n2
but maximizes the split between n3 and n4 independently at every update.
Here one update is composed with the next before that split is discarded.

## The full-input transition determines the intermediate histogram

Let w be the width and n=(n0,...,n4) the incoming age histogram. Put
alpha=3*lambda-tau and phi(n)=product_i h_i^n_i. The selection weights are

    d=(h1,h2,h3,h4*b^alpha,h4*b^lambda).

At a full-input update every quartet is allowed. The next histogram is

    (0,n0,n1,n2,n3+n4),

and the moment weight is b^(alpha*n3+lambda*n4). In particular, the next
third-age count is exactly the previous n2. It cannot be maximized anew
before the following update.

Let J_(l,k) count allowed-input masks of size l with exactly k positions
in fully allowed quartets. The factor J_(l,k)/binom(w,k) appears in the
existing exchangeable selection formula. Setup is still sampled once;
these counts describe its averaging, not per-message resampling.

## Full first, partial second

Suppose the second input count is l. Define

    E_n(z)=(1+d2*z)^n1 (1+d3*z)^n2 (1+d4*z)^(n3+n4).

If k positions are selected at the second update and v of them come from
its type1 group, the final first two cohorts are (l-k,v). Their exact
weighted coefficient, summed over k,v yielding that pair, is

    [b^(alpha*n3+lambda*n4)/phi(n)]
      * J_(l,k)/binom(w,k) * binom(n0,v)*d1^v
      * h1^(l-k) * [z^(k-v)]E_n(z).

The intermediate potential cancels. For fixed n1,n2 and a=w-n0, the only
dependence on the original n3/n4 split is its initial scalar factor.
Maximizing that factor gives

    max(b^alpha/h3,b^lambda/h4)^(a-n1-n2) / (h1^n1 h2^n2).

The remaining polynomial has three populations n1,n2,a-n1-n2. The same
coupled coefficient bound as before applies over the source-cell vertices.
It retains their common total a before choosing the coefficient marker.
Convolving its coefficients with truncated binomial factors for n0 gives
the outgoing second-cohort buckets. No additional state index is needed.

## Partial first, full second

Suppose the first input count is l. If its selection counts are v0,...,v4,
with k=sum v_i, the final first two cohorts are

    (w-l,l-k+v0).

The product of both rewards and the final potential is represented by

    e=(h2,h3,h4*b^alpha,h4*b^(alpha+lambda),h4*b^(2*lambda)).

The exact term is

    J_(l,k)/binom(w,k) * product_i binom(n_i,v_i)*e_i^v_i
      * h1^(w-l)*h2^(l-k) / phi(n).

The last two factors are outside the product. The first update's selected
type2 count contributes to the second reward through e2. This is the
dependency lost by separately maximizing the intermediate old-age split.
Extract the inactive selection v0 exactly, and bound the remaining coupled
active polynomial uniformly over the incoming cell.

## Checked local effect

Exact small tests compare785 coefficients with direct two-update enumeration.
Another777 comparisons check the uniform FLOAT bounds against exact bucket
sums. Fifteen exact counter rows verify that two ordinary updates and one
blocked update retain the same count labels and terminal vector.

Selected full-width row-sum comparisons use the ordinary potential:

| Source a; n1; n2 | Counts | Separate log2 bound | Coupled log2 bound |
| --- | --- | ---: | ---: |
| 312;64..79;16..31 | 312,512 | 659.84631 | 642.50263 |
| 312;64..79;16..31 | 512,312 | 703.56424 | 660.30587 |
| 408;192..207;96..111 | 312,512 | 667.56904 | 643.38828 |
| 408;192..207;96..111 | 512,312 | 646.61303 | 621.61045 |
| 408;192..207;96..111 | 408,512 | 733.72645 | 706.52719 |
| 408;192..207;96..111 | 512,408 | 688.30119 | 686.75927 |
| 512;0..15;0..15 | 312,512 | 642.95993 | 632.05652 |
| 512;0..15;0..15 | 512,312 | 877.50150 | 766.69060 |

These are local two-update factors, not security margins. The source
helper is `quartet_full_pair.py`; the initial check and screen take10.66s.

## Compose before assessing the improvement

`quartet_pair_regional.py` retains the same selected counts
{0,312,384,392,407,408,424,512} and width16 age buckets as the preceding
restricted diagnostic. It replaces each aligned two-update transition
containing count512 by the minimum of the composed ordinary bound and
the new bound. Other pairs retain their ordinary composition. The initial
row uses its actual zero histogram, and the terminal factor is applied once.

Each block carries the exact sum of its two input counts. Eight blocks
therefore retain the same sixteen-update regional quota, centered at6528.
The normalization still accounts for4096 actual updates. The regional
potential and balanced-coefficient diagnostic are then evaluated as before.

This diagnostic includes every path over its selected count set and every
regional payment allocation at r=2100,c=111572. It does not cover omitted
input counts, other covered/row counts, or directed arithmetic. A result
for this restricted transfer cannot be called an end-to-end certificate.
Runs are serial, remote, one-thread, file-free, with a1200s watchdog and
an8GiB conservative array guard.

The calculation completes in302.90s. After three regional-weight updates,
the FLOAT log2 upper expression is+255164.00196014, versus+263756.41563486
for the preceding separate-update diagnostic. Its balanced lower diagnostic
is+216759.24666338. Both concern the restricted relaxed transfer, not an
actual failure probability. The counter-floor contribution is over275bits
below every evaluated regional output.

Thus retaining this dependency helps, but tightening only pairs containing
count512 is insufficient. Ordinary partial/partial pairs still discard the
intermediate third age. The next [arbitrary-pair calculation](QUARTET_ARBITRARY_PAIR.md)
retains it for those pairs too. Further boundary-weight tuning of the same
paired transfer is not supported by its positive lower diagnostic.
