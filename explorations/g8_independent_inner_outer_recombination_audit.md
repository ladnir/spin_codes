# Audit of independent inner/outer witness recombination

The reported `257500`-bit improvement on the 181 frozen leaves is a useful
discovery result. It is not yet a certificate result. Recombination is sound
only if the inner bound has the pointwise interface stated below.

## Required bound interface

Fix a profile domain `P` and a profile `a in P`. Let `H(a)` be the logarithm
of the packet-profile orbit normalization. Write the stored affine parts as

```text
I_i(a)=c_i-u_i*a,                 O_j(a)=d_j-v_j*a.
```

The normalization `H` is not part of either affine vector. A valid pair gives
the exponent

```text
F_(i,j)(a)=I_i(a)+O_j(a)-H(a).
```

The outer lemma may count or average outer words with profile `a`. The inner
lemma must then hold for every concrete outer word in that conditioned set.
Equivalently, it may hold after conditioning on all information used by the
outer lemma. This uniform conditional statement permits multiplication of the
two bounds. It does not require independence between the outer and inner
randomness.

The present formulas have the required syntactic separation. The inner
shared-drive transfer uses a pole, packet fugacities, and a Collatz vector.
The conditioned-row outer uses packet variables, BL coefficients, a Cauchy
split, and exact graph and puncture spectra. Their constants and charges add,
and the orbit normalization is subtracted once.

The theorem still needs the semantic uniformity lemma above. Cross-pairing is
invalid if an inner row is only an average over outer words, graph syndromes,
or another random variable also used by the outer row. Two upper bounds on
expectations over shared randomness cannot generally be multiplied. A row is
also not cross-pairable if its Collatz vector, state kernel, or conditioned law
depends on the selected outer witness rather than only on the frozen
construction and the concrete packet profile.

Joint tuning at one anchor profile creates no dependency. Analytic tilts are
deterministic proof parameters, not sampled construction randomness. The term
"conditioned row" in the outer bound likewise creates no coupling by itself.

## Convex-hull identity

Regard each `I_i` and `O_j` as a vector of affine coefficients on `P`. Then

```text
conv{I_i+O_j : i in A, j in B} = conv{I_i : i in A}+conv{O_j : j in B}.
```

For the forward inclusion, fix pair weights `gamma_(i,j)>=0` whose sum is one.
Define their marginals by

```text
alpha_i=sum_j gamma_(i,j),        beta_j=sum_i gamma_(i,j).
```

Then

```text
sum_(i,j) gamma_(i,j)(I_i+O_j)
  =sum_i alpha_i I_i + sum_j beta_j O_j.
```

For the reverse inclusion, fix marginal weights `alpha` and `beta`. The
product coupling `gamma_(i,j)=alpha_i beta_j` has those marginals and gives
the same affine function. Thus pair weights contain no information relevant
to evaluation beyond their two marginals.

Every component is already an upper bound. A convex combination remains an
upper bound because `min_i I_i(a)<=sum_i alpha_i I_i(a)`, and similarly for
the outer family. Hence separate marginal mixtures are proof objects, not a
claim that witnesses are sampled independently.

## Vertex minimax problem

Let `V` be the exact rational vertices of one leaf. Let `A` and `B` contain
only inner and outer rows eligible on the complete leaf. Discovery may solve

```text
minimize    U
subject to  sum_i alpha_i I_i(x) + sum_j beta_j O_j(x) - H(x) <= U
            for every x in V,
            alpha_i>=0,  sum_i alpha_i=1,
            beta_j>=0,   sum_j beta_j=1.
```

This LP has `|A|+|B|` weights instead of `|A||B|` pair weights. Binary64
values may select and rationalize the marginals during discovery. The final
artifact stores fixed reduced rational weights. The independent verifier
replays those weights; it does not trust the discovery objective.

If the leaf has affine dimension `d`, each marginal affine function has a
representation using at most `d+2` source rows. This follows from
Caratheodory's theorem in the at-most `(d+1)`-dimensional space of affine
functions restricted to the leaf. If `m=|V|`, an optimal extreme LP solution
exists with at most `m+1` positive weights across both marginals. The latter
bound follows from the `m` vertex inequalities and two marginal-sum
equalities. These are existence bounds; a serialized selector may use more
components and remain valid.

## Support eligibility

A positive-weight inner row with fugacity zero in class `k` is finite only on
the face `a_k=0`. The same rule applies to any outer row that uses a zero
packet variable. Therefore a component is eligible on an exact-support shard
`S` only if every class in `S` has a strictly positive variable for that
component. Zero variables are permitted only outside `S`.

Eligibility is checked separately for the two marginals. The eligible pair
set is their Cartesian product. A zero mixture weight does not activate its
row. Positive mixture weights cannot average away an ineligible component or
an infinite affine charge.

## Certificate representation

The preferred leaf selector stores two sparse marginals:

```json
{
  "kind": "independent-sum-mixture-v1",
  "inner_marginal": [
    {"witness": {"source_id": "...", "source_sha256": "...",
                  "row": 0, "row_sha256": "..."},
     "weight": "p/q"}
  ],
  "outer_marginal": [
    {"witness": {"source_id": "...", "source_sha256": "...",
                  "row": 0, "row_sha256": "..."},
     "weight": "r/s"}
  ],
  "normalization": "packet-profile-orbit-v1"
}
```

Each marginal has nonnegative reduced rational weights summing exactly to
one. Component references bind the original inner or outer parameter row, not
a producer-computed combined affine row.

If a consumer requires explicit pairs, it may serialize the product coupling.
This uses at most `r*s` pairs for marginal support sizes `r` and `s`. A sparse
rational transportation coupling with the same marginals uses at most
`r+s-1` pairs. Either form evaluates identically, provided every cross-pair is
eligible. The verifier checks row and column sums exactly. It must not accept
pair coefficients without component provenance.

The immutable manifest needs three additions:

1. Separate content-addressed source roles for inner and outer parameter rows.
2. A recombination contract naming the uniform conditional inner lemma and
   the outer profile quantity that it multiplies.
3. A normalization rule declaring the combined evaluator as the sole owner of
   `H(a)`.

The manifest must also bind the outward algorithms and every transitive input
used by each component. Inner inputs include the transfer kernel, split caps,
and code spectra. Outer inputs include the conditioned-row implementation,
split spectra, graph spectrum, and puncture data. Because the current manifest
binds combined atlas rows, this change requires a new manifest version rather
than reinterpretation of an immutable source.

## Independent outward replay

For every leaf, the verifier performs these gates.

1. Resolve each component by source identifier and digest.
2. Reconstruct every inner and outer affine interval from frozen parameters.
3. Check component eligibility against the complete leaf support.
4. Check each marginal's rational weights and exact sum.
5. Reconstruct every exact leaf vertex.
6. At each vertex, form the nonnegative rational weighted interval sum of the
   two marginals.
7. Subtract one outward enclosure of `H(a)`, using its lower endpoint for the
   final upper bound.
8. Compare the resulting upper endpoint with the claimed leaf threshold and
   aggregate the verified leaf bounds through the existing union ledger.

Stored constants, charges, pair sums, diagnostic scores, and the reported
`257500`-bit improvement are nonbinding diagnostics. In particular, the
verifier must not add two producer floats or materialize a binary64 pair bank.

## Conclusion

Independent recombination is mathematically sound under the uniform
conditional inner interface. The convex-hull identity makes separate marginal
mixtures canonical and avoids a quadratic optimization bank. The only
potential obstruction is semantic correlation: if the inner theorem averages
over information retained by the outer theorem, arbitrary cross-pairing is
not justified. Freeze that interface as a hashed lemma dependency before
promoting the diagnostic replay to a certificate run.
