# Theorem integration for factorized selectors

The factorized certificate combines two kinds of deterministic analytic
witnesses. The combination is valid only after the inner statement is uniform
over the object counted by the outer statement.

## Uniform conditional inner lemma

Write the finite construction seed as `R=(R_out,R_in)`. The variable `R_out`
contains every random value read by the conditioned-row outer argument. In
particular, it includes the code and graph choices used to form the complete
outer word. The variable `R_in` contains the lane bijections, the global packet
permutation, and the recursive state permutations used by the inner transfer.
The construction must sample `R_out` and `R_in` independently.

Fix a message `m` and a value `r_out` in the support of `R_out`. Let
`x=x(m,r_out)` be the complete word before the randomization governed by
`R_in`. Suppose that `x` has packet profile `a`. Fix an inner row `theta`,
consisting of a pole, nine packet fugacities, and one positive 65-entry Collatz
vector. The required inner statement is

```text
Pr_R_in[low output | m, R_out=r_out, x(m,r_out)=x]
  <= 2^(I_theta(a)-H(a)).
```

Here `I_theta` is the unnormalized inner affine bound. The value
`H(a)=log2(M!)-sum_q log2(a_q!)` is the packet-profile orbit normalization.
The displayed inequality must hold for every feasible triple `(m,r_out,x)`.
Its right side depends on `x` only through `a`. Thus, the universal quantifier
over complete words precedes the probability over `R_in`.

To combine the bounds, let `nu` be the nonnegative counting or probability
measure used by the outer argument. The outer row `phi` must prove

```text
integral 1[profile(x(m,r_out))=a] d nu(m,r_out) <= 2^O_phi(a).
```

Multiply the pointwise inner inequality by `d nu(m,r_out)` and integrate over
the profile-`a` fiber. The resulting contribution is at most

```text
2^(I_theta(a)+O_phi(a)-H(a)).
```

This argument uses independence only to retain the fixed law of `R_in` after
conditioning on `R_out`. The analytic rows are deterministic proof choices,
not random variables. The argument fails if `R_in` shares randomness with
`R_out`. It also fails if `I_theta` holds only after averaging over a value
already integrated by the outer bound.

## Fixed marginal selectors

For a certified leaf, fix rational weights `alpha_i` on inner rows and
`beta_j` on outer rows. Both weight vectors are nonnegative and sum to one.
Convexity of upper bounds gives the affine exponent

```text
F(a)=sum_i alpha_i I_i(a)+sum_j beta_j O_j(a)-H(a).
```

The identity

```text
conv{I_i+O_j}=conv{I_i}+conv{O_j}
```

shows that separate marginals lose no pair-mixture power. The selector does
not sample witnesses. It records one fixed proof inequality for the leaf.

A component with a zero variable is eligible only when that profile class is
absent throughout the leaf. In the full-support shard, every positive-weight
inner fugacity and outer packet variable must therefore be strictly positive.

## Vertex and union argument

The function `F` is convex on a cumulative cell. Both component sums are
affine. Also, `-H(a)=sum_q log2 Gamma(a_q+1)-log2 Gamma(M+1)` is convex.
Therefore the maximum of `F` over one exact cumulative cell occurs at a cell
vertex.

The verifier reconstructs the 165 balanced cumulative roots and every prefix
split. It checks 197 terminal cells, exact integer counts, and child-count
conservation. At each reconstructed vertex, it outward-evaluates every
positive-weight component, forms the two rational marginal sums, and subtracts
one outward interval for `H`.

If terminal `l` has exact count `n_l` and verified upper endpoint `U_l`, the
expanded ledger is

```text
sum_l n_l 2^U_l.
```

The collapsed ledger is

```text
binom(262143,8) 2^(max_l U_l).
```

The theorem may use the smaller independently verified upper endpoint. The
expanded endpoint must be at most `-40` for the intended probability claim.

## Integration state

The materializer now freezes every positive-weight component in the current
197-leaf checkpoint. It emits 92 inner rows and 34 outer rows. Each inner row
contains the complete 65-entry Collatz vector obtained by deterministic replay
of its frozen pole, fugacities, and iteration limit. The materializer also
emits a 197-row selector binding with exact rational marginals.

`G8_FACTORIZED_MANIFEST_V2.json` binds both role-separated sources, the selector
binding, the exact geometry source, all five origin sources, the materializer,
and the factorized verifier. Its status is explicitly pending outward replay.
The materialization does not certify an affine endpoint.

The verifier supplies strict parsing for `independent-sum-mixture-v1`, exact
marginal-mass checks, role separation, normalization-once evaluation, and exact
cumulative geometry replay. Static validation accepts both component sources
and all 197 selector rows. The exact geometry still has total count

```text
553017927440720481221689864611271442433 = binom(262143,8).
```

The full outward replay is now complete. The receipt binds 126 hardened
components, 197 leaves, and 8,687 vertex inequalities. The expanded union is

```text
[580922.1202845244781846443146682297122334696524598473041607809543195495710766140174695247875903691387,
 580922.1203803731596459168906997473020538722094109015706139165025333502545517062066965018079483879131].
```

Thus, the frozen factorized selector ledger is outward-valid but does not
prove the `-40` target. The upper endpoint misses that target by
`580962.1203803731596459168906997473020538722094109015706139165025333502545517062066965018079483879131`
bits. The expanded aggregation is tighter than the collapsed aggregation.

The remaining theorem obligation is the uniform conditional inner lemma for
the stated split `R=(R_out,R_in)`. Proving that lemma would justify the
factorized composition, but it would not repair the numerical miss.

The next numerical proof step should inspect the receipt's worst outward leaf,
`h2:010/R`. Its upper vertex value is about `86546.959` bits before count
aggregation. First compare each hardened component interval against the
binary64 discovery affine row at that leaf's active vertices. This comparison
separates interval hardening loss from selector/geometry loss. If hardening
loss is small, retune or refine that leaf using outward component values. If
hardening loss is large, repair the affected inner or outer component bound
before any further geometry wave.

Discovery improvements do not enter this theorem. Only the frozen component
parameters, exact geometry, rational selectors, and independent outward
receipt do.
