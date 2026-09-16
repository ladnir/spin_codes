# Syndrome bounds beyond single-row occupancy

The transfer in [PROOF_REQUIREMENTS.md](PROOF_REQUIREMENTS.md) needs bounds on
weighted syndrome fibers for every epoch input weight. `general_occupancies.py`
implements those bounds with exact integer combinatorics, then evaluates the
transfer in binary64. `certify_no_constant.py` independently recomputes these
bounds in Arb for the balanced-image map. The full 256-bit certificate passes;
the consolidated status is in [BALANCED_RESULT.md](BALANCED_RESULT.md).

## Exact fiber-size bounds

Fix a weight j, let `C=binom(t,j)`, and write `C_h` for the number of weight-j
inputs with syndrome h under B. There are `m=2^s-1` nonzero syndromes. The
known kernel spectrum gives `C_0=k_j`; put `n=C-k_j`.

Let `K_j(w)` be the binary Krawtchouk polynomial. Fourier inversion and Parseval
give, when A=B^T,

```text
P := sum_h C_h^2 = (C^2 + sum_(w>0) a_w K_j(w)^2) / 2^s.
```

For the m nonzero fiber sizes, their sum is n and their squared sum is
`P-k_j^2`. If their largest value exceeds their mean by x, the other m-1
deviations sum to -x. Cauchy-Schwarz therefore gives

```text
max_(h != 0) C_h <= [n + sqrt((m-1)(m(P-k_j^2)-n^2))] / m.
```

The implementation rounds the square root upward using integer `isqrt`, then
takes the integer floor of the resulting rational upper bound. It also
intersects this bound with n and the Fourier triangle bound

```text
(C + sum_(w>0) a_w |K_j(w)|) / 2^s.
```

Finally, a fixed syndrome fiber is a constant-weight code whose pairwise
distance is at least the B-kernel distance d. Put `r=floor((d-1)/2)` and
`h=min(j,t-j)`. If `h<r`, a fiber has at most one member. Otherwise two distinct
members cannot contain the same `(h-r)`-subset, so its size is at most
`floor(binom(t,h-r)/binom(h,r))`. Complementing inputs preserves pairwise
distance, which justifies using h for this size bound.

Let `R_j` denote the minimum of these integer caps and set
`p_j=R_j/binom(t,j)`. Safe weighted bounds are

```text
c_D,j <= min(d_j, 1-beta_j, p_j z^min_w |w-j|),
c_w,j <= min(f_j(w), min(1-beta_j,a_w p_j) z^|w-j| / a_w).
```

For j=1 and j=2, the implementation instead enumerates all supports, groups
them by their exact syndrome, and records `wt(X+A BX)`. This does not assume
injectivity. The resulting weighted histograms tighten the general caps.

The exact small-map tests compare every cap to enumerated fibers and every
epoch-weight transfer to a rational full-state operator.

## Direct bound for dense Bernoulli inputs

The dense outer counting argument uses an iid Bernoulli(theta) input surrogate.
Fix its incoming state q, with `wt(Aq)=v`, and put `z=exp(-lambda)`.
The emitted moment is

```text
F_v = g0^(t-v) g1^v,
g0 = 1-theta+theta*z,       g1 = theta+(1-theta)*z.
```

Normalize the original input law after weighting it by `z^wt(X+Aq)`.
Under this tilted law, output bits remain independent. Their success
probabilities are `p0=theta*z/g0` outside the support of Aq and
`p1=(1-theta)*z/g1` inside it. This is a change of measure used only to bound
the moment; setup does not sample from the tilted law.

Set `rho_i=|1-2p_i|`. For a dual state a of image weight w, let b be the
overlap of Aa with Aq. Its syndrome Fourier coefficient has absolute value
at most `rho0^(w-b) rho1^b`. Since
`wt(A(q+a))=v+w-2b`, the known image spectrum restricts b. The implementation
maximizes over only those permitted overlaps. It handles `a=q` separately;
when an all-one image exists, it also handles its complementary state
separately. This avoids treating one exceptional difference as available
for every dual state in a shell.

Fourier inversion gives a pointwise syndrome-probability cap `H_v`. Shifting
a syndrome by a fixed state does not change that cap. The lazy branch thus
has weighted point mass at most `epsilon F_v H_v` at any next state. The
fresh branch adds at most `(1-epsilon)F_v/m`. These pointwise bounds define
a second valid shell transfer for the Bernoulli surrogate. The zero-source
row uses the exact weighted kernel mass and sends the remaining mass to D.

The dense calculation may take the smaller terminal moment from this transfer
and the fixed-weight transfer averaged over binomial input counts. It chooses
between complete valid moment bounds, not between incompatible state meanings.
Binary64 discovery selects the witnesses. The final balanced-image certificate
recomputes the selected transfer and its full moment in outward arithmetic.

## Current limitation

The old selected map passes the sparse diagnostic but has unresolved dense
boxes under the cheap mixer. Limited proposal/tilt retuning and periodic full
refresh probes did not close them with the retained box cover. Negative margin
bounds are inconclusive, not counterexamples to the construction.

The same-size no-constant map in [SPECTRUM_REDESIGN.md](SPECTRUM_REDESIGN.md)
removes the exceptional expansion state. Its sparse results are essentially
unchanged. Adaptive splitting and retuning closed all dense boxes; the full
256-bit outward certificate passes both targets. The final map is implemented,
dense-reference tested, and about 4% faster in serial comparisons. See
[BALANCED_RESULT.md](BALANCED_RESULT.md) for replay status and scope.
