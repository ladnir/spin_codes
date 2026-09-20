# A cheaper mixer with an exact state law

The candidate changes only the state mixer of the quarter-rate SPIN inner.
It retains the selected RM2Sub maps with t=128, s=19 and the fixed [128,32,32]
outer. We have proved the mixer law and a weighted-state envelope, and evaluated
the **single-active-outer-row contribution** with outward arithmetic for this
original map. Its dense cover remains unresolved. The later balanced-image map
has a full outward certificate and measured speedup; see
[BALANCED_RESULT.md](BALANCED_RESULT.md). The generic transfer argument below
applies to both maps using their respective spectra and cancellation bounds.

## Interface and probability space

For one binary message, an epoch takes `X` in F_2^t and incoming state `q` in
F_2^s. It emits `Y=X+Aq` and updates `q'=Mq+BX`. Initially `q=0`; the state
continues across regions, and no final state weight is charged. The selected
maps satisfy `B=A^T`, `BA=0`, rank(A)=19, minimum nonzero A-image weight 48,
and minimum B-kernel weight 6.

Setup samples the mixers independently across epochs and independently of
the permutations. All messages share that setup. For a fixed message, once
epoch input counts are conditioned on, each epoch support is uniform among
supports of its count, independently of the incoming state. The region
coefficient calculation performs this conditioning before averaging counts.
The first-moment union over messages does not assume independent encodings
of different messages.

The current uniform nonzero field multiplier sends every fixed nonzero state
to a uniform nonzero state. Its implementation computes a general binary state
matrix in each transposed epoch. We ask whether a much cheaper update can
retain enough of that marginal property for the distance proof.

## The rank-one mixer

Sample a nonzero vector `u` uniformly from F_2^s. Conditional on `u`, sample
`v` uniformly from `u^perp`, **including zero**. Define `T=I+uv^T`.
Orthogonality gives `(uv^T)^2=0`, hence `T^2=I`: every sampled map is invertible.
The case `v=0` is the identity. Otherwise this is a transvection.

For context, random transvections are a classical matrix-walk construction;
see [Hildebrand, Generating Random Elements in SL_n(F_q) by Random Transvections](https://deepblue.lib.umich.edu/items/7cc1acac-8e85-489a-97b1-5061cfa1a5b1).
We do not require mixing to a uniform random matrix. The elementary single-state
law below is derived directly for our sampling convention, which includes the
identity, rather than imported from a group-mixing theorem.

**Lemma.** Put `m=2^s-1`. For every fixed nonzero state `q`,

```text
Pr[Tq=y] = (1/2) 1[y=q] + 1/(2m),     y != 0.
```

**Proof.** If `y` differs from both zero and `q`, the only possible choice is
`u=q+y`. Then `u` and `q` are linearly independent over F_2. Exactly half the
vectors in `u^perp` satisfy `v^Tq=1`, giving probability `1/(2m)`. Invertibility
excludes output zero. The remaining probability at `q` is `(m+1)/(2m)`.

For a product of `r` independent rounds, the nonzero-state kernel is
`epsilon I + (1-epsilon) J`, where `epsilon=2^-r` and `J` is uniform refresh.
This follows from `J^2=J`. Thus

```text
Mq has law epsilon delta_q + (1-epsilon) Uniform(F_2^s \ {0}).
```

This identity describes a marginal law, not an implementation that chooses a
different fresh state for each message. Setup samples linear matrices once.

A uniform sampler for `v` starts with a uniform s-bit word and flips one fixed
set-bit position of `u` iff its dot product with `u` is one. Each output in
`u^perp` has exactly two preimages.

## What the old proof cannot simply reuse

When `BX` is nonzero, fresh refresh cancels the state with probability `1/m`.
The lazy branch can instead cancel with certainty conditional on `q=BX`.
Replacing a full multiplier by a transvection while retaining the old `1/m`
termination bound is invalid. Permutation spreading helps average the input
syndrome but does not remove this event.

There is a second issue: previous emitted weights bias the entering state law.
An unconditional mixing statement cannot be applied as though this weighted
law were uniform. Our envelope tracks weighted measures explicitly.

## Weight-class envelope

For each nonzero A-image weight `w`, let `S_w={q:wt(Aq)=w}` and `a_w=|S_w|`.
For the retained map the weights are 48, 56, 64, 72, 80 and 128, with counts
5,040; 110,848; 292,510; 110,848; 5,040; and 1.

Fix `0<z<=1`. The measure on the current state weights each setup by
`z^(weight emitted so far)`. Represent an upper envelope by coordinates
`(Z,D,(C_w)_w)`: mass `Z` is at zero; mass at most `D` has an arbitrary law
on nonzero states; and the remaining measure has point mass at most
`C_w/a_w` at each state in `S_w`. The coordinates are unnormalized upper
bounds. They do **not** assert uniform conditional distributions in the shells.

For a uniform weight-j input, define

```text
f_j(w) = E_X z^wt(X+Aq),     q in S_w,
d_j    = max_w f_j(w),
beta_j = Pr[BX=0],
c_D,j  >= max_(q != 0) E_X [1[BX=q] z^wt(X+Aq)],
c_w,j  >= (1/a_w) E_X [1[BX in S_w] z^wt(X+A BX)].
```

The value `f_j(w)` depends only on `w`, by the hypergeometric overlap formula.
Write `eta=1-epsilon`. The following nonnegative, row-to-column transfer is
valid for j>0; unspecified entries are zero:

| Source | Destination | Weight |
|---|---|---|
| Z | Z | `beta_j z^j` |
| Z | D | `(1-beta_j) z^j` |
| D | Z | `epsilon c_D,j + eta min(d_j,1-beta_j)/m` |
| D | D | `epsilon d_j` |
| D | C_w | `eta (a_w/m) d_j` |
| C_v | Z | `epsilon c_v,j + eta min(f_j(v),1-beta_j)/m` |
| C_v | D | `epsilon f_j(v)` |
| C_v | C_w | `eta (a_w/m) f_j(v)` |

For j=0 there is no cancellation. Replace the lazy `C_v -> D` entry by
`C_v -> C_v = epsilon z^v`, adding it to the fresh entry. Also `Z -> Z=1`.
The same retention is valid whenever `beta_j=1`: every supported input has
zero syndrome, so the lazy branch retains the state and its shell. In
particular the all-one input must not be discarded as an exceptional endpoint.

**Justification.** In the lazy branch, the state before feedback is unchanged.
The displayed cancellation definitions bound its zero mass. Nonzero surviving
mass is at most the full emitted moment; assigning this moment to D is safe.
On an empty epoch a shell is unchanged, so its pointwise envelope is multiplied
by `z^v` without forgetting the shell. In the fresh branch, conditional on each
old nonzero state and input, any specified output state has probability at most
`1/m`, independently of the emitted weight. Summing the old weighted measure
gives the displayed shell caps. A transition to zero additionally requires
nonzero syndrome, giving the intersections in the zero column. The zero and
live bounds may overlap in mass; overcounting is conservative. The construction
is closed under addition of incoming envelope components, so it iterates.

For Q=1, only j=0,1 occur. All 128 columns `b_p` of B are distinct and nonzero.
Put `h_p=wt(e_p+A b_p)`. Then the lazy cancellation bounds are exact:

```text
c_D,1 = max_p z^h_p / t,
c_w,1 = sum_(p: b_p in S_w) z^h_p / (t a_w).
```

Retaining shells on empty epochs is the crucial refinement. The original
three-class screen forgot that information and lost about sixteen margin
bits at one round. Maximizing separate shell-to-shell caps on active epochs
also introduces large slack; simply sending those lazy survivors to D is
tighter in this sparse regime.

## Q=1 result

There are L=32768 outer rows and 256 epochs per region. Let `T_0,T_1` be the
eight-state transfers. The empty region is `R_0=T_0^256`; the single-input
region is the degree-one coefficient of `(T_0+x T_1)^256`, divided by 256.
For an outer word of weight `w`, average its support by taking the degree-w
coefficient of `(R_0+x R_1)^128`, divided by `binom(128,w)`. Start at Z and
sum all final coordinates, without a state reset between regions. Multiply
by `L A_outer,w`, sum outer weights, and apply `z^(-bad_weight)`.

For one transvection round and K=2^20, the resulting outward Q1 union bounds are:

| Distance threshold | Bad-weight cutoff | Certified Q1 margin |
|---|---:|---:|
| 16.5% | 692,060 | 41.048247 bits |
| 19% | 796,917 | 30.034070 bits |

[TRANSVECTION_Q1_CERTIFICATE.json](TRANSVECTION_Q1_CERTIFICATE.json) stores exact
upward dyadics produced at 256-bit precision. A 512-bit replay is enclosed by
those bounds. The accompanying tests exhaustively check small-state mixer
laws and compare three-epoch envelopes against rational full-state operators.
These checks supplement the argument above; they do not prove higher occupancy
coverage by testing small instances.

**This is not a full SPIN certificate.** Subsequent binary64 calculations cover
Q=2..64 at 16.5% and Q=2..128 at 19%; their Q=2 margins are 55.27 and 42.43 bits.
The dense cover remains unresolved. At 19% the Q1 bound consumes about 97.7%
of the 30-bit budget. See [FIBER_BOUNDS.md](FIBER_BOUNDS.md) for the general
transfer and [SPECTRUM_REDESIGN.md](SPECTRUM_REDESIGN.md) for the next map.

## Remaining proof and engineering work

1. General `c_D,j` and `c_w,j` bounds are now implemented and tested using
   constant-weight packing, Fourier/Parseval, and exact small-input histograms.
   The all-occupancy calculation still needs an outward verifier.
2. Reuse the region coefficient and outer counting arguments with the new
   transfer, but recompute all witnesses. Cover sparse, intermediate and dense
   occupations, including full input and unflushed terminal cases. Old numeric
   certificates do not transfer automatically.
3. Preserve zero-state obstructions when varying t or s. Q active outer rows
   supply `32Q` message dimensions, while all feedback-zero conditions impose
   at most `s*N/t` linear equations. If `32Q > s*N/t`, there are nonzero inputs whose state stays
   zero. Mixing alone cannot repair a state trajectory that never activates.
4. The one-round implementations now reuse the emission table for `u^T state`.
   Sparse and fixed-masked updates pass the dense reference. A changed A/B
   map requires new generated circuits and verification.
5. Serial full-encoder measurements establish about a 4% gain with the original
   A/B maps. Re-measure the final proof-compatible maps and mixer together;
   routing/outer costs remain even if the mixer becomes free.

The successful Q1 envelope requires the A spectrum and distinct nonzero B
columns, plus the exact cancellation weights. The fresh-mixing argument itself
does not require `BA=0` or `B=A^T`; those identities simplify this implementation
and the retained spectrum reconstruction. They must not be mistaken for the
only possible proof-compatible inner design.
