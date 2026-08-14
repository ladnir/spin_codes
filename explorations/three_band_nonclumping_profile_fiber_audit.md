# Three-band incidence constraints on packet-profile fibers

The current construction already uses a fixed three-band tile map. Its
incidence restrictions can refine an outer profile count without changing the
layout distribution. They do not, by themselves, define a smaller set of
nine-entry packet profiles.

## Fixed incidence object

Let

```text
E = Z_256 x {0,...,63}
```

be the 16,384 data-block positions. For `e=(t,l)` and band `b`, define

```text
lambda_0(e) = t,
lambda_1(e) = t + 9l mod 256,
lambda_2(e) = t + 20l mod 256.
```

Regard each position as the hyperedge
`(lambda_0(e),lambda_1(e),lambda_2(e))` in a three-partite hypergraph `H`.
For `F subset E`, let `P_b(F)` partition `F` by equality of `lambda_b`.

The fixed map has four certified properties.

1. Each of the 256 tile labels in each band contains exactly 64 edges.
2. For distinct bands `b,c`, the map `e -> (lambda_b(e),lambda_c(e))` is
   injective. Thus, two data blocks share a tile in at most one band.
3. No four edges induce two disjoint collision pairs in every band. This is
   the Pasch, or intercalate, exclusion.
4. Exactly 281,856 four-edge sets have five of their six block pairs
   colliding. The counts are 79,616, 122,880, and 79,360 according to the band
   containing the single collision.

The second property follows algebraically because every nonzero difference
among `0,9,20` is odd and therefore invertible modulo 256. The exact tile-map
certificate also checks all four properties exhaustively.

## Fiberwise lemma

Let `L` be the 16,384 logical data blocks. The construction samples a uniform
bijection `Pi:L -> E`. Fix a message `m`, and let `S(m) subset L` be its
nonzero data blocks. For an outcome `omega`, define

```text
F(omega) = Pi(S(m)),
tau(omega) = (P_0(F(omega)),P_1(F(omega)),P_2(F(omega))).
```

Let `A(omega)` be the final packet profile. For each profile `a`, define the
profile fiber

```text
Omega_a = {omega : A(omega)=a}.
```

The exact incidence lemma is

```text
Omega_a = disjoint union over realizable tau of
          {omega in Omega_a : tau(omega)=tau}.
```

Every nonempty summand satisfies band balance, pair capacity one, and Pasch
exclusion. Its four-edge subpatterns belong to the certified affine-map
catalogue. Equivalently, for every positive-probability conditioning event
`C`,

```text
Pr[tau is not H-realizable | A=a, C] = 0.
```

This is a support statement. It remains valid after conditioning on a packet
profile because an impossible incidence type remains impossible. It makes no
claim that the surviving incidence types retain their unconditional
probabilities.

The lemma constrains an augmented fiber `(a,tau)`, not `a` alone. A packet
profile forgets which data blocks supplied its active lanes. Therefore, the
incidence restrictions yield a numerical improvement only when the outer
enumerator retains `tau`, or a certified relaxation of `tau`.

## Exact probabilities before profile conditioning

Fix `m` and every construction variable sampled before `Pi`. If `|S(m)|=s`,
then `F=Pi(S(m))` is a uniform `s`-subset of `E`. For any family `U` of
`s`-subsets,

```text
Pr[F in U | m, pre-Pi variables] = |U| / binom(16384,s).
```

For labeled active blocks, the denominator is the falling factorial
`(16384)_s`. The labeled and unlabeled conventions must not be mixed.

Band balance gives an exact one-band occupancy law. If `n_u` active positions
use tile `u` in one fixed band, then

```text
Pr[(n_u)_u] = product_u binom(64,n_u) / binom(16384,s).
```

The three band occupancies are correlated. Their marginal laws cannot be
multiplied.

For four fixed active blocks,

```text
Pr[Pasch] = 0,
Pr[K4-minus-one-edge] = 281856 / binom(16384,4).
```

The three sparse-band numerators are the certified counts above. For `s>4`, a
safe union count for subsets containing a dense four-block motif is

```text
281856 * binom(16380,s-4).
```

This count may overcount an `s`-subset containing several motifs. That
overcount is safe for an upper bound.

These formulas apply before conditioning on `A=a`. In general,

```text
Pr[F in U | A=a] != |U| / binom(16384,s).
```

To use a motif count inside one profile fiber, the verifier must count the
joint set `{omega:A(omega)=a, F(omega) in U}` or prove the required
independence. No current certificate proves such independence.

## Existing certificate inventory

`scripts/certify_three_band_tile_map.py` is the primary incidence certificate.
It fixes the affine map, checks exact band balance and pair capacity, excludes
Pasch configurations, and counts the K4-minus-one-edge family.

`scripts/certify_packet8_constant_min_support.py` already uses pair capacity
fiberwise. It combines the restriction with exact fixed-band outside distances
to exclude nonzero constant-packet words below 40 active packets. This is a
special profile-fiber theorem, not a general nonclumping law.

`scripts/certify_three_band_small_motifs.py` enumerates all pairwise-compatible
triples of set partitions for sizes two through six. It deliberately includes
patterns absent from the affine map. Its rational moment bounds are safe broad
envelopes, but they do not certify the exact affine-map fiber.

`scripts/certify_three_band_pair_core.py` assumes that every collision tile has
size two and relaxes the pattern to three independent matchings. It drops
connectivity and the prohibition on reusing a pair. It is a valid envelope
only within its stated size-two collision case.

`scripts/certify_three_band_star.py` handles all one-tile stars of sizes two
through 64. Pair capacity proves that the other two bands are collision-free.
The certificate includes exact rational puncture and graph factors.

`scripts/certify_three_band_peeling.py` certifies finite algebra used by the
cluster-peeling argument. Pair capacity gives width at most 17 through the
relevant component range. The affine map improves the width to seven for at
most 64 edges. The same certificate exhibits a 64-edge, 4-regular pattern, so
balanced dense patterns remain possible.

`scripts/certify_three_band_exact_length.py` supplies the distinct-hole,
puncture, and graph bookkeeping. These objects are not edges of `H`. Any
profile-fiber verifier must keep them separate or use a proved unpunctured-to-
physical replacement bound.

## Verifier contract

A theorem-facing verifier needs the following immutable inputs.

1. The construction manifest must bind `256` tiles, `64` lanes, slopes
   `(0,9,20)`, and the uniform block-assignment law.
2. A layout-conformance receipt must connect the production encoder to this
   edge map. A certificate for a hardcoded map does not establish that the
   implementation uses it.
3. The verifier must reconstruct all 16,384 edge triples canonically. It must
   check balance, pair-projection injectivity, and the four-edge catalogue.
4. A fiber record must identify the active logical data blocks, their labeled
   or unlabeled counting convention, and the retained incidence statistic.
5. A joint outer bound must connect local BCH band weights, coordinate
   bijections, tile clusters, graph holes, and the final packet profile.
6. Every multiplicity must be an exact integer. Every moment factor must use
   exact rational or outward interval arithmetic.
7. The receipt must bind the map, motif catalogue, exact-length inputs, profile
   fiber, and every local inequality used in the aggregation.

The smallest useful first contract retains the three compatible collision
partitions and per-cluster local BCH weight data. It then sums only compatible
partition triples. The sizes-two-through-six enumerator already provides a
safe relaxation of this contract. Incorporating the exact Pasch exclusion and
the affine four-edge catalogue can only decrease that sum.

## Proof obligations

Before this direction enters the main theorem, a proof must establish:

1. the sampled block assignment is independent of the message and of all
   variables fixed before the incidence count;
2. the active data-block set is defined before that assignment;
3. the production layout realizes the certified affine map for every setup;
4. the outer sum counts setup randomness with its normalized probability,
   rather than as unnormalized multiplicity;
5. the augmented fiber decomposition preserves every profile outcome exactly;
6. the local moment bound is uniform over all complete codewords represented
   by one incidence type;
7. punctured data rows and graph holes are included exactly once; and
8. the final marginalization over incidence types has no omitted patterns and
   no unjustified division by pattern multiplicity.

The conditional inner lemma remains separate. It must hold for every complete
outer word after the outer proof has used the incidence refinement.

## Fatal shortcuts

- Pair capacity one does not mean that every collision tile has size at most
  two. A one-tile star can have 64 data blocks.
- Tile nonclumping does not imply packet nonclumping. Coordinate bijections and
  lane grouping still determine which active bits share a packet.
- Band balance does not make the three band partitions independent.
- A fixed affine map cannot be replaced by three fresh random hash functions.
- Unconditional motif probabilities cannot be reused after conditioning on a
  packet profile.
- Pasch exclusion does not exclude dense four-block patterns. The map contains
  281,856 K4-minus-one-edge sets.
- A global motif count is not a deterministic per-fiber cap. Larger active
  sets require overlap-aware counting or a safe union bound.
- The small-motif certificate cannot be extrapolated beyond six blocks.
- Maximum tile degree does not bound connected-component size.
- The graph block is not a fourth part of the certified hypergraph.
- Sampling a new no-clumping layout, or conditioning the old layout on a
  no-clumping event, would change the construction and its probability law.

The recommended theorem step is therefore modest: add the augmented incidence
type to the outer profile fiber, replay the sizes-two-through-six sum with the
actual Pasch and four-edge exclusions, and measure the rigorous gain. A direct
constraint on the nine profile counts is not presently justified.
