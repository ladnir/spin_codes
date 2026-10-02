# Four-bit dense-range proof and hill-climbing

The first [complete four-bit certificate](FIRST_CLOSURE.md) now proves
relative distance greater than 0.5% with more than 43.744 bits of margin.
It was independently replayed and assembled at 256-bit precision.
The first [hill-climb to 1%](HILL_CLIMB.md) is still incomplete.
The [boundary diagnosis](BOUNDARY_DIAGNOSIS.md) finds removable shuffle
comparison slack: four selected terms that failed at 5% now pass a
sharper exact comparison. This special-composition result is not a
whole-code certificate.
The [posterior comparison](POSTERIOR_COMPARISON.md) extends the sharper
bound to mixed packet types and completely enumerated local cells. It is
available to the atlas through `--posterior-shuffle`; broad unresolved
cells still need further coverage work.

This directory first closes at a deliberately reduced distance, then
increases the cutoff. It does not change the encoder. A partial atlas, a selected
composition, or a floating-point screen is not a whole-code certificate.

The setup is the independent-row ensemble in [the bridge](../BRIDGE.md):
K=2^20, N=2^21, BCH[256,128], four adjacent rows per group, 2048 groups,
256 regions, independent row, regional packet, and packet-lane shuffles,
and two-update IMT(128,19). The state starts at zero and is not flushed.
The maps and BCH spectrum caps are authenticated by the existing loaders.
The bounds concern setup sampled once, simultaneously for every
nonzero message; the comparison measures below do not change that setup.

The existing [sparse proof](../LOW_OCCUPANCIES.md) covers 1 through 58
active groups at bad-weight cutoff 209715. Its conservative summed upper
is 6.786362e-14. It also applies to any smaller cutoff. The new dense
calculation covers 59 through 2048 at the smaller cutoff. For further
hill-climbing, a complete dense upper below 2^-60 suffices to combine
the two ranges with more than 40 bits of margin.

The [lower-cutoff sparse extension](SPARSE_HANDOFF.md) now covers every
occupancy through 64 at cutoff 104857 (5% relative distance). A dense
cover can therefore start at 65 for any cutoff at most 104857. This is
an alternate proof handoff, not a completed whole-code 5% certificate.

## Comparison measure

For a BCH row, `row_mixture.envelope` bounds the expected count at every
particular shuffled binary word by a positive mixture of five product
measures. Their bit probabilities are 0, 1, 1/2, theta, and 1-theta.
It checks the pointwise inequality for all 257 weight shells with exact
rationals. These are counting-measure coefficients, not probabilities of
choosing a new outer code.

Four row components give 70 unordered component types. A type's
coefficient includes the number of row orderings. Its packet-weight law
pi_i has five entries, for weights 0 through 4. Independent uniform lane
shuffles make words of a given packet weight equiprobable. Only the
four-zero-row component has inactive label; comparison inputs that happen
to be zero retain their component's active label.

Let n_i count group types, with sum n_i=G=2048. In one region, the type
laws are independent but not identical, before the uniform packet shuffle.
For positive category tilts t_j, set

    Z_i = sum_j pi_i(j) t_j,
    f_i(j) = pi_i(j) t_j / Z_i,
    x_j = sum_i (n_i/G) f_i(j),
    w_j = x_j / t_j.

The tilted heterogeneous-shuffle inequality from
[the two-bit comparison proof](../../two_bit/MIXTURE_BOUND.md) applies
with five categories. For a realized category-count vector a, its
pointwise density loss is bounded by

    R(a) = G^G / G! * product_j a_j! / a_j^a_j,

where the factor for a_j=0 is 1. We use either its exact maximum under
category-count caps or a separable tangent upper
`R(a) <= B product_j tau_j^a_j`. The latter multiplies reference weights
by tau_j and pays B per region. Its integer anchors are witnesses, not
assumptions on the actual category counts. There are 256 independent
regional shuffles, so all regional losses are paid 256 times.

The group coefficient becomes `c_i Z_i^256`. A multinomial generating
function sums all labeled component assignments in each cell with at
least 59 active groups. No sampling of this composition space is used
for a coverage claim.

## Inner moment

`kernel.py` bounds the output-weight Laplace transform under the iid
packet comparison measure. It enumerates all feedback characters and all
nonzero expansion-image packet-weight histograms for the actual maps.
Its two states indicate whether the IMT state is zero or nonzero.

For zero state, the zero-feedback contribution is computed by Fourier
inversion; the other contribution is the total moment minus that value.
For nonzero state, let h bound the output moment for every possible state.
Two independent updates give the mixture `(1/4) identity + (3/4) uniform
nonzero refresh`. The transition to zero is at most

    (1/4) c + (3/4) h / (2^19-1),

where c bounds the weighted feedback cancellation on the identity branch.
The verifier takes the best of the total moment, a second-moment atom
bound, and a tilted Fourier atom bound. The transition to nonzero is
bounded by h; double-counting the zero contribution here is safe.
Output is emitted before the refresh and feedback update.

Raising this nonnegative transition bound to 16384 steps, starting at
zero, bounds the full iid reference moment. The output Chernoff factor
is `exp(lambda * cutoff)`. Scalar reference-mass normalization is paid
per packet (524288 packets), not per bit.

## Cover and replay

`atlas.py` covers x_1 through x_4 with a binary split tree; x_0 is their
complement. Every terminal cell is either proved empty or bounded outward.
Float optimization proposes witnesses only. Replay reconstructs the
whole tree, rechecks every emptiness claim, and recomputes the bounds.

The refinements keep the cover from paying unnecessary slack:

- Supporting directions from a floating convex hull are rounded to
  rationals. Their intercepts are recomputed on every exact vertex before
  they can prune a cell. Exact edge clipping also enumerates the vertices
  of the box intersected with these verified halfspaces.
- Near the least-active component, `boundary.py` enumerates every feasible
  integer composition. If a work limit is hit, it returns no conclusion,
  never a partial list. Replay regenerates the complete list.
- Alternative positive category tilts can be selected cell by cell.
  Rational affine majorants of the new reference weights are checked
  against every component type. This avoids using a tilt meant for light
  packets in a cell containing a significant all-one component.
- An optional affine moment bound couples the inner moment to the outer
  count. The exact unnormalized moment is a polynomial with nonnegative
  coefficients. After tangent upper bounds on log reference weights, its
  logarithm is convex in the cell coordinates. Bounds at every vertex of
  the clipped box (or all 16 un-clipped vertices for an earlier witness)
  therefore give a valid affine majorant on the whole cell. Convexity is
  asserted for this exact polynomial, not for the two-state envelope;
  the envelope is used only to bound the vertex values.

### Coupling alternate comparison weights

The optional `--coupled-alternates` search applies the affine moment
bound to alternate category tilts as well as the base tilt. This avoids
maximizing the inner moment and the outer count at unrelated compositions.
It does not change the code or its setup distribution.

Fix a cell in the base coordinates x=(x_1,...,x_4), and an alternate
positive tilt t'. For group type i, define Z'_i=sum_j pi_i(j)t'_j and
v_ij=pi_i(j)/Z'_i. For a composition n, its alternate reference weight is
w'_j=sum_i (n_i/G)v_ij. The base features f_i and active labels a_i were
defined above. Let rho=q_min/G. Given rational eta_j in Q^4 and mu_j>=0,
the verifier computes

    h_j = max_i (v_ij - eta_j dot f_i(1:4) + mu_j a_i).

Every composition with at least q_min active groups then satisfies

    w'_j <= h_j + eta_j dot x - mu_j rho =: W_j(x).

The verifier checks the intercept against all 70 component types. It
does not trust a floating linear-program optimum. Earlier alternate
witnesses maximized W_j separately over the cell. The coupled witnesses
instead retain W_j as an affine function of the same x used for counting.

For any positive rational anchor b_j, the exponential tangent gives
W_j(x) <= b_j exp((W_j(x)-b_j)/b_j). The exact unnormalized input moment
has nonnegative coefficients. After this substitution, its logarithm
is a log-sum-exp of affine functions of x, hence convex. The existing
vertex argument bounds it by an affine function over the clipped cell.
The outer generating function retains that affine function and uses
the matching alternate coefficients c_i (Z'_i)^256. Thus both parts
refer to the same composition.

Saved witnesses contain the alternate tilt, rational weight duals,
positive anchors, moment slope, and outer-count dual. Replay reconstructs
all weight majorants and bounds every vertex with outward arithmetic.
The search flag affects proposals only; replay recognizes the witness
without requiring that flag. Earlier base-tilt witnesses remain valid.

The same split-tree engine is reused from the two-bit proof, but its
geometry, packet law, kernel, and group count are not reused. The two-bit
certificate itself is unchanged and is not evidence of four-bit closure.

## Commands

Run from the repository root with Python, NumPy, SciPy, and python-flint:

```sh
python -B -m unittest discover -s research/workstreams/permutation_locality/independent_rows/dense_closure -p 'test_*.py'
python -B research/workstreams/permutation_locality/independent_rows/dense_closure/atlas.py --threshold 10485 --input-tilt 1/32 1/64 1/128 1/256 --max-cells 300 --max-depth 64 --output tmp/four-bit-dense-plane-half.json
python -B research/workstreams/permutation_locality/independent_rows/dense_closure/atlas.py --resume tmp/four-bit-dense-plane-half.json --max-cells 300 --max-depth 64 --output tmp/four-bit-dense-plane-half-next.json
```

Replay accepts a complete, outward-verified atlas only:

```sh
python -B research/workstreams/permutation_locality/independent_rows/dense_closure/atlas.py --replay tmp/complete-four-bit-dense.json --precision 256
```

The last filename is illustrative, not an assertion that such a completed
artifact currently exists. Raw atlas data stay under `tmp/`.

Once a complete dense atlas exists, `assemble.py` independently replays it
and adds the conservative upper from the existing sparse theorem. It
rejects incomplete/screening-only records, an occupancy gap, a different
schema, and a cutoff beyond the sparse theorem's range. To raise the
distance, use `atlas.py --retarget OLD --threshold NEW ...`; every old
accepted witness is checked again at the new cutoff before being reused.
The same operation may raise `--minimum-groups` to restrict the domain;
the assembler accepts that handoff only when a matching sparse theorem
covers all smaller occupancies.

The regression tests cover exhaustive eight-bit input/state transitions,
uniform packet-lane shuffles, all row-component multiplicities, bounded
integer composition enumeration, exact four-dimensional polytope clipping,
and alternate-weight duals. A passing test suite does not substitute for
complete coverage and outward replay on the actual maps.
