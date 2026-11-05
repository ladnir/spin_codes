# Directed bounds for explicit integer-weight bands

One directed coefficient bank can bound several active-group weights.
The bound below evaluates normalized coefficients at every integer weight;
it does not interpolate between selected weights. It groups the subsequent
transfer calculation and uses the midpoint rectangle counts.

Fix the full-field gamma construction with width 64, state dimension 30,
and 256 regions. Its setup distribution and persistent state are those in
[the all-pair composition proof](HD_FIELD128_ALL_PAIR_COMPOSITION.md).
Put `P=4096`, `q=2^128`, and fix positive proof parameters `0<z<1` and
`0<alpha<=1` for the whole band. For target distance 40%, put `L=838860`.
Define

    beta = q^160/(q-1)^32,
    sigma = alpha log2(beta),
    c = -alpha [log2(q-1)+L log2(z)].

All probabilities below are over this construction's setup randomness.
Choosing the proof parameters or coefficient markers does not sample or
condition that setup.

## The profile bound used as input

For an indexed binary support of length `2P`, let r be its weight and
let a and b count full pairs in the even and odd matchings. Write
`N(P,r,a,b)` for the exact number of supports with this profile, as defined
in [the exact-count note](HD_FIELD128_EXACT_PROFILE_COUNTS.md).

For each original weight r and feasible full count f, let `R(r,f)` be a
nonnegative 2-by-2 regional upper matrix. It must bound the regional
operator from the all-pair proof after normalization at

    single = r-2f,  full = f,  empty = P-r+f.

Here `e0=(1,0)` and `one=(1,1)^T`. The existing proof gives

    Pr[there is a forbidden codeword with profile (r,a,b)]
      <= N(P,r,a,b) 2^(c+sigma*r)
         e0 [R(r,a) R(r,b)]^128 one.                 (1)

The outer counting factor and projective divisor occur only in
`2^(c+sigma*r)`. Supports are counted outside the fractional power.
The matrix order is even then odd, with no state reset.

## Cells in deficit coordinates

First fix an integer band `1<=first<=last<=P`. Split its weights by
parity: write `r=2u+eps`, where `eps` is 0 or 1. Define

    d=u-a,  e=u-b.

The full marginal domain is `0<=d,e<=u`. In particular, d and e may be
zero. For even r, these boundary values include the singleton-free
matching faces; they must not be omitted.

A cell consists of a fixed parity eps and three inclusive integer
intervals `U,D,E`. Its integer triples satisfy

    u in U,  first<=2u+eps<=last,
    d in D intersect [0,u],  e in E intersect [0,u].    (2)

For a permitted u, write `D_u=D intersect [0,u]` and similarly `E_u`.
Skip that u if either intersection is empty. Its corresponding full-count
intervals are

    I_u = [u-max(D_u), u-min(D_u)],
    J_u = [u-max(E_u), u-min(E_u)].

These intervals include marginally feasible profiles with no supports.
Including such profiles is harmless. The raw midpoint count is zero on
impossible profiles, including the doubly singleton-free profile at a
nonconstant even weight. The alternating profile `r=P,a=b=0` instead
has count 2 and is handled explicitly by the midpoint helper.

Let `C_r(I,J)` be the exact integer sum of raw midpoint upper counts
over the ordered rectangle `I x J`, including its exceptional cases.
It is computed by [the rectangle-prefix formula](HD_FIELD128_MIDPOINT_BOXES.md).
Define the cell support-count upper bound

    C_cell = sum_permitted u C_(2u+eps)(I_u,J_u).      (3)

This is an exact sum of count upper bounds, not the exact support count.
A valid alternative cap on the same cell may be combined with (3) by
taking their minimum. No extra number-of-weights or profile multiplier
is needed: (3) already sums every relevant integer weight and profile.

## Absorbing the weight-dependent outer factor

At every permitted integer pair `(u,d)`, form

    W(u,d) = 2^(sigma*(2u+eps)/256) R(2u+eps,u-d).

For the cell in (2), take componentwise maxima

    M_D = max W(u,d) over u with D_u,E_u nonempty and d in D_u,
    M_E = max W(u,e) over u with D_u,E_u nonempty and e in E_u.

Return a zero bound before taking maxima if the cell has no integer
triples, or if `C_cell=0`. Otherwise both maxima are over nonempty finite
sets. The cell contribution is at most

    2^c C_cell e0 (M_D M_E)^128 one.                 (4)

To prove (4), fix an arbitrary triple `(u,d,e)` in the cell. Put
`r=2u+eps`, `a=u-d`, and `b=u-e`. Scalar extraction gives

    e0 [W(u,d) W(u,e)]^128 one
      = 2^(sigma*r) e0 [R(r,a) R(r,b)]^128 one.

There are 256 regional factors, so the exponent is exactly `sigma*r`.
Nonnegative matrix multiplication is monotone in every entry. Thus the
left side is at most `e0 (M_D M_E)^128 one`. Apply (1), sum over the
cell's integer triples, and use (3). This proves (4).

This proof permits the maxima to come from different weights or matrix
entries. Such choices can weaken the bound but cannot invalidate it.
The support count is symmetric in a and b; the matrix product need not
be. Retain the ordered even and odd factors.

## Required coefficient evaluation

Fix a positive coefficient marker `(x,y)`. For singleton count s and
full count f, its histogram normalizer is

    H(s,f) = P!/[s! f! (P-s-f)!].

The directed alias bank supplies a nonnegative upper numerator at the
actual residues `(s mod grid,f mod grid)`. Its positive fallback may also
be used. Divide each entry by `H(s,f) x^s y^f` at that same integer
histogram. Entrywise minima across valid markers are permitted.
There is no factor `2^s` in H: source orientations are already counted
by N, and the regional operator averages the fair pair swaps.

The implementation must then apply the weight-dependent scalar and take
the cell maxima. It must not replace this sequence by separate maxima
of alias numerators and inverse normalizers and claim the same bound.
That replacement can only be used as a separately justified, generally
looser upper bound. Endpoint samples or a monotonicity assumption in r
cannot replace evaluation at every integer histogram.

There is a useful exact normalization recurrence. For fixed eps and d,
put `s=2d+eps`, `f=u-d`, and `h=P-u-d-eps`. Define

    S(u,d) = 2^(sigma*(2u+eps)/256) / [H(s,f) x^s y^f].

When both u and u+1 are in the marginal domain,

    S(u+1,d)/S(u,d) = 2^(sigma/128) (f+1)/(h*y).     (5)

The denominator h is then positive. Equation (5) follows from
`H(s,f+1)/H(s,f)=h/(f+1)`. It reuses normalization work, not coefficient
values: the residue numerator at each new u must still be selected.
Initialization and updates require upward-directed arithmetic.

## Coverage, extension, and checks

The two parity domains in (2) partition the requested integer weights.
A certificate must check that its cells partition each domain without
missing a triple or counting one twice. Inclusive lattice splits use
`[lo,k]` and `[k+1,hi]`. An exhausted budget leaves uncovered cells
pending; it cannot produce a complete-band claim. Sum the completed
cell bounds with directed arithmetic only after the requested partition
is complete. Events in setup space may overlap; this final sum is a
union bound, not an assertion that those events are disjoint.

The same argument extends to `last<2P` by replacing `[0,u]` in (2) with
`[0,min(u,P-u-eps)]`. Every regional matrix and outer scalar remains at
the original weight r. Only the count helper uses the bijection
`(r,a,b) -> (2P-r,P-r+a,P-r+b)` when r exceeds P. No symmetry of the
encoder or transfer matrices is inferred. The all-active weight `r=2P`
may remain in the separately certified terminal range.

Before production use, independent tiny tests should check:

- Exhaustive indexed supports map into exactly one parity cell, including
  singleton-free faces, alternating words, and upper-half weights.
- Cell counts match direct sums of raw midpoint caps and dominate exact
  support counts; splitting a cell preserves its count sum.
- With arbitrary positive, noncommuting matrices, (4) dominates the
  explicit integer-profile sum of (1), and singleton cells agree.
- The recurrence (5) agrees with direct rational histogram normalization,
  including the last valid update; residue numerators vary independently.
- Empty or zero-count cells return zero, while partial covers return no
  complete-band bound. Directed scale factors and the final union are
  checked at more than one precision.

For B weights, m markers, and O(r) full-count values per weight, direct
preparation uses O(Bmr) normalized coefficient operations and O(Br)
count preparation. It needs no new coefficient FFT per weight. After
preparation, a cell needs O(B) exact rectangle queries and one positive
matrix power. These are arithmetic-operation counts; integer sizes and
directed precision affect their cost. The grouped bound may be loose,
so no particular band width or numerical success is asserted here.
