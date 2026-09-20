# Proof plan

The proof uses state classes `Z` and `U`.  Class `Z` contains the zero state.
Class `U` contains a uniform element of `GF(2^64)^*`.  A fresh nonzero field
scalar refreshes class `U` whenever the pre-randomizer state is nonzero.

Fix an epoch input `X`.

From class `Z`, the output is `Y=X`.  The next state remains zero exactly
when `B(X)=0`.  Inputs that fail to activate the state therefore belong to
the code `ker(B)`.  The proof needs its low-weight spectrum after the packet
permutation.

From class `U`, the output ranges over the affine code

\[
X+\operatorname{im}(A).
\]

The affine maximum has an exact simplification.  Let
`C=im(A)` and fix `0 <= z <= 1`.  Define

\[
F_X(z)=\sum_{c\in C}z^{\operatorname{wt}(X+c)}.
\]

Then

\[
F_X(z)\leq F_0(z)
\]

for every `X`.  To see this, expand the coset sum over the dual code:

\[
F_X(z)=\frac{|C|}{2^{256}}
\sum_{w\in C^\perp}(-1)^{\langle X,w\rangle}
(1+z)^{256-\operatorname{wt}(w)}
(1-z)^{\operatorname{wt}(w)}.
\]

Every coefficient after the sign is nonnegative.  Removing the signs gives
`F_0(z)`.  Thus an adversarial epoch input cannot have a larger exponential
weight moment than the zero coset.  The live-state distribution omits
`Q=0`; consequently its moment is at most `F_0(z)/(2^64-1)`.  This differs
from the exact nonzero-codeword moment only by one codeword term.

The next state is zero only when

\[
Q=B(X).
\]

For fixed nonzero `B(X)`, this event has probability `1/(2^64-1)`.  Its
output equals

\[
(I+AB)X.
\]

Since `BA=0`, the map `I+AB` is invertible: `(AB)^2=0` and
`(I+AB)^{-1}=I+AB`.  A nonzero input cannot cancel both the output and the
next state.

The exact epoch transfer therefore depends on three constituent quantities:

1. the weight spectrum of `ker(B)`;
2. the ordinary weight enumerator of `im(A)`, which bounds every affine
   coset by the lemma above;
3. the weight transfer of the invertible map `I+AB`.

The first exploration may average over sampled constituent maps.  Such an
average is evidence about a design ensemble, not a certificate for a fixed
map.  A final proof must either certify one fixed pair `(A,B)` or include the
sampling of that pair in setup and account for all setup randomness.

## Current proof boundary

The affine input and the live-state constituent are no longer conceptual
obstructions.  The exact sampled-ensemble enumerator proves the existence of
a fixed nested pair with `d(im(A)) >= 40`, `d(ker(B)) >= 4`, and 1,294,555
bits of all-live suppression at the 9% global threshold.  The generic
minimum-distance MILP remains unsuitable for certifying the particular
sampled matrices used by an implementation.

The zero-state transition is now exact as a function of the two coordinate
weights.  The selected parity-lane wiring also proves that every input
supported on at most three packets activates.  The modeled weight-38 prefix
calculation is positive through four active outer blocks.

The occupation-ladder calculation now combines:

1. outer multiplicities and placement with zero-state inputs in `ker(B)`;
2. live-state coset moments, using the zero-coset upper bound;
3. termination probability `1/(2^64-1)` and the nonzero output transfer
   `(I+AB)X`;
4. the transition back from a terminated state to the next activating input.

The resulting common recurrence is positive through the tested
medium-occupation range once the first three live-output moments are used.
The remaining boundary is combinatorial: shared four-block packet groups,
the all-one word, and the transition to a dense-regime argument.

## Deterministic near-one moment

The one-active shell chooses a Chernoff tilt close to one.  At that tilt, a
generic conditional-expectation factor is too expensive because the factor
is repeated in every live epoch.  Full output support removes that loss.

Let `C=im(A)`, let `d(C)>=d`, and assume that every output coordinate of `A`
is a nonzero linear form.  For a uniform nonzero `c` in `C`,

\[
\mu:=\mathbb E[\operatorname{wt}(c)]
=\frac{256\cdot2^{63}}{2^{64}-1}.
\]

For `0<z<1`, the function `w -> z^w` is convex.  Its secant on `[d,256]`
therefore gives

\[
\mathbb E[z^{\operatorname{wt}(c)}]
\leq
\frac{256-\mu}{256-d}z^d
+\frac{\mu-d}{256-d}z^{256}.
\]

This bound is deterministic.  It does not condition on a sampled matrix and
does not pay a multiplicative factor in every epoch.  For an epoch input of
weight `h<d`, the triangle inequality also gives

\[
\operatorname{wt}(X+A(q))\geq d-h
\]

for every nonzero state `q`.  A terminating transition has the same weight
floor and probability `1/(2^64-1)`.

Using `d=40`, the exact one-active zero/live recurrence has 79.671 bits of
modeled aggregate margin at 9%.  Weight 38 is dominant.  The calculation is
floating-point, and the outer spectrum remains modeled.  The fixed sampled
pair has full output support, but its global distance 40 is not yet
certified.  See `receipts/one_active_outer_model_secant.json`.

## Three-moment live bound

The degree-7 fixed map has no linear dependency among one, two, or three
output-coordinate forms.  If `W` is the weight of `A(Q)` for uniform
nonzero `Q`, then, for `j=1,2,3`,

\[
\mathbb E[(W)_j]
=
(256)_j\frac{2^{64-j}}{2^{64}-1},
\]

where `(x)_j` denotes the falling factorial.  These equalities are
deterministic consequences of the coordinate audit.

For every Chernoff parameter `0<z<1`, maximize
`sum_w p_w z^w` over `40 <= w <= 256`, subject to nonnegative `p_w`, unit
mass, and the three factorial-moment equalities above.  This linear program
is a valid upper bound on the live-output moment for any fixed constituent
that also has distance at least 40.  The affine-coset lemma then applies the
same live bound to every epoch input.

At `z=0.95`, the mean-only bound has base-2 logarithm `-3.7149`; the
two-moment and three-moment bounds give `-8.7051` and `-8.9585`,
respectively.  This is enough to remove the artificial medium-occupation
failure of the mean-only recurrence.

## Occupation recurrence

For a fixed number `a` of active outer blocks, the current calculation uses:

1. the regular outer spectrum envelope for every modeled even block weight
   from 38 through 218;
2. an exact coefficient dynamic program for placing the resulting active
   packets among the 8,192 epochs;
3. the exact zero-state activation table and packet-support-four audit; and
4. one common two-by-two zero/live transfer with the three-moment live bound.

Under the current distinct-packet-group assumption, the recurrence has
positive 9% margin at every tested occupation from 1 through 128.  The
occupation-128 margin is 9,929.636 bits.  Lower occupations were already
positive with the weaker mean-only moment, and the three-moment feasible set
is contained in the mean-only feasible set.

The remaining proof cases are now explicit:

1. prove a dominance rule or extend the recurrence to active blocks that
   share a four-block packet group;
2. handle the all-one outer word, which the regular spectrum envelope omits;
3. cover occupations above 128, preferably with a dense-regime argument; and
4. certify distance 40 for a fixed degree-7 constituent, or prove joint
   existence of distance 40, three-coordinate independence, and the packet
   compressor properties.

## Shared-group reduction attempt

For a fixed profile `(n1,n2,n3,n4)`, let `nr` count packet groups containing
`r` active outer blocks.  The four counts determine the exact packet-placement
transfer.  No individual block identities are required.

The exact profile enumerator closes every profile at occupations 2, 4, 8,
and 16.  The distinct-group profile is worst at each of those occupations.
The maximally packed profiles remain positive throughout the tested ladder
after profile-specific optimization where needed.

The natural next reduction would transform any profile into the maximally
packed profile.  A local transformation acts on two packet groups.  If the
groups lie in separated epochs, their comparison contains the transfer
between those epochs.  Entrywise domination must therefore hold for every
reachable intervening transfer.

This domination fails pathwise.  A split pair supports the state sequence

\[
U\longrightarrow Z\longrightarrow U.
\]

The first packet terminates the live state.  The second packet reactivates
it.  After merging, one of those two transitions has no input and the same
state path has zero contribution.  Hence no finite local scalar compares all
state paths.

The obstruction is narrower than arbitrary profile enumeration.  Every
unsupported path contains a termination transition and therefore a factor
`1/(2^64-1)`.  A promising continuation is a termination-stratified
recurrence:

1. prove packet-packing domination for paths with no termination;
2. expose the first termination epoch;
3. restart the packing comparison after the next activation; and
4. sum the resulting segments with the explicit termination factor.

This route preserves the four-count compression and isolates the only state
transition that defeats it.  The termination-fugacity audit, however, shows
that this route is not uniform: the optimized all-quad tilted measure has
about 123 termination transitions.  Termination cannot be charged as a small
exception in the packed regime.

## Region-level packing lemma

The next proof target moves the comparison outside the state-path expansion.
Fix `p`, `z`, and a total occupation `a`.  For a profile
`n=(n1,n2,n3,n4)`, let `R_n(p,z)` be the exact two-by-two transfer after:

1. sampling the Bernoulli outer bits;
2. applying the uniform packet permutation in one region; and
3. composing all 32 epoch transfers, including termination.

Let `pack(a)` contain `floor(a/4)` full groups and one group for the
remainder.  The desired statement is

\[
R_n(p,z)\leq R_{\operatorname{pack}(a)}(p,z)
\]

entrywise for every profile with total occupation `a`.

The quantifiers over the Chernoff parameters are deliberately restricted.
The certificate uses the symmetric envelope `p=1/2`; the optimized
occupation-128 point has surprisal `-log(z)=0.0288`.  Finite grid tests support
the order at `p=1/2` throughout the much wider interval
`0.01 <= -log(z) <= 1` for every profile through occupation 16.  The order is
false for arbitrary `p` and for sufficiently large surprisal, so no proof
effort should target that stronger statement.

This statement is stronger than necessary but has a large payoff.  It permits
the proof to multiply the packed inner moment by all `C(8192,a)` active block
sets.  At occupation 128, that universal substitution retains 6,393.067 bits
at the optimized packed point.

The lemma holds in every exact low-occupation scan and throughout the exact
near-packed scan.  A proof should compare the six elementary transformations
of two packet widths only after averaging their positions among all 32
epochs.  The averaging is essential: the analogous pathwise transformation
is false.

The six transformations are

\[
1+1\to2,\quad 1+2\to3,\quad 1+3\to4,\quad
2+2\to4,\quad 2+3\to1+4,\quad 3+3\to2+4.
\]

They generate a path from every profile to `pack(a)`.  Exact recorded data
contain 270 such source-target edges: 240 in the complete profile lattices
through occupation 16 and 30 in the occupation-128 packed corner.  Every edge
has the desired entrywise order.  The next proof step is therefore one local
region lemma covering these six moves uniformly in the surrounding packet
profile; proving it would replace the profile enumeration by induction.

The phrase "uniformly in the surrounding packet profile" cannot be
strengthened to conditioning on every surrounding packet position.  A direct
conditional-insertion audit produces counterexamples in the live-to-live
entry for four merge moves.  The proof object must instead be the normalized
coefficient of the complete region polynomial.  With `A_m` denoting the
Bernoulli-averaged epoch transfer, define

\[
P(x_1,x_2,x_3,x_4)=
\sum_{c_1+\cdots+c_4\leq64}
\frac{(64)_{c_1+\cdots+c_4}}{c_1!c_2!c_3!c_4!}
A_{c_1+2c_2+3c_3+4c_4}x_1^{c_1}\cdots x_4^{c_4}.
\]

For profile `n`, the exact region matrix is

\[
R_n=
\frac{[x^n]P(x)^{32}}
{(2048)_{|n|}/(n_1!n_2!n_3!n_4!)}.
\]

The next lemma should compare these normalized coefficients under the six
profile moves.  This formulation keeps precisely the joint permutation
average that repairs the conditional counterexamples.

An equivalent formulation exposes an induction on the number of epochs.  Let
`R_E(n)` be the normalized matrix after `E` epochs, let
`k=sum_r n_r`, and for `c<=n` let `ell=sum_r c_r` and
`m(c)=sum_r r*c_r`.  Then

\[
R_E(n)=\sum_{c\leq n}
\frac{\left(64(E-1)\right)_{k-\ell}(64)_\ell}
{(64E)_k}
\prod_{r=1}^4\binom{n_r}{c_r}
R_{E-1}(n-c)A_{m(c)}.
\]

The scalar coefficients form the multivariate hypergeometric law for the
groups assigned to the last epoch.  Complete occupation-8 and occupation-16
audits show the six packing inequalities at every tested depth
`E in {1,2,4,8,16,32}`.  This suggests proving a strengthened order that is
preserved by the displayed hypergeometric convolution.  Such an induction
would use the joint permutation symmetry without conditioning on a fixed
background placement.
