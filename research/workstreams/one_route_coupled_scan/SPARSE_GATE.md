# Two bounded sparse proof gates

This note specifies two finite calculations for width four under the actual
paired route. Neither calculation covers every one-group word, let alone
the full code. Their purpose is to test the continuation before a larger
profile search. Section 6 records their first bounded evaluations;
Section 7 records the subsequent closures and the larger classes now covered.

Use the lookahead transfer from [RANK_CONTINUATION.md](RANK_CONTINUATION.md),
with the proposed pair-antipodal state cycle. The formulas also apply to
any other fixed four-cycle, but no such change is implicit. Fix rational
1 <= tau <= 256 and 0 < z <= 1. Write e=floor(0.35N)=734003.

## 1. A region with exactly one nonzero input

One active outer group supplies one coordinate in each region where its
source support is nonzero. The other member of its source pair is zero.
A uniform permutation of the 4096 pairs and the independent swap therefore
place this coordinate uniformly among all 8192 regional positions.

There are m=2048 width-four steps per region. The active step is uniform
among these m steps, and its coordinate is independently uniform among
the four positions. This statement is exact for either source matching
parity and any selected group. It samples no replacement coordinates.

Let E be the augmented empty-input transfer matrix. Let A be the average
of the four augmented singleton-input matrices. Their states are pairs
(H,H_1), where H and H_1 are current and proposed next state supports.
The exact regular-region matrices are

    B_0 = E^m,
    B_1 = (1/m) sum_{ell=0}^{m-1} E^ell A E^(m-1-ell).

B_0 represents an empty region; B_1 represents exactly one nonzero input.
The latter is an insertion average, not a Chernoff bound on the number of
inputs. It cannot introduce an all-empty input history. It is also the
specialization of the paired coefficient-extraction formula with
(n_0,n_1,n_2)=(4095,1,0).

### Free terminal state

The last physical step uses the ordinary, non-lookahead terminal vector
f_I from the continuation note. Put f_0=f_empty and let f_1 be the average
of the four singleton terminal vectors. The two last-region vectors are

    b_0 = E^(m-1) f_0,

    b_1 = (1/m) [
        sum_{ell=0}^{m-2} E^ell A E^(m-2-ell) f_0
        + E^(m-1) f_1 ].

The last summand places the nonzero input in the final physical step.
All terminal states are summed. No final reset or output beyond N is
introduced. Regular region matrices retain lookahead across region
boundaries; only the final physical step loses its lookahead charge.

Let l select every augmented initial state (empty,H_1), with weight one.
Compatibility at the first step enforces its actual next support.

## 2. Gate A: one group with full source support

Fix one outer group and require all its 256 output coordinates to be
nonzero. Its shortened message dimension is 128 over F, for every
invertible realization of its local outer mixers. Other groups are zero.

The rank continuation gives the following bound for this exact support:

    tau^127 z^(-e) l^T B_1^255 b_1.

There are 8192 choices of the group. Thus

    U_full = 8192 tau^127 z^(-e) l^T B_1^255 b_1

bounds the setup probability of a forbidden word in this class. The
uniformity of the continuation makes the bound valid conditional on any
outer realization. It does not give a conditional guarantee for other
source profiles or the whole code.

A first successful checkpoint would be U_full < 2^-40. This covers only
one group with exact full source support. In particular, it says nothing
about a one-group message whose outer word has zeros.

## 3. Gate B: 33 RS columns, all four components

Choose one group and a set C of 33 RS output columns. Require all four
components in each column of C to be nonzero, and every other column to
be zero. At a zero column, invertibility of its 4-by-4 mixer forces all
four constituent RS symbols to be zero.

Each RS[64,32] row shortened to the 33 selected columns has dimension one.
The group's shortened dimension is therefore exactly four over F. The
nonzero coordinate requirements do not add zero equations; impossible
branches may be discarded. For this exact support, the initial rank
factor is tau^3.

Each selected RS column consists of four consecutive regions with one
nonzero input each. An unselected column consists of four empty regions.
Define regular-column matrices and last-column vectors

    C_0 = B_0^4,           C_1 = B_1^4,
    c_0 = B_0^3 b_0,       c_1 = B_1^3 b_1.

Summing every choice of the 33 columns while preserving their order gives

    U_33 = 8192 tau^3 z^(-e)
      l^T [u^33] (C_0 + u C_1)^63 (c_0 + u c_1).

Matrix products need not commute. The coefficient extraction retains
column positions and the intervening empty intervals. It already counts
all binomial(64,33) choices: do not multiply by that binomial again.
The only external selection multiplicity here is the 8192 group choices.

This class excludes supports with zero components in their selected
columns, fewer or more active columns, and multiple active groups.
It is a proof gate, not an enumeration of all minimum-weight outer words.

The coefficient can be evaluated by a degree-33 positive recurrence over
the 63 regular columns. If row-vector polynomials are used, start V_0=l^T
and update, in chronological order,

    V'_j = V_j C_0 + V_(j-1) C_1,

with missing indices zero. After 63 columns the desired scalar is

    V_33 c_0 + V_32 c_1.

All operations preserve nonnegativity. Any truncation discards only
degrees above 33; it does not replace the position-dependent coefficient
by a maximum over column profiles.

## 4. Bound long empty stretches without losing the quiet state

The matrices above are finite but contain long powers of E. A positive
potential can bound these powers while retaining the quiet state exactly.

An empty-input step has zero rows unless its incoming augmented state is
either (empty,empty) or a pair of nonempty masks. Let P project onto those
states. Then E=P E. Define K=P E P. Its quiet state has one self-loop of
weight one. It cannot connect to the live block, whose states have both
masks nonempty. Write that live block as K_live. Consequently

    E^n = K^(n-1) E,       n>=1,

and the quiet and live parts of K evolve separately.

Choose a positive rational vector v on the live states and a positive
rational lambda such that

    K_live v <= lambda v

componentwise. This is a finite inequality to check with outward rounding
or exact rational arithmetic. Do not infer it from an approximate
eigenvalue. A useful choice normally has lambda<1; the inequalities below
remain valid without that assumption.

Embed v and its coordinatewise reciprocal on the live states, with zeros
elsewhere. Let P_0 project onto (empty,empty), and define

    U = P_0 E,
    V = v (1/v)^T E.

Positivity and the potential inequality give the componentwise bound

    E^n <= U + lambda^(n-1) V,       n>=1.

The proof uses K_live^k[i,j] <= lambda^k v[i]/v[j], including k=0.
Keep E^0=I exactly. The constant U is the genuinely quiet contribution;
it must not be merged into a spectral radius at least one for the entire
empty-input operator.

The singleton activation rule gives

    U A U = 0.

Indeed, U enters A with empty actual state. A singleton input activates a
nonempty state, so it cannot reach the quiet row accepted by the right U.
This removes the spurious all-quiet leading term from an active region.

For m>=3, define S_n(lambda)=sum_{i=0}^{n-1} lambda^i. Substituting the
empty-run bound into the insertion average gives the explicit positive
domination

    B_1 <= (1/m) [
        A(U+lambda^(m-2)V) + (U+lambda^(m-2)V)A
        + S_(m-2)(lambda) (UAV + VAU)
        + (m-2)lambda^(m-3) VAV ].

The omitted term (m-2)UAU is zero. The two boundary insertions are kept
separately; they must not be divided by an extra factor m. Empty regions
satisfy B_0 <= U+lambda^(m-1)V.

One can similarly dominate b_0 and b_1 from their exact formulas, keeping
zero-length empty powers and the final-step vector separate. Substituting
any componentwise upper bounds into the sparse gates is valid because
all subsequent products and coefficients are nonnegative. These bounds
are optional: if the rank-one live envelope loses too much, tighten that
specific interface rather than altering the source multiplicity.

The dominated calculation need not use dense regional matrix powers. Write
U=u r^T and V=v s^T. The displayed upper bound for B_1 has left factors in
the four-dimensional span of Au, Av, u, and v, and right factors in the
span of r^T, s^T, r^T A, and s^T A. The B_0 bound uses subsets of the same
spans. Thus regional and column products can be composed through positive
4-by-4 matrices of inner products after the local potential has been
checked. This is an algebraic compression of the stated upper bounds,
not a new relaxation or a reason to skip boundary vectors.

## 5. What to record from the first calculation

For each attempted gate, retain the chosen tau, z, positive vector and
verified lambda, the resulting upward-rounded log2 bound, and whether the
exact or dominated regional operator was used. Record the last-region
boundary treatment explicitly. A numerical failure is a failure of the
chosen bound, not a low-distance witness.

If both gates close, their union costs at most U_full+U_33. They are still
only two exact source-profile classes. The next missing classes should
be selected by a structural explanation, particularly partially supported
columns and unions of neighboring groups. No full certificate follows
from passing these gates.

## 6. Selected finite checks, 2026-10-20

The support-only lookahead screen did not close Gate A. At tau=256, the
five tried markers z=0.99, 0.994, 0.996, 0.998, 0.999 all gave an upper
bound larger than one. The smallest diagnostic log2 bound was about
2264.59 at z=0.994. This was floating-point evaluation of a relaxation,
not a low-distance witness or an optimization over every marker.

The [finite-state-rank envelope](EMPTY_SPAN_RANK.md) closes Gate A without
changing the encoder. For width-four antipodal wiring, use

    tau=256, z=499/500, lambda=993/1000.

The checker proposes a positive vector by a bounded float iteration and
then verifies all 32 live-state inequalities with exact fractions. It also
constructs the singleton active-step envelope exactly. Positive products in
the final two-state calculation round upward at both 80 and 120 decimal
digits. Both runs give

    P_full < 1.897 * 10^(-256) < 2^(-849).

This sums all 8192 choices of the active group and permits every terminal
state. It is conditional on any allowed outer realization. It covers only
the group's exact full source support, not its partially supported words.

The same two-state envelope does not close Gate B. Its computed upper bound
is approximately 1.77316 * 10^31, which is uninformative as a probability.
Here empty columns use E_env^(4m), active columns use B_env^4, and a
degree-33 positive recurrence sums the 64 column positions. The input rank
factor is tau^3. No additional binomial factor is used.

Reproduce both calculations with
`theory/check_sparse_full_gate.py`. Despite its name, that command reports
the still-open 33-column screen as well. It generates no result files.

The next discriminating task is this 33-column class. Its global message
dimension is four, but the current envelope retains only the local
state-span budget and a scalar global rank potential. It can be useful
to retain the remaining global rank explicitly across active steps, or
to sharpen those steps' paired zero/nonzero choices. This is a proposed
refinement, not an established explanation of the entire remaining loss.

## 7. Separate rank and probability costs: wider sparse closures

The [separated-weight lemma](RANK_PRICE_SEPARATION.md) permits tau>256 in
the empty-span Bellman rule. Independent equations cost 1/tau, whereas
exceptional coefficient roots still cost 1/256. The ordinary active matrix
must therefore use charge parameter min(tau,256). This is a proof change,
not a change to the code or setup distribution.

At tau=65536, z=499/500, lambda=993/1000, all 32 positive-vector inequalities
pass exact checking. Use the two-state B_1 from [EMPTY_SPAN_RANK.md](EMPTY_SPAN_RANK.md),
B_0=diag(1,lambda^2048), and its free terminal vector t_free. For every
h=33,...,64, a one-group support with h fully populated columns has
shortened dimension 4(h-32), independently of the allowed outer realization.
The union over all such supports is bounded by

    8192 z^(-e) sum_(h=33)^64 tau^(4(h-32)-1)
      (1,0) [u^h] (B_0^4 + u B_1^4)^64 t_free.

The chronological coefficient recurrence counts each column set once.
Upward-rounded positive arithmetic at 80 and 120 decimal digits gives

    P_full_columns < 2.349 * 10^(-88) < 2^(-291).

The h=33 term is largest. The full-support h=64 term is still covered more
tightly by the earlier tau=256 calculation. Reproduce this combined bound
with `theory/check_full_columns.py`; no new source-multiplicity assumption
or global-rank operator is used.

The [shared-outer pair moment](PARTIAL_COLUMN_MOMENTS.md) extends the minimum
h=33 class to every partial scalar mask. It uses tau=65536, z=999/1000,
and lambda=4981/5000. Its directed bound is

    P_33_all_masks < 1.084 * 10^(-66) < 2^(-219).

This averages the actual independent column mixers; it is not conditional
on every outer realization. Both bounds include every group choice and
allow a free terminal state. They overlap and do not cover h>33 partial
columns, multiple active groups, or the complete code.
