# Retain a short empty tail before an active event

A paired input can reset a sparsely supported state with fewer coefficient
charges than a full state. The current two-state envelope maximizes that
reset cost before processing the preceding empty gap. It can therefore
forget the cost of reaching the required sparse endpoint.

This note retains a bounded number of empty steps immediately before each
event. It changes only the proof envelope. The final section records local
diagnostics and selected two-group bounds, not a whole-code certificate.

## 1. Keep functions of the incoming state until the last possible step

Use the ordinary, non-lookahead continuation from
[EMPTY_SPAN_RANK.md](EMPTY_SPAN_RANK.md), with the separated rank price in
[RANK_PRICE_SEPARATION.md](RANK_PRICE_SEPARATION.md). Fix tau>=1, 0<z<=1,
q=256, and the width-four antipodal wiring. Let T be the empty-step Bellman
operator on valid live states (r,H), where 1<=r<=|H|. Extend T to the quiet
state by T(f)(quiet)=f(quiet).

Fix a positive live potential v and lambda>0 whose inequality

    T(v) <= lambda v

has been rigorously checked. Here v is zero on the quiet state. The two
terminal basis functions are

    phi_q = 1[quiet],
    phi_l(quiet)=0,       phi_l(r,H)=v(r,H).

The ordinary empty-step matrix on their nonnegative coefficients is
E_env=diag(1,lambda). All potentials continue to use the same global
tau^(-dim W); neither r nor the gap boundary replaces W.

Fix a nonempty active input mask I, and put k=|I|. Before maximizing over
the incoming state, retain these two terminal functions for the event:

    f_(I,q)(quiet)=0,
    f_(I,q)(r,H)=M_I(z;psi)[H,empty],

    f_(I,l)(quiet)=z^k max_(1<=s<=k) v(s,I),
    f_(I,l)(r,H)=sum_(H' nonempty) M_I(z;psi)[H,H']
                     max_(1<=s<=min(r+k,|H'|)) v(s,H'),

where psi=min(tau,q). The ordinary matrix M_I already includes the output
weight and coefficient charges. These functions bound the event followed
by phi_q and phi_l, respectively, for every retained incoming quotient and
every permitted collection of actual input forms.

In particular, f_(I,q) retains which supports can reset and their costs.
Replacing it immediately by max(f_(I,q)/v) discards that information.
Any stronger uniformly proved event functions may be substituted. If they
use a freshness hypothesis, that hypothesis must remain valid throughout
the preceding empty gap.

## 2. A linear envelope for a fixed preceding gap

For an integer L>=0 and j in {q,l}, define

    c_(I,j,L)=max_live T^L(f_(I,j))/v.

Then the matrix

    C_L(I) = [[f_(I,q)(quiet), f_(I,l)(quiet)],
              [c_(I,q,L),     c_(I,l,L)]]

bounds L empty steps followed by the active event I. Its row and column
conventions are the same as E_env: rows describe incoming quiet/live
states; columns describe the terminal basis functions.

To prove the bound, first apply the uniform active-event inequality for
each terminal basis function. Apply the empty-gap Bellman bound to that
function, and only then divide by v and maximize over live states. Quiet
states stay quiet during the gap, so the top row is unchanged.

For an arbitrary terminal continuation c_q phi_q+c_l phi_l, with nonnegative
coefficients, sum the two bounds with these coefficients. Actual weighted
expectations are linear in the terminal function. This proves a linear
matrix envelope even though T itself is nonlinear. The same coefficients
and the same potential v are used at every boundary.

For L=0, C_0(I) is the earlier active envelope. For L>0, it can be strictly
tighter than E_env^L C_0(I). No such improvement is assumed without checking
the resulting constants.

## 3. Compute only a bounded tail

Choose a small integer ell>=0. Compute T^s(f_(I,j)) for 0<=s<=ell, using
exact arithmetic or verified upper bounds. For L>=ell, monotonicity and
positive homogeneity give

    T^L(f_(I,j)) <= c_(I,j,ell) lambda^(L-ell) v

on live states. Hence a valid fixed-gap envelope is

    C_hat_L(I) = C_L(I)                         if L<ell,
                 E_env^(L-ell) C_ell(I)         if L>=ell.

This accounts for exactly L empty steps and one active step. In particular,
do not use the long-gap formula with a negative exponent when L<ell.

For the reset column, compare

    c_(I,q,ell)/lambda^ell

with the old constant c_(I,q,0). A reduction quantifies the endpoint gain
for long gaps. The short-gap terms still have to be included in the route
average. A better long-gap constant alone is not a regional certificate.

For a suffix of L empty steps with no following event, the original
E_env^L remains valid. Empty regions likewise use E_env^m. The terminal
state remains free, with coefficient vector (1,1/min_live v)^T.

## 4. Average the actual one- and two-event route after bounding each gap

Each region contains m=2048 width-four batches. Define averaged fixed-gap
envelopes for the actual input masks:

    Cbar_L(P) = (C_hat_L(3)+C_hat_L(12))/2,

    Cbar_L(C) = (C_hat_L(5)+C_hat_L(6)
                 +C_hat_L(9)+C_hat_L(10))/4,

    Cbar_L(1) = (C_hat_L(1)+C_hat_L(2)
                 +C_hat_L(4)+C_hat_L(8))/4.

Numeric masks are four-bit coordinate masks. P fills one physical pair,
C is a collision of two singleton source pairs, and 1 is a singleton.
Compute each mask's tail bound before taking these averages. No nonlinear
Bellman operator is interchanged with a route expectation.

For one uniformly placed event of type A in {P,C,1}, set

    B_tail(A) = (1/m) sum_(L=0)^(m-1)
                         Cbar_L(A) E_env^(m-1-L).

This retains zero-length prefixes and suffixes. It contains m-1 empty
steps and one event, not m empty steps followed by an event.

For two singleton events in distinct batches, put

    B_tail(1,1) = 1/binomial(m,2)
      sum_(L0,L1,L2>=0; L0+L1+L2=m-2)
          Cbar_(L0)(1) Cbar_(L1)(1) E_env^L2.

The three gaps identify the two occupied batch positions uniquely.
Conditional on distinct batches, their positions are uniform among the
binomial(m,2) choices. Their singleton coordinates are independent uniform
choices, so the two separately averaged Cbar matrices are valid. The
matrix product remains chronological and need not commute.

For two adjacent active groups, the regional envelopes are therefore

    B_P = B_tail(P),

    B_S = p_col B_tail(C) + (1-p_col) B_tail(1,1),
    p_col = 1/(2m-1).

These are the same exact pair-slot probabilities as in
[TWO_GROUP_ROUTE.md](TWO_GROUP_ROUTE.md). Pairing regions contain one full
pair; the other regions contain two singleton pairs sampled without
replacement. The collision term is not discarded or replaced by another
independent singleton insertion.

The source forms may be correlated. Each fixed-mask envelope is uniform
over those forms and every incoming quotient, so their ordering requires
no extra independence assumption. Route averaging concerns only the
specified positions and swaps. Regions can then be composed in their
actual parity order, with the same initial quiet state and free terminal
vector as before. State is not reset at a region or column boundary.

## 5. Positive convolution without enumerating every gap pair

Fix an event type and write K_L=Cbar_L for its gap envelopes. Introduce
the matrix series

    K(x)=sum_(L>=0) K_L x^L,
    E(x)=sum_(L>=0) E_env^L x^L.

The one-event numerator is [x^(m-1)] K(x)E(x). The two-singleton numerator
is [x^(m-2)] K(x)K(x)E(x), with K from the singleton mask average.
The series products preserve the chronological matrix order from Section 4.

For any nonnegative matrix series J(x)=sum_(n>=0) J_n x^n, compute
the coefficients of K(x)J(x) as follows. Set S_n=0 for n<ell. For n>=ell,
define

    S_n=K_ell J_(n-ell)+E_env S_(n-1).

Then

    [x^n] K(x)J(x)
      = S_n + sum_(L=0)^min(ell-1,n) K_L J_(n-L).

Indeed, S_n is exactly the sum of E_env^(L-ell) K_ell J_(n-L)
over ell<=L<=n. Its recurrence separates the term L=ell and shifts all
remaining terms by one. The multiplication by E_env is on the left;
commuting it past K_ell would generally be wrong. The formula also permits
ell=0, with an empty short-gap sum and S_(-1)=0.

Apply this recurrence once to E(x), then again for a second event.
Each application uses O(m(ell+1)) operations on two-by-two matrices and
only nonnegative arithmetic. Directed upper rounding therefore preserves
the bound. No closed-form subtraction or quadratic gap enumeration is
required.

## 6. Recorded diagnostics and remaining coverage

The root agent ran [check_endpoint_tail.py](theory/check_endpoint_tail.py)
with tau=65536, z=.999, and lambda=.9962. It computed the bounded tails
with exact rational arithmetic. The following displayed ratios are rounded
diagnostics of c_(I,j,L)/(c_(I,j,0) lambda^L), not additional proof constants.

| Input mask I | Terminal column | Tail L | Normalized ratio |
| --- | --- | ---: | ---: |
| 3: physical pair | reset | 2 | about 4.568e-5 |
| 3: physical pair | live | 4 | about .258 |
| 5: cross-pair collision | reset | 4 | about 1.788e-7 |
| 1: singleton | reset | 4 | about 1 |
| 1: singleton | live | 4 | about .508 |

These results identify a large endpoint loss for pair resets. They do not
assume that a reset starts from a typical state: the maximization still
ranges over every permitted incoming state and common quotient.

The root agent then ran [check_endpoint_route.py](theory/check_endpoint_route.py)
with ell=4, retaining both ordinary events and the proved fresh-input
prefix from [TWO_GROUP_ROUTE.md](TWO_GROUP_ROUTE.md). The calculation uses
the positive convolution above, actual batch collisions, both adjacent-pair
parities, and a free terminal state.

Its 80- and 120-digit upward-rounded evaluations give about 5.423407141617e-86
for the h=33 shared-complete-column class. Summing the h=33,34,35 bounds
gives less than 5.510e-36. Both precisions explicitly verify that this
selected union is below 2^-117. It includes every shared column set and
all 8192 cyclic adjacent group pairs, with the terminal state free.
The union over all h=33,...,64 remains uninformative, about 6.992819e689.
These values cover the adjacent-group classes defined in that source note,
not arbitrary two-group supports. No independent execution was performed
for this note. Independent read-only audits found no discrepancy in the
endpoint lemma, fresh-prefix hypothesis, positive convolution, route
normalizers, chronology, or directed arithmetic.

The bounded exact tests compare six small route fixtures against 72
full-pair and 784 ordered singleton-pair slot/swap samples. They also
check 24 ordered insertion coefficients, 72 endpoint convolutions, and
24 endpoint regional averages against exact rational enumeration. Eight
zero-tail fixtures recover the original regional formula. These test
the arithmetic and route accounting; the analytic lemmas justify the
state envelopes for the actual field.

The next proof step is to examine the remaining large-h loss. The endpoint
refinement alone does not close all shared-column supports, partial columns,
nonadjacent groups, or larger active-group sets. Failed upper bounds are
not low-distance counterexamples.
