# Random state pairings: exact local law and inverse proof interface

This research branch changes only the inner pairing schedule. Its purpose
is to prevent a proper state support from remaining aligned with several
successive fixed matchings. The RS outer, regional routing, state width,
and accepted SL2 arithmetic remain the controls from
[UNIFORM_PAIR_PINS.md](UNIFORM_PAIR_PINS.md).

The [projection continuation](RANDOM_PAIRING_PROJECTION.md) now certifies
two constant-support families; whole-code closure and timing remain open.
In particular,
the previous eight-defect certificate belongs to the cyclic bit-pairing
construction and does not transfer. Uniform random matchings can preserve
a proper support for arbitrarily many steps with positive probability;
the old deterministic eight-step proper-run bound no longer applies.

## Setup and update

Use J=GF(2^32), K=2^20, N=2^21, state width n=512, and L=4096 updates.
At setup, independently sample a uniform perfect matching M_t of the 512
state coordinates for each t. Each matched pair receives an independent
accepted SL2 matrix A as before: both A and A+delta I are superregular.
All matchings, matrices, and the existing regional route are then fixed
for every message. The matchings are independent of the outer and route.

The block-diagonal operator in physical coordinate order is A_t. Keep

    Q_0=0,    Q_(t+1)=A_t Q_t+X_t,    Y_t=Q_(t+1)+delta Q_t.

Only one pair layer is applied per update, and every update emits a full
512-symbol batch. This adds no global input permutation, but arbitrary
local partners may cost state shuffling or indexed accesses. No claim
about that cost is made before measuring an actual kernel.

## Exact matching law

Fix an incoming support H of size w before exposing the current matching.
Let C count pairs with one endpoint in H and one outside H. For
0<=k<=min(w,n-w) with k congruent to w modulo2,

    Pr[C=k] = binom(w,k) binom(n-w,k) k!
              (w-k-1)!! (n-w-k-1)!! / (n-1)!!.

The convention (-1)!!=1 handles empty residual sets. Choose the k active
and k inactive crossing endpoints, biject them, then match each remaining
side internally. These cases partition all perfect matchings.

Let F be the coordinates in live pairs. Then |F|=w+C and H is contained
in F. Conditional on C=k, the set F minus H is a uniform k-subset of
the inactive coordinates. In particular,

    E[C]=w(n-w)/(n-1).

For even w, the probability of no crossing is

    Pr[C=0]=binom(n/2,w/2)/binom(n,w).

At n512,w256 this is approximately2^-255.5007, and E[|F|]=384.2505.
At w2 the no-crossing probability is1/511, approximately2^-8.9972.
These are local fixed-support probabilities, not simultaneous assertions
over supports chosen after seeing the matchings.

The mean describes live pair coordinates, not actual next-state weight.
Input cancellation and two-active-coordinate cancellation must still be
counted. For X=0, superregularity guarantees that a crossing pair produces
two active coordinates, but an active-active pair can lose one coordinate.

## Mark the inactive positions that matter

Fix a marked set B outside H, of size b. Let J_B=|(F minus H) intersect B|.
Conditional on C=k, its exact distribution is

    Pr[J_B=j | C=k]
      = binom(b,j) binom(n-w-b,k-j) / binom(n-w,k).

For example, suppose all coordinates in B must be newly active in Q_next
but have X=0. Every such coordinate must belong to F; an empty pair has
Q_next=X. The necessary event J_B=b can be charged with this joint law.
At w256 its probability is approximately2^-0.9972 for b1,2^-15.9535 for
b16, and2^-63.7927 for b64. It is not valid to charge that event when the
outer input is allowed to supply the births instead.

This is the first interface to the outer: retain which new-state positions
are forced to have zero input, instead of averaging the live mass alone.
The matching event and its field constraints must be counted jointly.

## Fix outputs to remove the setup-dependent-support ambiguity

Fix a complete candidate output y before sampling the inner setup.
Reconstruct Q_0=0 and Q_(t+1)=delta Q_t+y_t. Every state, including the
free terminal state Q_L, is now fixed independently of all matchings.
Thus the preceding local law applies separately at every t. This avoids
assuming that a bad message chosen after setup has an independent state.

Fix also a realization of the regional route. Let C_pi be the routed
outer code and let H_pi be a parity-check matrix for C_pi. Conditional on
all matchings, let F be the union of the live-pair coordinates across time
and D its complement. In the inverse experiment,

    X_t=Q_(t+1)+A_t Q_t.

On D, the input is fixed: X_D=(Q_next)_D. On F, the accepted-map law is
pointwise dominated by the uniform law on J^F with the retained single
global density factor less than2. This uses independent maps per pair,
not fresh randomness each time an encoder is called.

Consequently, for this y and these route/matching realizations,

    Pr_maps[X in C_pi]
      <= 2 * q^(-rank((H_pi)_F))
         * 1[(H_pi)_D (Q_next)_D belongs to image((H_pi)_F)].

The indicator tests whether the fixed complementary inputs are compatible
with the code. For a uniform vector on F the expression without2 is exact:
the affine parity equations have either no solutions or q^(|F|-rank) solutions.

A valid whole-code first-moment starting point is therefore

    Pr[exists nonzero output of weight<=W]
      <= 2 sum_(0<wt(y)<=W) E_(route,matchings)[q^(-rank((H_pi)_F)) * indicator].

The sum is over fixed output vectors, not messages selected after setup.
Injectivity of the outer and invertibility of the inner ensure every
nonzero codeword appears in this union. The density factor is paid once,
not at every temporal step after it has already been bounded globally.

This formula is an interface, not an evaluated certificate. Summing it
efficiently while preserving the indicator is the next central obligation.
The RS rank depends on coordinates across regions, not only on |F|.

## Why expected spreading alone does not close the proof

Even in the favorable case where all selected parity columns are independent,
the one-step inverse factor is q^-w E[q^-C], not q^(-w-E[C]).
At q=2^32,w256, the additional exponent -log2 E[q^-C] is only about255.5007
bits. The C=0 term contributes more than0.99999999999999 of that moment.
The typical128 extra coordinates would suggest roughly4096 additional
bits, which is not justified. Rare matchings dominate the inverse moment.

The same effect occurs at smaller weights: at w16 and64 the matching-only
gain is about50.8668 and138.6544 bits respectively. This does not refute
the construction. It shows that marginalizing the matching before retaining
input compatibility discards precisely the constraint needed to penalize
the rare low-spreading histories.

Do not classify all singleton pairs as a few exceptional batches in this
new model. A half-active state normally has about128 singleton pairs in
every update. Their number is a useful spreading statistic, not a rare
defect of the old kind. A new transfer count should retain crossing count,
next-support overlap, and outer-input zeros explicitly.

## Checks, preserved controls, and next step

`diagnose_random_pairing.py --test` exhausts11463 matching/support cases
through width10, checks160 marked-set profiles, and checks256 affine
parity-probability cases over GF2. Conditional external subsets are checked
uniform, not inferred from the crossing mean. All tests pass.

The default command computes exact rational laws at width512 and prints
floating-point logs only for interpretation. These logs are not directed
certificate margins. The test and diagnostic took about0.05s and0.15s,
ran serially under30s limits, and wrote no files. The pre-run disk guard
requires20GiB free. No kernel, paper, library, frozen proof, or benchmark
was changed. CWC was used for the proof interface and scope statements.

The [joint pin continuation](RANDOM_PAIRING_PINS.md) now derives a local
transfer polynomial and checks its composition with a tiny outer code.
It retains the old pin identity together with
the new matching law. Empty pairs fix inputs to Q_next; births on those
pairs pin new temporal components. Test the half-support persistence and
activation-to-reset examples first, with RS-supported input constraints.
Only then evaluate whether a precomputed arbitrary-partner kernel is worth
benchmarking. The goal remains active; no whole-code distance is claimed.
