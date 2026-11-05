# Replace support histories by covered input sets

2026-10-30. Construction unchanged; the full45%/40-bit goal remains open.
The [filtered-count obstruction](QUARTET_FILTERED_COUNT_OBSTRUCTION.md)
shows that local history filters leave too much multiplicity. This note
eliminates support histories exactly before bounding the geometric event.
It gives a new deterministic formulation and checked fixed-geometry bounds,
but not the missing ensemble tail.

Continuation: [the quota audit](QUARTET_WINDOW_QUOTA.md) proves that even
an exact Bernoulli moment cannot close the middle case within the multiplier
family below. The independent-column diagnostics miss an impossible dense
stream introduced by the comparison. The replacement retains the unpaid
count exactly and has a checked five-type temporal transfer. The earlier
fixed-geometry allocation certificates remain valid.

## Optimize the support after fixing which inputs may be covered

Fix the route, quartet partitions, and an inactive-row set T of size2100.
Let I be its537600 inactive input positions. A state coordinate touches its
own input and all four inputs of its quartet at the next update. At the
last update it touches only its own input.

For P contained in I, permit covering only positions in P. A state coordinate
is then allowed precisely when none of the inactive inputs it touches lies
in I outside P. Write z_(t,j)(P) for this allowed-coordinate indicator.
For lambda,tau>0, define

    H_(lambda,tau)(P)
      = max_(S uses only allowed coordinates) [lambda*A(S)-tau*D(S)].

The initial support is empty; terminal deaths are not charged. Once P is
fixed, optimization separates into maximal allowed runs in each coordinate.
A nonterminal run of length l contributes max(0,lambda*l-tau); a terminal
run contributes lambda*l.

For lambda=1/4,tau=37/50, let W_l(P) count all windows of l consecutive allowed
coordinates, and E_l(P) count columns whose last l coordinates are allowed.
The exact identity is

    H(P)= W_3(P)/100 +24 W_4(P)/100
          +25 E_1(P)/100+25 E_2(P)/100+24 E_3(P)/100.       (1)

The last three terms handle the free terminal state exactly. A nonterminal
window of length l depends on4l+1 distinct input positions: its first own
input and l subsequent quartets. A terminal window depends on4l-3 positions.
Thus the windows in (1) involve at most17 input positions.

For each window w, let B_w be the inactive inputs among those positions,
and let a_w be its weight in (1). The window is enabled precisely when
B_w is contained in P. We therefore have the exact optimization identity

    max_S [A(S)/4-C(S)-(37/50)D(S)]
      = max_(P subset I) [sum_w a_w 1[B_w subset P]-|P|]. (2)

For the forward inequality choose P to be S's actual covered set. For the
reverse inequality choose an allowed support attaining H(P). Its actual
covered set is contained in P, so charging |P| only lowers its score.
No independence assumption is used in (1) or (2).

## Every allocation gives a checked upper bound

For each window with nonempty B_w, distribute its weight among those labels:

    x_(w,i)>=0,  sum_(i in B_w) x_(w,i)=a_w.

Let L_i=sum_w x_(w,i), and let H_empty=sum_(B_w empty) a_w. For every P,

    sum_w a_w 1[B_w subset P]-|P|
      <= H_empty+sum_(i in P)(L_i-1)
      <= H_empty+sum_i max(0,L_i-1).                      (3)

This is a deterministic certificate for all supports on the fixed geometry.
It does not require that the allocation be optimal. The window weights
total N/4 exactly, including the boundary terms. A corresponding bipartite
flow has one source of capacity a_w per window, unit capacity per inactive
input, and edges to its blockers. Max-flow/min-cut gives (2); an arbitrary
allocation gives the explicit upper certificate (3).

The helper uses integer units of1/100. It initially processes low-degree
windows first, assigning each unit to a least-loaded blocker. Subsequent
passes remove a window's old assignment and redistribute its mass. Every
reported bound independently reconstructs all loads and checks conservation.
The update order is only a diagnostic algorithm, not an assumed random law.

On the SAME seed20261027 geometry and fixed row set as the earlier min-cut
test, H_empty=3816.84 and the exact optimum is3988.44. Results are:

| Allocation | Deterministic upper bound | Excess above H_empty |
|---|---:|---:|
| Equal fractional shares | 46354.2 | about42537.3 |
| Initial least-load pass | 8545.79 | 4728.95 |
| One reassignment pass | 4099.26 | 282.42 |
| Two reassignment passes | 4072.10 | 255.26 |
| Three reassignment passes | 4059.83 | 242.99 |

The final bound is only71.39 above the exact optimum. It is NOT a tail
bound over setups or row sets. Sequential assignments can depend on distant
windows; their dependence radius is not bounded merely because B_w is small.
Do not apply a local read-degree inequality to this algorithm's output.

`quartet_window_allocation.cpp` checks524288 exhaustive toy histories and128
complete paid-label sets, including equality in (2) and both boundaries.
Optimized and UBSan checks pass. The full-size optimized diagnostic takes
about3.24 seconds, uses fewer than256MiB of owned arrays, and writes no files.
These are proof-diagnostic runtimes, not encoder timings.

## Keep the field constraints on the actual covered set

This representation preserves c=|P| when P is an actual codeword's covered
set. It does not replace P by the set of a geometric maximizer. On the
previously certified good events for this row class,

    A-C>=713983,  A>=743225,  C+D<=255463,
    A<=604080+5C.

Thus every remaining codeword has27476<=c<=255463. The directed joint-prior
coefficients also imply the following slightly rounded necessary condition:

    (39/512)A-C-(381/512)D+111127 >0.                      (4)

The constant includes the maximum contribution from the free terminal weight.
The helper verifies these rounding directions against `selected_prior_weights`;
the displayed decimals in older notes are not used as certificate inputs.

For arbitrary mu,nu,eta>=0, add nonnegative multiples of the three field
constraints to (4). Set

    lambda=39/512+mu+nu, tau=381/512+eta,
    B(c)=(1+mu+eta)c+713983mu+743225nu-255463eta-111127.

Every remaining low-output history then requires

    H_(lambda,tau)(P)>B(c), |P|=c.                        (5)

This keeps the actual activity, rank, and turnover restrictions BEFORE
maximizing the allowed runs. It avoids imposing codeword restrictions on
an unrelated canonical geometric witness. The multipliers are proof
parameters; they do not modify the encoder.

## A quota-preserving count, and its remaining loss

Let R=537600, p=1-2100/8192. For a fixed original row set, choose each member
of I independently into P with probability t/(1+t), where t>0. Then

    sum_(P subset I) t^|P| b^H(P)
      =(1+t)^R E_P[b^H(P)], b>=1.                        (6)

For a fixed c, upper-bound its coefficient by t^-c times (6). The newly
allowed input bits have marginal probability p'=1-(1-p)/(1+t). They are
not independent under the actual quotas. Nevertheless the positive-window
polynomial comparison applies: generate independent clearing bits at EVERY
input, and OR them with the original active bits. Expand in the original
active indicators with nonnegative coefficients. Each monomial is bounded
by its independent-bit counterpart, region by region. This proves the
comparison without a conditioning penalty.

Let M(p',b^4;lambda,tau) be the one-column moment from the generalized run
recurrence, excluding the terminal correction bounded by512tau. The read4
inequality therefore gives the valid fixed-c tail bound with log2 at most

    -log2(b)*(B(c)-512tau) + R log2(1+t)-c log2(t)
      +128 log2 M(p',b^4;lambda,tau)+log2 binom(8192,2100). (7)

A union over c is still required. Three bounded evaluations of (7) remain
positive, even after the field constraints. The initial unweighted version
also failed. No full count range is certified by these tests.

To diagnose the loss, replace read4 by independent whole columns ONLY as an
explicitly false idealization. The following rounded log2 values compare
selected expressions, with their stated moment parameters:

| c | (mu,nu,eta) | t | Valid read4, log2 b=1 | Hypothetical independent columns |
|---:|---|---|---:|---:|
| 27476 | (0,1/4,0) | 6869/127531 | +109939 | -49872 at log2 b=4 |
| 111572 | (1/5,0,0) | 1/6 | +287109 | -3437 at log2 b=3 |
| 255463 | (0,0,1) | 255463/282137 | +399923 | -30406 at log2 b=4 |

The negative numbers are NOT bounds for the construction. They show why
retaining actual quartet mixing could be useful, rather than assuming
independence or spending more time evaluating the same read4 relaxation.
For the middle case only three deterministic t choices were checked.
The strongest diagnostic still has limited room for a dependence correction.

## Exact spatial dependence for one window start

We can measure one source of that correction without assuming independent
columns. Work under independent active-input bits of probability p, after
the justified comparison above. For a specified k-subset of512 coordinates,
let G_k(p) be the probability that all quartets meeting it are active. The
quartet partition is uniform. Its exact expression is

    G_k(p)=[x^k](1+p^4*((1+x)^4-1))^128 / binom(512,k).   (8)

Consider one nonterminal start time with window rewards alpha W3 and
(lambda-alpha)W4, where alpha=3lambda-tau is nonnegative. Put

    a=b^alpha-1, d=b^lambda-b^alpha.

Expanding the product over columns gives the exact one-start moment

    sum_(k=0)^512 binom(512,k) p^k G_k(p)^3
      * sum_(h=0)^k binom(k,h) a^(k-h) d^h G_h(p).        (9)

The first three future partitions must cover the union of the marked
three- and four-window columns; the fourth partition only covers the latter.
Those partitions are independent. Neighboring start times SHARE partitions,
so their moments in (9) cannot be multiplied.

Integer coefficient extraction computes every G_k exactly. The final
logarithmic evaluation of (9) is currently FLOAT diagnostic arithmetic:

| Parameters | Actual one-start log2 moment | Independent columns | Read4 bound |
|---|---:|---:|---:|
| p=1523/2048, lambda=1/4,tau=.74,log2 b=1 | 1.00215 | .98443 | 1.27362 |
| p=799/1024,lambda=707/2560,tau=381/512,log2 b=3 | 13.89185 | 11.88818 | 26.96341 |

The second row has about2.004 bits of spatial excess at ONE start. It is
not negligible, and cannot be extrapolated as an additive whole-trajectory
correction. Formula (8) passes315 exact small-partition checks; (9) passes
256 exact single-quartet category checks.

## Reproduction and next discriminating task

The new C++ helper includes the existing geometry diagnostic, without changing
the geometry distribution or any production code. `quartet_run_components.py`
adds `--coverage-sets` and `--field-coverage`; the latter imports the existing
directed joint-prior coefficients. The thinning test checks448 paid subsets
and256 product-measure masks. `quartet_window_spatial.py` reproduces (8)-(9).
All numerical work is serial, remote, watchdog-bounded, and file-free except
temporary C++ builds, which are removed. Prior frozen helpers are unchanged.

Next retain short temporal overlap together with the actual partition law
in (8), or prove a tail for a suitably analyzable allocation in (3). First
test the middle covered-count class111572: it has the narrowest diagnostic
room. Do not multiply one-start moments or assume that local reassignment
passes have local dependence. A successful bound must then cover the entire
covered-count interval and the other unresolved row counts. No new allowance
has been added to the whole-code failure budget.
