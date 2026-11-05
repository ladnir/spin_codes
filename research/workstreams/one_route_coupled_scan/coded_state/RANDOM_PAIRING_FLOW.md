# Local capacity sharing and an exact pairing process

This continues the cut-scope diagnostic. The encoder is unchanged. The
objects below are auxiliary flows used to lower-bound coverage plus deaths;
their five-state variables are not the encoder's field-valued state.

## Why share capacity?

Fix a target set T of 1974 outer rows. Each pair gate has capacity equal
to the number of its two routed positions belonging to T. Equal splitting
can waste capacity on a coordinate timeline that does not need it.
Instead, assign that capacity between the two incoming timelines.

The standalone helper `random_pairing_demand_flow.py` computes exact
timeline cuts and optimizes one gate allocation at a time. Its selected
32-update window, seed 20261028, gives these lower bounds:

| Target-label model | Equal split | Demand-dependent split |
|---|---:|---:|
| IID labels | 1303 | 1513.5 |
| Exactly 1974 labels per 8192-position region | 1280.75 | 1483.75 |

These are sampled-window bounds, not full-geometry bounds or setup tails.
The allocation search is coordinate ascent, not a global optimum claim.

## A causal rule with five auxiliary states

Scale capacities by four. A timeline has unit supply and temporal edge
capacity two. Its cut recurrence needs only the clipped difference
c in {-2,-1,0,1,2} between its two conditional cut costs.
Start every timeline at c=0.

For incoming auxiliary states a,b and a gate with j target positions,
set k=4j and allocate

    l = clamp(floor((k+b-a)/2), 0, k),  r = k-l.
    a' = clamp(a+l-1, -2, 2),  b' = clamp(b+r-1, -2, 2).

Write neg(c)=max(-c,0). The increment in the sum of timeline cut values is

    gain = 2 + neg(a)+neg(b) - neg(a')-neg(b') - loss,
    loss = 1[j=0] (1[a=-2]+1[b=-2]).

Thus the total flow is the number of coordinate-updates minus accumulated
loss and the final sum of neg(c). This is an exact identity.
It turns the flow certificate into a tail question about overflow counts.

The local rule is also monotone and nonexpansive. Increasing one incoming
state by one increases both outgoing states weakly, and

    outgoing L1 increase + loss reduction <= 1.

Increasing j by one gives the analogous upper bound four. Telescoping
over a common future schedule bounds the cumulative loss change caused
by a one-unit initial perturbation by one, independently of the horizon.
These statements follow by exhausting the finite local domain. They do
not by themselves give a sufficiently strong setup tail.

For the R=256, half-full-history diagnostic, the full sequence has 4095
random service updates followed by one zero-service terminal update.
The first input batch supplies no service to a preceding timeline.
The final loss plus final negative potential is at most 3*512.
To certify C_T+D >= 168785, it therefore suffices that random-update loss
is at most 371900. Failure of this sufficient flow certificate requires
loss at least 371901, about 17.738% of the 4095*512 coordinate-updates.
This is not an equivalence with a bad codeword or a distance failure.

## Re-pairing must remain in the model

An independent-partner mean-field calculation predicts approximately
12.34% loss and a full-sequence dual lower bound 197313. This is only
a diagnostic approximation, not the finite-population distribution.

As a control, keep partners fixed forever and use IID target labels with
p=1974/8192. Equal initial states within a pair stay equal. The resulting
five-state Markov chain has an exactly checked rational stationary law,
with loss approximately 15.827%. This is not a counterexample to uniform
re-pairing. It shows why a proof that forgets re-pairing can lose useful slack.

## Exact finite-population transition

Uniform re-pairing permits an exact summary by the histogram h of the
five auxiliary states. Let n=sum h_i be even. For each unordered state
pair i<=j, let e_ij count matched pairs of that type. Feasibility requires
2e_ii + sum_{j!=i} e_min(i,j),max(i,j) = h_i. Then

    Pr[e | h] = (product_i h_i!) /
      ((n-1)!! product_i (2^e_ii e_ii!) product_{i<j} e_ij!).

This follows by partitioning each color class among its partner colors,
matching cross-color vertices bijectively, and pairing same-color vertices.
Divide the resulting count by the total (n-1)!! labeled perfect matchings.

Conditional on e, assign each pair j=0,1,2 target positions, with
multiplicity binomial(2,j), then apply the causal rule above. Retain jointly
the new histogram, accumulated loss, and total number J of target labels.
If a region has m positions and r target labels remaining, each specified
label vector for the next n positions has probability

    binomial(m-n, r-J) / binomial(m,r).

Infeasible binomial coefficients mean probability zero. Update the remaining
counts by (m,r) -> (m-n,r-J). Regional boundaries reset these counts to
(8192,1974). Account separately for the ignored initial input batch.
Thus the histogram together with regional remaining counts is an exact
Markov summary; independence of partner states is not assumed.

The formula applies at n=512. Direct enumeration is deliberately capped
at n=8 in `random_pairing_flow_kernel.py`. Exact rational checks compare
36 matching profiles with recursive matching counts and eight complete
histogram/loss/label transitions with labeled brute force. All pass.
The same run passes 2916 threshold checks, 1215 allocation checks,
81 exhaustive timeline cuts, and 170 local stability checks.

## Scope and next test

This closes the local probabilistic interface, not the global tail.
Next derive a rigorous exponential-moment envelope for this actual
histogram process, retaining regional counts. Test the 1974-row class
first and account for the union over target sets. A mean-field average
or a fixed-target sampled flow cannot replace that step.
Even a successful structural test still needs the field-cancellation
argument and coverage of other state-weight profiles.

Validation ran serially and file-free on Peach because local free space
was below the 20 GiB numerical guard. The final exact suite took 0.293 s,
under a 30 s timeout and 256 MiB memory cap. No native encoding benchmark,
construction change, raw output file, paper edit, or promotion was made.

## The deterministic interface also covers varying state weights

The flow argument does not require every state to have weight256. Its
constant-weight presentation hid a telescoping endpoint term. Retaining
that term gives a simultaneous bound for arbitrary support histories,
including births, deaths, and complete resets. This extends the geometric
interface only; it is not a new distance certificate.

Fix the route, all pairings, a target row set T, and the resulting causal
service allocation. Write T_t for the target coordinates in update t.
Let H_1,...,H_L be any state supports, with H_0 empty. Define

    w_t=|H_t|,   S=sum_(t=1)^L w_t,
    D=sum_(t=1)^(L-1) |H_t minus H_(t+1)|,
    V=sum_(t=1)^(L-1) |H_t symmetric_difference H_(t+1)|.

Let F_t be the union of pairs meeting H_(t-1) under update t's matching.
The target coverage is C_T=sum_t |T_t intersect (F_t union H_t)|.
Here D counts no initial deaths because H_0 is empty; the terminal state
is free. Summing births minus deaths gives the exact identity

    2D = V + w_1 - w_L.

The allocation assigns service s_(t,i) to coordinate i of H_t. For t<L,
each pair {a,b} of update t+1 satisfies
s_(t,a)+s_(t,b)=4|T_(t+1) intersect {a,b}|. At t=L service is zero.
For any H, its charged service is at most4C_T: a pair with no active
incoming endpoint charges zero, and every live pair charges at most its
full target capacity. Outgoing-state coverage is only additional coverage.

Let F be the sum of the exact timeline minimum-cut values computed with
unit supply and temporal capacity two. Evaluate each timeline cut on the
bits indicating membership in H_t. The cut value gives

    F <= nL - S + sum_(t,i in H_t) s_(t,i) + 2V
      <= nL - S + 4C_T + 2V.

Consequently, for every support history simultaneously,

    C_T + D >= (F-nL+S+2(w_1-w_L))/4.

The causal loss identity remains unchanged. If L_random is the loss in
the first L-1 service updates, the final zero-service loss plus final
negative auxiliary potential is at most3n. Therefore

    C_T + D >= (S+2(w_1-w_L)-L_random-3n)/4.

This bound is deterministic after setup and allocation. The allocation
does not depend on the chosen code-state history. A setup event bounding
L_random for all target sets of a given size therefore controls every
history satisfying the stated activity condition; it needs no extra union
over those histories for this geometric conclusion.

### Consequences for the selected structural threshold

At n512 and L4096, put A0=nL/2=1,048,576. The existing loss condition
L_random<=371900 implies C_T+D>=168785 for every history with

    S+2(w_1-w_L) >= A0.

This includes histories of varying weights whose first and last weights
agree and whose total activity is at least A0. It also includes every
history with S>=A0+1024, regardless of endpoint weights.

If only S>=A0 is known, the worst endpoint penalty is1024 in the numerator.
The sufficient loss condition becomes L_random<=370876. The corresponding
failure threshold is370877 instead of371901. At theta=1/10, this reduces
the per-update logarithmic allowance by1024/(10*4095), without changing
the auxiliary transition kernel or potential. This is an optional stronger
tail requirement, not an established global tail bound.

An80-digit directed calculation gives allowance greater than
7.654289121736758 for that stronger requirement. The current checkpoint's
worst upper bound,7.651481649163582 rounded upward, is still below it.
This comparison applies only to the already checked auxiliary histogram
regions. The unverified complement prevents either global tail claim.

The helper's `--tests-only` mode exhausts all4096 histories on each of36
width-four, three-update geometries. It checks737354 identities and
inequalities, including135000 nonconstant-weight histories. The cases
include empty/full tags, asymmetric labels, three different matchings,
births, and resets. It independently evaluates coverage, allocated service,
timeline cuts, and the causal loss identity. The original4382 local tests
also pass. The analytical inequalities above establish arbitrary-width
validity; these finite tests check the implementation, not that generality.

Remaining obligations are substantial: complete the setup tail over all
auxiliary histograms and relevant target counts, then combine coverage
with actual field cancellation and output-weight accounting. In particular,
the endpoint identity alone does not certify histories below the activity
threshold or prove whole-code minimum distance.

    python -B research/workstreams/one_route_coupled_scan/coded_state/random_pairing_demand_flow.py --tests-only
