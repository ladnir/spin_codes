# An obstruction to the unweighted one-variable comparison

Improving the current global envelope alone cannot close the second-moment
argument. This note proves that limitation for the comparison in
`CONVEX_OVERLAP_REDUCTION.md`, using the sharper law nu_J from
`JOHNSON_OVERLAP_LAW.md`. It does not lower-bound actual setup failure.

## A term that every such comparison must contain

For each possible total core overlap A, allow a separately chosen positive
log-convex envelope h_A. Require the previous uniform region bound:

    R(x,y,k+d) <= h_A(k)

for every covered even marginal pair x,y in [740,900] and correction
increment d in {0,1,2,3}. In particular,

    h_A(740) >= R(740,740,740) = 1/beta_740.

The equality holds because identical region inputs pass the same kernel
test with probability beta_740, rather than beta_740 squared.

Consider identical core supports, so their intersection has size 2610.
This support event has probability 1/binom(8189,2610). In the artificial
independent nu_J row law, retain just the following composition:

    1592 rows of overlap80,
    1000 rows of overlap61,
      18 rows of overlap60.

Its total overlap is A=189440=256*740. Its probability is

    (2610!/(1592!1000!18!))
      * nu_J(80)^1592 nu_J(61)^1000 nu_J(60)^18.

All three masses are positive in the exact certificate. In the
multinomial comparison conditional on A, retain just the count vector
whose every coordinate equals 740. Its probability is

    189440! / ((740!)^256 * 256^189440).

The envelope product at this vector is at least beta_740^(-256).
Multiplying these four factors proves that the comparison is larger
than 2^23000. The stored log2 lower bound is approximately 23189.64854.

`certify_convex_route_obstruction.py` evaluates the bound with outward
256-bit arithmetic; its 512-bit replay passed. The proof applies to
different h_A for different A, not merely to the original shared h.
It also applies when the probability of A is upper-bounded by Chernoff
tilts of the nu_J generating function: every such bound is at least
the actual coefficient of that artificial generating function.

The required normalized second-moment upper bound is below (9/16)*2^40.
Thus this comparison cannot provide the desired failure-probability result,
even with exact coefficient evaluation and ideal envelope optimization.

## What was lost, and what the obstruction does not say

The argument above combines two relaxations: a comparison distribution
for row overlaps and a separate worst case over each region's marginals.
Those choices need not describe one actual message pair.

For example, the worst-case diagonal type (740,740,740) repeated in all
regions has total marginal weight 189440 in each message. Every actual
core message has weight 208800 before the correction rows are added.
Repeating that worst-case type therefore violates the fixed marginal totals.
The unweighted envelope is nevertheless forced to pay its cost in each
region separately. The artificial comparison law also saturates some
rare-overlap bounds that may substantially exceed actual probabilities.

Neither the composition above nor the multinomial event establishes an
actual family of bad setups. The fixed SPIN target remains unresolved.

## A repair that retains the fixed marginal totals

The following weighted region bound avoids the particular forced term.
Set mu=208800/256=815.625. For each total overlap A, choose b_A<=0
and a positive log-convex h_A satisfying

    R(x,y,k+d) <= exp(b_A*(x+y-2mu)) h_A(k)

for all retained region types. Every completed candidate has total input
weight at least 208800, since correction rows only add bits. Consequently,

    b_A * sum_j (x_j+y_j-2mu) <= 0,

and multiplying the weighted bounds still gives

    product_j R(x_j,y_j,k_j+d_j) <= product_j h_A(k_j).

The remaining subset comparison and actual-generating-function bounds
from the preceding notes therefore remain applicable. The factor involving
marginal weights introduces no positive global cost.

At the diagonal type (740,740,740), the required lower value of log h_A(740)
changes from -log beta_740 to

    -log beta_740 + 151.25*b_A.

This decreases when b_A<0. Other region types can become more expensive,
so b_A must be optimized jointly with the envelope, separately for each A.
The argument does not claim that such optimization will succeed.

Next, screen this weighted bound while retaining every overlap and the
exact support-intersection law. If it remains too loose, stronger joint
row bounds or a comparison that keeps both marginal count vectors is needed.
