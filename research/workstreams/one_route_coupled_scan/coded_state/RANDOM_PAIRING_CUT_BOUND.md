# A verified cut bound for one selected row set

The [rank search](RANDOM_PAIRING_RANK_SEARCH.md) did not find an obstruction,
but local-search convergence could not rule one out. This continuation
bounds all half-full support histories for one selected target-row set
in the sampled geometry. It does not cover all row sets or prove distance.
The construction and accepted pair arithmetic are unchanged.

## What is fixed and what is bounded

Recreate the uniform regional route and independent uniform state
matchings from seed20261027. For the constant first-256-coordinate state
mask, compute each outer row's coverage. Fix T to the1974 least-covered
rows, breaking ties by row index. T therefore depends on the sampled
geometry, but is fixed throughout the following deterministic bound.

For arbitrary half-full masks H_t, let C_T count T's coordinates touched
by a live incoming pair or the next state. Let D count deaths. At the
46% target and R256, the contribution of T to the structural gap is

    g_T = 128|T|-C_T-D-83887.

The full parity deficit maximizes over T; bounding this one g_T does not
bound that maximum. For a message whose exact inactive row set is T,
the same expression also appears in the conditional field-count bound.
No resulting field-failure certificate is claimed here.

## Directed-cut representation

Partition4096 updates into disjoint32-update windows. Introduce one binary
variable for each next-state coordinate in a window. A target coordinate
after the first update has coverage cost equal to the OR of its outgoing
variable and the two incoming pair variables. Represent that OR by an
auxiliary node with unit cost and implications from those variables.
A death has cost x_prev*(1-x_next), represented by a directed unit edge.

At the first update of each window, keep only the outgoing-variable
coverage cost and drop the unknown incoming terms. Also drop deaths
across window boundaries. These omissions lower the cost. Thus summing
window lower bounds gives a lower bound on the full C_T+D.

Let n_var be the number of state variables in a window. Every half-full
history has n_var/2 active variables. For any nonnegative price mu,

    min_x(cost(x)-mu*|x|) + mu*n_var/2

is a lower bound on its cost. The graph adds source edges of capacity mu
to implement this price. Integer scaling16 allows prices3/16,4/16,5/16.
The corresponding lower bound is flow_value/16-mu*n_var/2.

Every returned flow is checked for antisymmetry, capacity feasibility,
and conservation. Residual reachability and cut capacity also match the
reported flow value. Feasibility alone suffices for the lower bound;
the cut checks additionally verify the computed optimum. Implication
capacities exceed the sum of finite capacities. No floating optimization
or assumed max-flow accuracy enters the quoted bounds.

## Use half-fullness before relaxing

The first relaxation allows its minimizing assignment to have sparse
early layers and dense late layers. This cannot happen in the actual
family. For consecutive masks of equal weight,

    |H_prev minus H_next| = |H_prev symmetric-difference H_next|/2.

Replace each directed death edge by two opposite edges of half capacity.
This preserves the cost for every admissible half-full history, but gives
a different relaxation outside that family. It removes the artificial
boundary advantage of increasing state weight through a window.

This replacement does not change the encoder or its pairing distribution.
It uses a constraint already present in the restricted proof family.

## Completed results

Both32-update-window replays process all128 windows. For each window they
take the largest verified lower bound among the three fixed prices.

| Death-cost representation | Lower bound on C_T+D | Upper bound on g_T |
| --- | --- | --- |
| Directed deaths | 194833.5 | -26048.5 |
| Symmetric changes, using half-fullness | 216727 | -47942 |

These are deterministic bounds for every half-full history in this
sampled route and matching geometry, for the stated T. They are stronger
evidence than the local search for that fixed set. They do not exclude a
different target set, a different setup, changing state weights, or a
bad word whose existence is not detected by the sufficient rank-gap test.

The first replay took15.60s; the symmetric replay took23.68s. Exhaustive
tests on small windows pass72 cut-value and half-full lower-bound checks,
including both death representations. The gate coefficients do not enter
this geometric computation.

## Larger-window attempt and resource limits

A64-update-window replay at price17/64 processed windows0 through27 before
its25s budget. Their verified lower bounds sum to97239.875. This is only
partial coverage, not an estimate for the remaining windows.

A separate continuation started at window28 and terminated with an
allocation error while forming a sparse residual matrix. That process
exited; no unreported partial contribution from it is used. The final
read-only check showed17.55GiB free disk, below the20GiB guard, and about
8.38GiB free physical /11.45GiB free virtual memory. Those later readings
do not establish the cause of the earlier allocation failure. No files
were deleted and no other process or system setting was changed.

The helper now checks disk space before each flow, catches allocation
failure, and reports only already-verified windows. That reporting-only
change has not been numerically rerun because the disk guard is currently
below threshold. The72 exact tests and the completed32-window replays
precede that change; the cut formulation itself is unchanged.

## Next step

The fixed-set cut bound gives stronger information for this sampled
geometry, but does not remove the union over row sets or exclude rare
setups that affect the earlier moment. Before extending the computation,
identify how to bound a family of target sets or the setup probability
of low cut cost. Merely checking more selected T values is not sufficient.

Another useful refinement is to retain the specified turnover B in the
cut relaxation. If a weighted cut proves C_T+eta*D>=L_eta, then histories
with D=B satisfy C_T+D>=L_eta+(1-eta)*B. This may sharpen the hardest
profile without changing the construction. It still needs the all-target
set and field-count interfaces before becoming a whole-code proof.

Prior restricted certificates and frozen helpers are unchanged. No native
encoding benchmark, construction/library/paper edit, commit, or promotion.
