# One Shared Budget for the Lazy Branch

The class allocation in [BIRTH_DENSITY.md](BIRTH_DENSITY.md) bounds each
outgoing class separately. Those bounds can collectively spend more
mass than the lazy branch has. We can retain the total-mass constraint
when bounding a future cost, instead of constructing another linear
comparison matrix for measures.

Fix a nonnegative function f on inner states. Let v contain upper
bounds on f(0), its supremum on nonzero states, its uniform average on
nonzero states, and its supremum on each class E_l. These are the Z, M,
U, and F_l coordinates, now interpreted as bounds on a future function.
The uniform-average coordinate is a property of f, not an assumption
that the actual entering state is uniform.

For one source coordinate and exact local packet count j, suppose the
weighted nonzero lazy output measure nu satisfies

    total mass <= H,
    nu(b) <= D/L for each b!=0,

where L=2^s-1. Both H and D include alpha=2^(-r). The previously proved
GF fiber bound supplies D; the entering U coordinate averages it over
the exact expansion profiles. For U, the earlier density bound gives
another valid cap on D. The source's existing output-moment bound
supplies H.

If n_l=|E_l|, the lazy contribution to the future cost is at most

    min(H v_M, D v_U, max sum_l x_l v_F_l),

where the maximum ranges over

    0 <= x_l <= D n_l/L,    sum_l x_l <= H.

The maximum is a small linear program. For nonnegative costs, fill the
classes in decreasing order of v_F_l until their caps or the shared
budget are exhausted. All three terms bound the same scalar integral
of f. Taking their minimum is valid; taking componentwise minima of
incompatible measure allocations would not be.

Keep the existing return-to-zero, refresh, and quiet-class terms. Add
the new lazy bound to those terms. Also retain the scalar bound obtained
from each previously valid linear operator, if it is smaller. Applying
this rule to every source coordinate defines a map F_j(v). It is
monotone and positively homogeneous. It need not be linear.

For any future function f bounded by v in the stated coordinates,
F_j(v) bounds the corresponding coordinates after one backward step.
Induction therefore permits composition. `mass_density_potential.py`
implements both outward local evaluation and floating proposals.
Exact small-state tests check three consecutive steps, all local packet
counts, two and three updates, and ordinary, rank-deficient, or zero
feedback. The tests use nonconstant terminal costs.

## Conditional Regional Placement

A region has E steps, each containing W packet positions. Conditional
on k active packets, their locations are uniform among the EW slots.
Let V_(e,k) be an unnormalized backward bound after e steps. Set
V_(0,0)=v. The recursion is

    V_(e,k) = sum_j binom(W,j) F_j(V_(e-1,k-j)).

Terms outside 0<=j<=W or 0<=k-j<=(e-1)W are omitted. Positive
homogeneity absorbs the previous count normalization. Thus
V_(E,k)/binom(EW,k) bounds the future-cost vector conditional on k.
This retains the ordered steps and exact regional count. It does not
replace regional placement by independent packet activity.

For a fixed mean/variance cell, let a_k be the existing checked bound
on the regional count's unnormalized mass, including the activity tilt.
Summing a_k times these conditional maps gives a monotone homogeneous
regional map R. If a strictly positive vector v satisfies R(v)<=rho v,
then induction gives R^B(1)<=rho^B v/min(v), for B regions. Starting at
zero, the moment is therefore at most rho^B v_Z/min(v). Each variance
part may choose its own v. The existing outer-count bound and the
output Chernoff factor must still be applied, and all cells and sparse
occupancies must still be covered.

`regional_potential.py` currently evaluates this recursion and proposes
potentials in floating arithmetic. It checks reconstruction of the saved
linear baseline first. Small tests reproduce all conditional products
of noncommuting linear operators, including values below floating
underflow represented with separate logarithmic scales. Regional
tests also compare the nonlinear bound with the exact two-step,
count-conditioned calculation for a small binary encoder. Regional
outward replay is not yet implemented. These files are not accepted
as complete certificates.

## Initial Evidence

At K=2^20 with three updates, an inner-only iid screen at the saved
output tilts improves the log2 bound by about 7105, 7465, and 4873 bits
at reference activities 0.35, 0.40, and 0.45. The comparison uses the
same local maps and 16384 steps. These gains omit regional count
constraints and outer multiplicities; they are not failure margins.
The first actual-regional calculation improves the activity-0.40 score
from +3594.5604 to +2449.4363 after six potential iterations. At
activity 0.45, it improves +11139.1129 to +10970.7775. Both still fail.
The inner-only gains therefore do not transfer in full. For both
activities, variance part 4 (zero-based) dominates the tested sum; its
average-variance interval is approximately [0.02778,0.03472] at activity
0.40 and [0.03325,0.04156] at 0.45. These are diagnostic values for
one fixed output tilt at each activity, not optimized limits.

Increasing variance resolution from 16 to 64 parts gives +2261.5565
and +10681.9476, respectively. This modest improvement does not justify
implementing an outward regional version of this bound yet.

`composition_probe.py` then fixes nearby integer compositions of the
positive outer envelope and computes their Poisson-binomial count laws
directly. The representative three-update scores remain positive,
including after nonlinear budget optimization. These are not actual
BCH words, and rounding can move the composition outside the narrow
source cell. The diagnostic records those membership checks explicitly.
It is evidence about the comparison bound, not the code's distance.

A subsequent ablation locates a more specific loss: the old weighted
zero-to-zero entry was still a coarse upper bound. It can dominate even
when the diagnostic deletes all returns from nonzero states. The
[exact-zero refinement](ZERO_TRANSITION.md) computes that entry with
the existing Fourier machinery. It is being tested before further work
on the nonlinear regional implementation.
