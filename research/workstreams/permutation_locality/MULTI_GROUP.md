# Extending one-column routing beyond two active groups

This analysis concerns the ideal uniform one-column distribution in
[ONE_COLUMN.md](ONE_COLUMN.md), at K=2^20 with BCH[256,128] and IMT(128,19).
The bad event is output weight at most 209715. The established one- and
two-group bounds remain unchanged. This extension is not yet a full-code
certificate. The subsequent [state-memory refinement](STATE_MEMORY.md)
closes the full three- and four-group support ranges; this note preserves
the diagnostics that motivated it.

## Average the placement without replacing it by independent slots

Fix a message and the set of active groups in a region. There are 64 epochs
with 32 four-bit windows each. A uniform permutation places the active
groups into distinct windows. Let T_j be an entrywise upper transfer matrix
for an epoch containing j active groups, including their lane shuffles and
the IMT update. The positive seven-coordinate representation is the one
used by the existing two-group verifier.

For r active groups in a region, a valid transfer envelope is

    R_r = [x^r] (sum_j binom(32,j) x^j T_j)^64 / binom(2048,r).

Matrix products follow epoch order; they need not commute. Each slot subset
is counted once. Labels can be removed because T_j bounds every ordered
tuple of nonzero column weights. The recurrence in `occupancy_model.py`
checks against exact enumeration of small noncommuting examples and against
the existing R_0, R_1, and R_2 implementation.

For three and four groups in an epoch, `window_feedback.py` supplies exact
zero-feedback counts. A separate atom bound controls cancellation from an
arbitrary incoming state. Fix the other j-1 windows and their lane masks.
Every pair of four-bit window images intersects only at zero. Thus a fixed
nonzero target syndrome is represented in at most one remaining window,
by at most one mask. If the final group's column weight is a, the conditional
probability is at most

    1 / ((32-j+1) binom(4,a)).

Choose the final group with largest binom(4,a). This bound does not assert
uniform feedback and does not prevent larger collections from having zero
feedback. Their exact zero counts are retained separately.

## Keep the support coefficients

For three groups with supports u, v, w, their shared coordinate permutations
give independent uniform support subsets of the 256 regions. Define

    P(x,y,z) = R_0 + R_1(x+y+z) + R_2(xy+xz+yz) + R_3 xyz.

The weighted moment is bounded by the coefficient of x^u y^v z^w in
e_0 P^256 1, divided by binom(256,u) binom(256,v) binom(256,w).
Multiply by exp(lambda*209715) for the Chernoff failure bound.

`occupancy_screen.py` bounds these coefficients by evaluating P at positive
arguments. Its numerical optimizer only chooses witnesses. It does not
certify a global optimum or cover unsampled supports.

`occupancy_exact.py` instead propagates all coefficients in a truncated
lower cube. Nonnegative polynomial degrees ensure that discarded larger
indices cannot affect retained coefficients. Exact toy enumeration checks
the recurrence, and selected coefficients are compared against the Cauchy
bound at the same lambda. The current implementation uses binary64 and is
only a diagnostic; a successful screen still needs outward replay and a
bound for the omitted supports.

## Use cumulative outer counts without inventing a spectrum

Let A(u) count the nonzero ordered four-row tuples whose union support has
size at most u. Let U(u) be a nondecreasing authenticated upper bound on
A(u). For a nonnegative decreasing function f on 0,...,m, summation by parts
gives

    sum_{u=0}^m (A(u)-A(u-1)) f(u)
      <= U(m) f(m) + sum_{u=0}^{m-1} U(u)(f(u)-f(u+1))
       = sum_{u=0}^m (U(u)-U(u-1)) f(u),

where A(-1)=U(-1)=0. All coefficients multiplying U are nonnegative.
This remains valid for a lower support prefix; extend f by zero above m.

Apply this one-dimensional inequality successively to each group. The
group message counts form a product, so a coordinatewise decreasing
multivariate f is sufficient. Mixed differences need not be nonnegative.
The differences of U are a dominating measure for these functions, not
upper bounds on the actual shell counts.

The pointwise minimum over several lambda witnesses remains an upper
bound. `occupancy_cdf.py` then takes suffix maxima along every support
coordinate to obtain a decreasing majorant. It checks the resulting
inequality on 243 integer examples, including a separate example with a
negative mixed difference. The calculation includes binom(2048,3) group
locations and all ordered support vectors; no extra permutation factor is
needed when summing the full cube.

## Diagnostic status

The optimized point screen with the basis-lattice outer refinement gives
log2 contributions -29.0121 at (80,80,80), -34.8582 at (76,76,76), and
-49.1649 at (96,96,96). These include all group locations but concern only
the indicated support vector. They are upper-bound diagnostics, not
observed failure rates.

The first direct-coefficient run used lambda at most .001. That range was
insufficient: the above difficult points prefer .00125 or .0016. Its vacuous
aggregate is therefore not evidence against coefficient extraction. The
screen now reports the selected lambda explicitly, and the coefficient
diagnostic checks the Cauchy comparison at matched lambda values.

The corrected direct run uses lambda in {.00064,.001,.00125,.0016} and
covers every ordered support triple in 38,...,100. Its aggregate log2 bound
is -31.3806 when treating cumulative counts as individual support caps, or
-31.8565 using the decreasing-majorant CDF calculation. The largest latter
terms occur near (83,83,84), at log2 contribution -40.8479 per ordered
triple. A full-code claim cannot use those individual terms in place of
their sum. Supports above 100 and larger occupancies are also still open.

Coefficient extraction improves the conditional bound at (80,80,80) from
log2 -702.0544 to -714.6603. It removes substantial slack but does not close
the aggregate. A tighter transfer envelope or outer count is still needed.

An additional outer-count experiment, `hyperplane_caps.py`, counts every
rank-h subspace through its 2^h-1 hyperplanes. For an extension of support
at most u and a hyperplane of support v, the outside coset has average
weight at most u-v/2. Counting short coset representatives and dividing by
their guaranteed multiplicity gives an extension cap. The CDF inequality
above sums these caps over hyperplanes. Exact checks pass on 87 small-code
inequalities. This refinement does not improve the existing rank-four
caps at the difficult supports 72--96, so it is not the missing gain here.

## Locate the inner-state loss

`occupancy_ablation.py` changes selected transfer entries only as a
counterfactual diagnostic. Those changed matrices are not bounds for the
implemented code. At (83,83,83), the point log2 contributions are:

| Diagnostic transfer | Log2 contribution |
|---|---:|
| Existing envelope, with a finer lambda grid | -33.7766 |
| Remove the three-group zero-feedback term | -33.7766 |
| Remove active cancellation except uniform-refresh cancellation | -214.1059 |
| Replace the lazy outgoing component at each active epoch by uniform refresh | -220.1695 |

The rare zero-feedback triple is not the dominant loss in this calculation.
The sensitivity is to cancellation after activation and the lazy state
component. This does not distinguish actual cancellation probability from
slack in its bound.

There is a concrete reason to retain more state information. For two
independent one-column feedback samples with fixed column weights, exact
convolution counts give maximum nonzero atom at most 1/512. A single sample
can have an atom of 1/32. The largest two-sample zero atom is still 1/32;
zero must be removed before using the stronger nonzero-density bound.
`feedback_convolution.py` checks all ten unordered weight pairs and exact
total mass and zero-coefficient identities.

This 16-fold density improvement is a local count, not an additional
16-fold improvement of the full distance bound. It suggests an envelope
that distinguishes a single activation from a state surviving at least
two feedback additions. Such an envelope must propagate both nonzero mass
and a pointwise density cap, include tilted output weights, and retain the
same-epoch placement constraints. Applying independent-window counts
directly to distinct windows in one epoch would be incorrect.

The next proof experiment is this state-memory refinement. The underlying
encoder and its distribution need not change for that experiment.

Reproduce the screen and coefficient diagnostic from the repository root:

    python -B research/workstreams/permutation_locality/occupancy_model.py
    python -B research/workstreams/permutation_locality/occupancy_cdf.py
    python -B research/workstreams/permutation_locality/occupancy_screen.py --groups 3 --fine --basis-lattice
    python -B research/workstreams/permutation_locality/occupancy_exact.py --maximum 100

No encoder timings change in this analysis. The measured 1.88--1.90x
speedup and the remaining 2x target are recorded in ONE_COLUMN.md.
