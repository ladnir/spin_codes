# Balanced shell measures

This evaluator changes the witnesses in the existing adaptive occupation
bound. It preserves its three-state transfer and its one-reused-constituent
model. All reported values use nearest binary64 arithmetic.

Partition the nonzero outer spectrum into bands g. For a band with midpoint
weight fraction m_g, choose a nonnegative scale a and set

    logit(p_g) = a logit(m_g).

The default bands are singleton supported weights. Scale zero assigns
probability one half to every band except the all-one singleton, which
retains its exact Bernoulli(1) law. Each p_g is fixed across all regions.
It is not optimized separately for different transfer entries. Compute

    Gamma_g = max_{w in g} A_w / (choose(B,w) p_g^w (1-p_g)^(B-w)),
    rho_g = Gamma_g^(1/B).

For a random outer, A_w is replaced by a simultaneous integer cap for one
full-rank subspace. The cap event uses the smaller of Markov and exact
variance Chebyshev thresholds from `random_spectrum_variance.py`.
Its failure term is charged once in the final union. A stronger event may
reuse an earlier conditional bound only after entrywise cap containment,
outer parameters, source hashes, and failure budgets have been checked.

Starting from exact region support coefficients R_j, retain the existing
entrywise adaptive recurrence

    K_j^(q+1) = max_g rho_g ((1-p_g) K_j^q + p_g K_(j+1)^q).

For G bands and occupation Q, the resulting first-moment upper bound is

    choose(L,Q) G^Q exp(lambda H) e_Z (K_0^Q)^B 1.

The maximum is a relaxation that dominates each fixed band assignment.
The G^Q factor counts all assignments. This is not a product of expected
random spectra. The evaluator also retains the valid message-count bound
when that bound is smaller.

The implementation removes dominated affine functions from
`rho_g ((1-p_g) x + p_g y)` and builds their upper hull as a function of
`y/x`. Binary searches in log space then evaluate the same maximum without
scanning every shell at every recurrence step. Tests compare this hull
against every original measure, including zero entries, and compare the
complete recurrence with the unpruned implementation.

`run_balanced_occupation_grid.py` records every evaluated integer occupation,
the selected scale and tilt, the spectrum event, and authenticated source
dependencies. A coefficient ceiling below L remains partial coverage. A
coarse witness grid can leave sharp variations between neighboring Q;
refined results remain additional bounds and do not replace old receipts.

The second producer version shares each power-of-two region ladder across
all requested message lengths and batches the terminal matrix powers over
occupations and scales. It uses the same inequality. Tests compare its
region coefficients and final CSV observations with the first producer.
